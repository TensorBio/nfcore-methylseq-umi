include { FASTQ_ALIGN_BWA               } from '../fastq_align_bwa/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_FIXMATEINFORMATION     } from '../../../modules/nf-core/picard/fixmateinformation/main'
include { PICARD_MARKDUPLICATES         } from '../../../modules/nf-core/picard/markduplicates/main'
include { PICARD_UMIMARKDUPLICATES      } from '../../../modules/nf-core/picard/umimarkduplicates/main'
include { SAMTOOLS_FIXMATE              } from '../../../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_INDEX                } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_UMITAGTRANSFER       } from '../../../modules/local/samtools/umitagtransfer/main'

workflow FASTQ_ALIGN_DEDUP_BWAMEM {

    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_fasta             // channel: [ val(meta), [ fasta ] ]
    ch_fasta_index       // channel: [ val(meta), [ fasta index ] ]
    ch_bwamem_index      // channel: [ val(meta), [ bwamem index ] ]
    skip_deduplication   // boolean: whether to deduplicate alignments
    umi                  // boolean: whether to use UMI-aware deduplication
    fixmate              // string:  tool to add MC tags before UMI deduplication: 'samtools', 'picard', or null (skip)

    main:
    ch_alignment       = channel.empty()
    ch_alignment_index = channel.empty()
    ch_flagstat        = channel.empty()
    ch_stats           = channel.empty()
    ch_idxstats        = channel.empty()
    ch_picard_metrics  = channel.empty()
    ch_multiqc_files   = channel.empty()
    ch_versions        = channel.empty()

    FASTQ_ALIGN_BWA (
        ch_reads,
        ch_bwamem_index,
        true, // val_sort_bam hardcoded to true
        ch_fasta
    )
    ch_alignment        = FASTQ_ALIGN_BWA.out.bam         // channel: [ val(meta), [ bam ] ]
    ch_alignment_index  = FASTQ_ALIGN_BWA.out.bai         // channel: [ val(meta), [ bai ] ]
    ch_stats            = FASTQ_ALIGN_BWA.out.stats       // channel: [ val(meta), path(stats) ]
    ch_flagstat         = FASTQ_ALIGN_BWA.out.flagstat    // channel: [ val(meta), path(flagstat) ]
    ch_idxstats         = FASTQ_ALIGN_BWA.out.idxstats    // channel: [ val(meta), path(idxstats) ]
    ch_versions         = ch_versions.mix(FASTQ_ALIGN_BWA.out.versions.first())

    if (!skip_deduplication) {
        /*
         * Run Picard AddOrReplaceReadGroups to add read group (RG) to reads in bam file
         */
        PICARD_ADDORREPLACEREADGROUPS (
            ch_alignment,
            ch_fasta,
            ch_fasta_index
        )
        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.first())

        if (umi) {
            /*
             * Optionally add MC (mate CIGAR) tags required by UmiAwareMarkDuplicatesWithMateCigar.
             * bwa mem does not emit MC tags; without this step Picard silently skips all pairs
             * (SKIP_PAIRS_WITH_NO_MATE_CIGAR=true default), making UMI deduplication a no-op.
             * Use --fixmate samtools  →  samtools sort -n | fixmate -m | samtools sort
             * Use --fixmate picard    →  picard FixMateInformation (sorts internally)
             */
            if (fixmate == 'samtools') {
                SAMTOOLS_FIXMATE (
                    PICARD_ADDORREPLACEREADGROUPS.out.bam
                )
                ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions_samtools.first())
            } else if (fixmate == 'picard') {
                PICARD_FIXMATEINFORMATION (
                    PICARD_ADDORREPLACEREADGROUPS.out.bam,
                    ch_fasta,
                    ch_fasta_index
                )
                ch_versions = ch_versions.mix(PICARD_FIXMATEINFORMATION.out.versions.first())
            }

            /*
             * Transfer UMI from read name to RX BAM tag before Picard
             */
            SAMTOOLS_UMITAGTRANSFER (
                fixmate == 'samtools' ? SAMTOOLS_FIXMATE.out.bam :
                fixmate == 'picard'   ? PICARD_FIXMATEINFORMATION.out.bam :
                                        PICARD_ADDORREPLACEREADGROUPS.out.bam
            )
            ch_versions = ch_versions.mix(SAMTOOLS_UMITAGTRANSFER.out.versions.first())

            /*
             * Run UMI-aware Picard deduplication
             */
            PICARD_UMIMARKDUPLICATES (
                SAMTOOLS_UMITAGTRANSFER.out.bam,
                ch_fasta,
                ch_fasta_index
            )
            ch_versions = ch_versions.mix(PICARD_UMIMARKDUPLICATES.out.versions.first())

            SAMTOOLS_INDEX (
                PICARD_UMIMARKDUPLICATES.out.bam
            )
            ch_alignment       = PICARD_UMIMARKDUPLICATES.out.bam
            ch_alignment_index = SAMTOOLS_INDEX.out.bai
            ch_picard_metrics  = PICARD_UMIMARKDUPLICATES.out.metrics
            ch_versions        = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
        } else {
            /*
             * Run Picard MarkDuplicates to mark duplicates
             */
            PICARD_MARKDUPLICATES (
                PICARD_ADDORREPLACEREADGROUPS.out.bam,
                ch_fasta,
                ch_fasta_index
            )
            ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

            /*
             * Run samtools index on deduplicated alignment
             */
            SAMTOOLS_INDEX (
                PICARD_MARKDUPLICATES.out.bam
            )
            ch_alignment       = PICARD_MARKDUPLICATES.out.bam
            ch_alignment_index = SAMTOOLS_INDEX.out.bai
            ch_picard_metrics  = PICARD_MARKDUPLICATES.out.metrics
            ch_versions        = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
        }
    }

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_picard_metrics.collect{ _meta, metrics -> metrics }
                        .mix(ch_flagstat.collect{ _meta, flagstat -> flagstat })
                        .mix(ch_stats.collect{ _meta, stats -> stats  })
                        .mix(ch_idxstats.collect{ _meta, stats -> stats  })

    emit:
    bam               = ch_alignment                     // channel: [ val(meta), [ bam ]       ]
    bai               = ch_alignment_index               // channel: [ val(meta), [ bai ]       ]
    samtools_flagstat = ch_flagstat                      // channel: [ val(meta), [ flagstat ]  ]
    samtools_stats    = ch_stats                         // channel: [ val(meta), [ stats ]     ]
    samtools_idxstats = ch_idxstats                      // channel: [ val(meta), [ idxstats ]  ]
    picard_metrics    = ch_picard_metrics                // channel: [ val(meta), [ metrics ]   ]
    multiqc           = ch_multiqc_files                 // channel: [ *{html,txt}              ]
    versions          = ch_versions                      // channel: [ versions.yml             ]
}
