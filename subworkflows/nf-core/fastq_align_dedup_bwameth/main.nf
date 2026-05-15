include { BWAMETH_ALIGN                                 } from '../../../modules/nf-core/bwameth/align/main'
include { PARABRICKS_FQ2BAMMETH                         } from '../../../modules/nf-core/parabricks/fq2bammeth/main'
include { SAMTOOLS_SORT                                 } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGNMENTS   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT                             } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_STATS                                } from '../../../modules/nf-core/samtools/stats/main'
include { PICARD_MARKDUPLICATES                         } from '../../../modules/nf-core/picard/markduplicates/main'
include { PICARD_UMIMARKDUPLICATES                      } from '../../../modules/nf-core/picard/umimarkduplicates/main'
include { SAMTOOLS_UMITAGTRANSFER                       } from '../../../modules/local/samtools/umitagtransfer/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUPLICATED } from '../../../modules/nf-core/samtools/index/main'

workflow FASTQ_ALIGN_DEDUP_BWAMETH {

    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_fasta             // channel: [ val(meta), [ fasta ] ]
    ch_fasta_index       // channel: [ val(meta), [ fasta index ] ]
    ch_bwameth_index     // channel: [ val(meta), [ bwameth index ] ]
    skip_deduplication   // boolean: whether to deduplicate alignments
    use_gpu              // boolean: whether to use GPU or CPU for bwameth alignment
    umi                  // boolean: whether to use UMI-aware deduplication

    main:
    ch_alignment         = channel.empty()
    ch_alignment_index   = channel.empty()
    ch_samtools_flagstat = channel.empty()
    ch_samtools_stats    = channel.empty()
    ch_picard_metrics    = channel.empty()
    ch_multiqc_files     = channel.empty()
    ch_versions          = channel.empty()

      /*
        * Align with parabricks GPU enabled fq2bammeth implementation of bwameth
        */
        PARABRICKS_FQ2BAMMETH (
            ch_reads,
            ch_fasta,
            ch_bwameth_index,
            [] // known sites
        )
        ch_alignment = PARABRICKS_FQ2BAMMETH.out.bam
        ch_versions  = ch_versions.mix(PARABRICKS_FQ2BAMMETH.out.versions)

    /*
     * Sort raw output BAM
     */
    SAMTOOLS_SORT (
        ch_alignment,
        [[:],[]], // [ [meta], [fasta]]
        ''
    )
    ch_alignment = SAMTOOLS_SORT.out.bam

    /*
     * Run samtools index on alignment
     */
    SAMTOOLS_INDEX_ALIGNMENTS (
        ch_alignment
    )
    ch_alignment_index = SAMTOOLS_INDEX_ALIGNMENTS.out.bai
    ch_versions        = ch_versions.mix(SAMTOOLS_INDEX_ALIGNMENTS.out.versions)

    /*
     * Run samtools flagstat
     */
    SAMTOOLS_FLAGSTAT (
        ch_alignment.join(ch_alignment_index)
    )
    ch_samtools_flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
    ch_versions          = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    /*
     * Run samtools stats
     */
    SAMTOOLS_STATS (
        ch_alignment.join(ch_alignment_index),
        [[:],[]] // [ [meta], [fasta]]
    )
    ch_samtools_stats = SAMTOOLS_STATS.out.stats

    if (!skip_deduplication) {
        if (umi) {
            /*
            * Transfer UMI from read name (last colon-delimited field) to RX BAM tag
            * so that Picard UmiAwareMarkDuplicatesWithMateCigar can use it.
            * Read name format: readname:UMI1+UMI2  →  RX:Z:UMI1+UMI2
            */
            SAMTOOLS_UMITAGTRANSFER (
                ch_alignment
            )
            ch_versions = ch_versions.mix(SAMTOOLS_UMITAGTRANSFER.out.versions)

            /*
            * Run Picard UmiAwareMarkDuplicatesWithMateCigar (UMI-aware deduplication)
            */
            PICARD_UMIMARKDUPLICATES (
                SAMTOOLS_UMITAGTRANSFER.out.bam,
                ch_fasta,
                ch_fasta_index
            )
            /*
             * Run samtools index on deduplicated alignment
            */
            SAMTOOLS_INDEX_DEDUPLICATED (
                PICARD_UMIMARKDUPLICATES.out.bam
            )
            ch_alignment       = PICARD_UMIMARKDUPLICATES.out.bam
            ch_alignment_index = SAMTOOLS_INDEX_DEDUPLICATED.out.bai
            ch_picard_metrics  = PICARD_UMIMARKDUPLICATES.out.metrics
            ch_versions        = ch_versions.mix(PICARD_UMIMARKDUPLICATES.out.versions)
            ch_versions        = ch_versions.mix(SAMTOOLS_INDEX_DEDUPLICATED.out.versions)
        } else {
            /*
            * Run Picard MarkDuplicates
            */
            PICARD_MARKDUPLICATES (
                ch_alignment,
                ch_fasta,
                ch_fasta_index
            )
            /*
             * Run samtools index on deduplicated alignment
            */
            SAMTOOLS_INDEX_DEDUPLICATED (
                PICARD_MARKDUPLICATES.out.bam
            )
            ch_alignment       = PICARD_MARKDUPLICATES.out.bam
            ch_alignment_index = SAMTOOLS_INDEX_DEDUPLICATED.out.bai
            ch_picard_metrics  = PICARD_MARKDUPLICATES.out.metrics
            ch_versions        = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
            ch_versions        = ch_versions.mix(SAMTOOLS_INDEX_DEDUPLICATED.out.versions)
        }
    }

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_picard_metrics.collect{ _meta, metrics -> metrics }
                        .mix(ch_samtools_flagstat.collect{ _meta, flagstat -> flagstat })
                        .mix(ch_samtools_stats.collect{ _meta, stats -> stats  })


    emit:
    bam               = ch_alignment                     // channel: [ val(meta), [ bam ]       ]
    bai               = ch_alignment_index               // channel: [ val(meta), [ bai ]       ]
    samtools_flagstat = ch_samtools_flagstat             // channel: [ val(meta), [ flagstat ]  ]
    samtools_stats    = ch_samtools_stats                // channel: [ val(meta), [ stats ]     ]
    picard_metrics    = ch_picard_metrics                // channel: [ val(meta), [ metrics ]   ]
    multiqc           = ch_multiqc_files                 // channel: [ *{html,txt}              ]
    versions          = ch_versions                      // channel: [ versions.yml             ]
}
