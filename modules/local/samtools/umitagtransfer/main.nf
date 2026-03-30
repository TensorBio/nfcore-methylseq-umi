process SAMTOOLS_UMITAGTRANSFER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args   ?: ''
    def rx_tag  = task.ext.rx_tag ?: 'RX'
    """
    samtools view -h ${args} ${bam} \\
        | awk -v tag=${rx_tag} '
            /^@/ { print; next }
            {
                # UMI is the last colon-delimited field of the read name (QNAME)
                # e.g. readname:GTTCGGA+ACGTATC -> RX:Z:GTTCGGA+ACGTATC
                n = split(\$1, a, ":")
                umi = a[n]
                printf "%s\\t%s:Z:%s\\n", \$0, tag, umi
            }
        ' \\
        | samtools view -bS -o ${prefix}.umi_tagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_tagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
