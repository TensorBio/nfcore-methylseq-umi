process PICARD_UMIMARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.bam")             , emit: bam
    tuple val(meta), path("*.metrics.txt")     , emit: metrics
    tuple val(meta), path("*.umi_metrics.txt") , emit: umi_metrics
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard UmiAwareMarkDuplicatesWithMateCigar] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    if ("$reads" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    picard \\
        -Xmx${avail_mem}M \\
        UmiAwareMarkDuplicatesWithMateCigar \\
        $args \\
        --INPUT $reads \\
        --OUTPUT ${prefix}.bam \\
        $reference \\
        --METRICS_FILE ${prefix}.UmiAwareMarkDuplicates.metrics.txt \\
        --UMI_METRICS_FILE ${prefix}.UmiAwareMarkDuplicates.umi_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard UmiAwareMarkDuplicatesWithMateCigar --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$reads" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.bam
    touch ${prefix}.UmiAwareMarkDuplicates.metrics.txt
    touch ${prefix}.UmiAwareMarkDuplicates.umi_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard UmiAwareMarkDuplicatesWithMateCigar --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
