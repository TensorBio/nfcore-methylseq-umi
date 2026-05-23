process PICARD_FIXMATEINFORMATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data' :
        'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard FixMateInformation] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    if ("$reads" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    picard \\
        -Xmx${avail_mem}M \\
        FixMateInformation \\
        $args \\
        --INPUT $reads \\
        --OUTPUT ${prefix}.bam \\
        $reference \\
        --SORT_ORDER coordinate

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard FixMateInformation --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$reads" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard FixMateInformation --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
