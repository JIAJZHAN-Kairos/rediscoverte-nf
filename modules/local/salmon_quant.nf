process SALMON_QUANT {
    tag "${meta.id}"
    label 'salmon'
    publishDir "${params.outdir}/salmon", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}"),                 emit: quant_dir
    tuple val(meta), path("${meta.id}/quant.sf"),        emit: quant_sf
    tuple val(meta), path("${meta.id}/lib_format_counts.json"), emit: libstats, optional: true

    script:
    def reads_arg = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    salmon quant \\
        --seqBias \\
        --gcBias \\
        --validateMappings \\
        -i ${index} \\
        -l ${params.salmon_libtype} \\
        ${reads_arg} \\
        -o ${meta.id} \\
        -p ${task.cpus}
    """
}
