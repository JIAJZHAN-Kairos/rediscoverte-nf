process ROLLUP {
    tag "rollup"
    label 'rollup'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path quant_dirs, stageAs: 'salmon/*'
    path rollup_annotation

    output:
    path "rollup_out",        emit: rollup
    path "rollup_metadata.tsv", emit: metadata

    script:
    def nozero = params.rollup_nozero ? '--nozero' : ''
    """
    set -euo pipefail

    # Build metadata TSV consumed by rollup.R
    echo -e "sample\\tquant_sf_path" > rollup_metadata.tsv
    for d in salmon/*/; do
        name=\$(basename \$d)
        echo -e "\${name}\\t\${d}quant.sf" >> rollup_metadata.tsv
    done

    cat rollup_metadata.tsv

    Rscript /opt/rollup/rollup.R \\
        --metadata=rollup_metadata.tsv \\
        --datadir=${rollup_annotation} \\
        ${nozero} \\
        --threads=${task.cpus} \\
        --assembly=${params.assembly} \\
        --outdir=rollup_out
    """
}
