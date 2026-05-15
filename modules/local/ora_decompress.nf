process ORA_DECOMPRESS {
    tag "${meta.id}"
    label 'ora'

    input:
    tuple val(meta), path(ora_files)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq

    script:
    """
    set -euo pipefail
    for f in ${ora_files.join(' ')}; do
        case "\$f" in
            *.ora)
                orad --gz -t ${task.cpus} -P . "\$f"
                ;;
            *.fastq.gz|*.fq.gz)
                # passthrough: link (avoid copying)
                ln -sf "\$(readlink -f \$f)" "\$(basename \$f)"
                ;;
            *)
                echo "Unsupported input extension: \$f" >&2
                exit 1
                ;;
        esac
    done
    """
}
