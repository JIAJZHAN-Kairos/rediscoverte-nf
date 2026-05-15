process BUILD_SALMON_INDEX {
    tag "REdiscoverTE"
    label 'salmon'

    input:
    path transcriptome_fasta

    output:
    path "salmon_index", emit: index

    script:
    """
    set -euo pipefail

    fasta_in="${transcriptome_fasta}"
    case "\$fasta_in" in
        *.xz)
            xz -dc "\$fasta_in" > transcriptome.fa
            fasta=transcriptome.fa
            ;;
        *.gz)
            gzip -dc "\$fasta_in" > transcriptome.fa
            fasta=transcriptome.fa
            ;;
        *)
            fasta="\$fasta_in"
            ;;
    esac

    salmon index \\
        -t "\$fasta" \\
        -i salmon_index \\
        -k ${params.index_kmer} \\
        --keepDuplicates \\
        -p ${task.cpus}
    """
}
