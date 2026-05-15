#!/usr/bin/env nextflow

/*
 * REdiscoverTE-nf
 * Minimal pipeline: (ORA decompress) -> Salmon quant -> rollup.R
 * Input:  samplesheet.csv with .ora or .fastq.gz reads
 * Output: per-sample quant.sf + rolled-up RE expression tables
 */

nextflow.enable.dsl = 2

include { ORA_DECOMPRESS } from './modules/local/ora_decompress'
include { BUILD_SALMON_INDEX } from './modules/local/build_salmon_index'
include { SALMON_QUANT   } from './modules/local/salmon_quant'
include { ROLLUP         } from './modules/local/rollup'

// ---------- Workflow ----------
workflow {
    // Parameter checks
    def required = ['input', 'outdir', 'rollup_annotation']
    required.each { p ->
        if (!params[p]) {
            exit 1, "Missing required parameter: --${p}"
        }
    }
    if (!params.salmon_index && !params.transcriptome_fasta) {
        exit 1, "Provide either --salmon_index or --transcriptome_fasta"
    }

    // Parse samplesheet: sample,fastq_1,fastq_2.
    // Multiple rows with the same sample are treated as separate lanes/read groups
    // and are quantified together by Salmon.
    ch_samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def f1 = file(row.fastq_1, checkIfExists: true)
            def f2 = (row.fastq_2 && row.fastq_2.trim()) ? file(row.fastq_2, checkIfExists: true) : null
            meta.single_end = (f2 == null)
            tuple(row.sample, meta.single_end, f2 ? [[f1, f2]] : [[f1]])
        }
        .groupTuple()
        .map { sample, single_end_values, reads_nested ->
            def single_end_set = single_end_values.unique()
            if (single_end_set.size() != 1) {
                exit 1, "Sample ${sample} mixes single-end and paired-end rows"
            }
            def meta = [id: sample, single_end: single_end_set[0]]
            def reads = reads_nested.flatten()
            tuple(meta, reads)
        }

    // Split: .ora needs decompression, .fastq.gz passes through
    ch_branched = ch_samples.branch { meta, reads ->
        ora:   reads.any { it.name.endsWith('.ora') }
        ready: true
    }

    ch_decompressed = ORA_DECOMPRESS(ch_branched.ora).fastq

    ch_quant_in = ch_decompressed.mix(ch_branched.ready)

    // Salmon quant against a prebuilt index, or build one from FASTA for this run
    if (params.salmon_index) {
        ch_index = Channel.value(file(params.salmon_index, checkIfExists: true))
    } else {
        ch_fasta = Channel.value(file(params.transcriptome_fasta, checkIfExists: true))
        BUILD_SALMON_INDEX(ch_fasta)
        ch_index = BUILD_SALMON_INDEX.out.index
    }
    SALMON_QUANT(ch_quant_in, ch_index)

    // Rollup all quant.sf via rollup.R
    ch_anno = Channel.value(file(params.rollup_annotation, checkIfExists: true))
    ROLLUP(SALMON_QUANT.out.quant_dir.collect(), ch_anno)
}
