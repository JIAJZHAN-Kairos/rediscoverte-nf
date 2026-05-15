#!/usr/bin/env nextflow

/*
 * REdiscoverTE-nf
 * Minimal pipeline: (ORA decompress) -> Salmon quant -> rollup.R
 * Input:  samplesheet.csv with .ora or .fastq.gz reads
 * Output: per-sample quant.sf + rolled-up RE expression tables
 */

nextflow.enable.dsl = 2

include { ORA_DECOMPRESS } from './modules/local/ora_decompress'
include { SALMON_QUANT   } from './modules/local/salmon_quant'
include { ROLLUP         } from './modules/local/rollup'

// ---------- Parameter checks ----------
def required = ['input', 'outdir', 'salmon_index', 'rollup_annotation']
required.each { p ->
    if (!params[p]) {
        exit 1, "Missing required parameter: --${p}"
    }
}

// ---------- Workflow ----------
workflow {

    // Parse samplesheet: sample,fastq_1,fastq_2
    ch_samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def f1 = file(row.fastq_1, checkIfExists: true)
            def f2 = (row.fastq_2 && row.fastq_2.trim()) ? file(row.fastq_2, checkIfExists: true) : null
            meta.single_end = (f2 == null)
            def reads = f2 ? [f1, f2] : [f1]
            tuple(meta, reads)
        }

    // Split: .ora needs decompression, .fastq.gz passes through
    ch_branched = ch_samples.branch { meta, reads ->
        ora:   reads.any { it.name.endsWith('.ora') }
        ready: true
    }

    ch_decompressed = ORA_DECOMPRESS(ch_branched.ora).fastq

    ch_quant_in = ch_decompressed.mix(ch_branched.ready)

    // Salmon quant against prebuilt index
    ch_index = Channel.value(file(params.salmon_index, checkIfExists: true))
    SALMON_QUANT(ch_quant_in, ch_index)

    // Rollup all quant.sf via rollup.R
    ch_anno = Channel.value(file(params.rollup_annotation, checkIfExists: true))
    ROLLUP(SALMON_QUANT.out.quant_dir.collect(), ch_anno)
}

workflow.onComplete {
    log.info "Pipeline ${workflow.success ? 'succeeded' : 'FAILED'}: ${workflow.duration}"
}
