# rediscoverte-nf

Minimal Nextflow DSL2 pipeline for [REdiscoverTE](https://github.com/ucsffrancislab/REdiscoverTE)
(transposable / repetitive element quantification from RNA-seq).

```
.ora or .fastq.gz  ->  [ORA_DECOMPRESS]  ->  [SALMON_QUANT]  ->  [ROLLUP]  ->  RE expression RDS
```

Three processes, no overhead. Designed for Seqera Platform on AWS Batch but
runs anywhere Nextflow + Docker runs.

## What it does

1. **ORA_DECOMPRESS** (optional, per-sample): if a fastq input is `.ora`,
   decompress with Illumina `orad` to `.fastq.gz`. `.fastq.gz` inputs are
   passed through.
2. **SALMON_QUANT** (per-sample): quasi-mapping + EM with the prebuilt
   REdiscoverTE Salmon index (GENCODE basic transcripts + ~5M RE copies).
3. **ROLLUP** (gather): runs Genentech's `rollup.R` to aggregate per-RE-copy
   counts to `repName` / `repFamily` / `repClass` levels, split by genomic
   context (intergenic / intronic / exonic).

## Prerequisites

Build / prepare these **once** before running:

### A. Salmon index (run on EC2 / HPC, ~30 min, ~90 GB output)

```bash
./scripts/build_index.sh \
    /path/to/REdiscoverTE_whole_transcriptome_hg38.fa.xz \
    s3://YOUR_BUCKET/refs/REdiscoverTE/salmon_index/
```

The FASTA comes from `REdiscoverTE_1.0.1.tar.gz` (Genentech):
`rollup_annotation/REdiscoverTE_whole_transcriptome_hg38.fa.xz`

The script enforces `--keepDuplicates` (critical — without it Salmon collapses
near-identical TE copies and destroys repeat-level signal).

### B. rollup annotation

Upload the `rollup_annotation/` directory (from same Genentech tarball) to S3:

```bash
aws s3 cp --recursive rollup_annotation/ \
    s3://YOUR_BUCKET/refs/REdiscoverTE/rollup_annotation/ \
    --exclude "*.fa*"
```

### C. Containers

Two custom images plus one off-the-shelf Salmon container.

```bash
# orad (Illumina DRAGEN ORA decompression v2.7.0; ~1.5 GB image)
cd containers/orad
docker build -t ghcr.io/jiajzhan-kairos/rediscoverte-orad:2.7.0 .
docker push    ghcr.io/jiajzhan-kairos/rediscoverte-orad:2.7.0

# rollup (R + edgeR + rollup.R)
cd ../rollup
docker build -t ghcr.io/jiajzhan-kairos/rediscoverte-rollup:1.0 .
docker push    ghcr.io/jiajzhan-kairos/rediscoverte-rollup:1.0
```

If you push images under a different GHCR namespace, override
`--orad_container` and `--rollup_container` at launch.

Salmon container is bioconda's: `quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2` — no build needed.

## Samplesheet

`samplesheet.csv` columns:

| sample | fastq_1 | fastq_2 |
|---|---|---|
| sample id (no spaces) | path/URI to R1 (`.ora` or `.fastq.gz`) | path/URI to R2 (same; empty for SE) |

See `assets/samplesheet.example.csv`.

Use multiple rows with the same `sample` for multi-lane or rerun libraries;
the pipeline quantifies all rows for that sample together in one Salmon call.

## Run

### Locally (smoke test with Docker)

```bash
nextflow run main.nf \
    --input          assets/samplesheet.example.csv \
    --outdir         results \
    --salmon_index   /local/path/to/REdiscoverTE_salmon_index \
    --rollup_annotation /local/path/to/rollup_annotation \
    -profile docker
```

### Seqera Platform (AWS Batch)

1. Push this repo to GitHub.
2. Seqera → Add Pipeline → paste GitHub URL → revision `main`.
3. At Launch, set parameters:

```yaml
input:              s3://YOUR_BUCKET/samplesheets/run1.csv
outdir:             s3://YOUR_BUCKET/results/run1
salmon_index:       s3://YOUR_BUCKET/refs/REdiscoverTE/salmon_index
rollup_annotation:  s3://YOUR_BUCKET/refs/REdiscoverTE/rollup_annotation
```

4. Pick the `seqera` profile. The compute environment should provide the AWS
   Batch queue and work directory. This profile uses regular AWS Batch Docker
   execution and does not require Fusion.

For `.fastq.ora` input, use the same samplesheet schema as `.fastq.gz`; the
pipeline automatically routes `.ora` reads through `ORA_DECOMPRESS` before
Salmon.

## Outputs

```
results/
├── salmon/
│   └── <sample>/
│       ├── quant.sf            # transcript-level counts (incl. all RE copies)
│       ├── lib_format_counts.json
│       └── aux_info/...
├── rollup_out/
│   ├── GENE_x_counts.RDS       # gene-level (GENCODE basic)
│   ├── RE_all_*.RDS            # all RE
│   ├── RE_intron_*.RDS         # intronic RE
│   ├── RE_exon_*.RDS           # exonic RE
│   ├── RE_intergenic_*.RDS     # intergenic RE
│   └── rmsk_annotation.RDS
├── rollup_metadata.tsv
└── pipeline_info/              # nextflow report, timeline, trace
```

The `RE_*` RDS files feed directly into edgeR for differential expression
(ecDNA vs non-ecDNA, etc.).

## Sanity checks after first run

```bash
# 1. quant.sf should have ~5M rows (GENCODE basic ~95k + RE copies ~5M)
wc -l results/salmon/<sample>/quant.sf

# 2. TE entries should be present (look for L1, Alu, LTR, HERV)
grep -E "L1MA|AluY|LTR|HERV" results/salmon/<sample>/quant.sf | head
```

If `quant.sf` has only ~100k rows, the Salmon index was built without
`--keepDuplicates` — rebuild with the provided script.

## Tuning notes

- `salmon_libtype` defaults to `A` (auto). Override to `ISR` / `IU` etc. via `--salmon_libtype` if needed.
- `rollup_nozero=true` (default) drops zero rows to speed up downstream.
- Salmon `--seqBias --gcBias --validateMappings` are hardcoded to match the Genentech REdiscoverTE protocol.

## License

Pipeline code: MIT.
REdiscoverTE / rollup.R: MIT (Genentech).
Illumina DRAGEN ORA Decompression binary inside the orad container: see Illumina EULA.
