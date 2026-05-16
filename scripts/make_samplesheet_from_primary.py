#!/usr/bin/env python3
"""
Create a REdiscoverTE samplesheet from mounted Cancer Atlas primary WTS ORA data.

The script scans run/hash/WTS_NebRNA folders, pairs R1/R2 .fastq.ora files, and
writes a Nextflow samplesheet with columns: sample,fastq_1,fastq_2.
Multiple rows with the same sample are allowed and are quantified together by
the pipeline.
"""

import argparse
import csv
import re
import sys
from pathlib import Path


READ_RE = re.compile(r"^(?P<sample>.+)_R(?P<read>[12])_001\.fastq\.ora$")
LANE_SUFFIX_RE = re.compile(r"_S\d+_L\d{3}$")


def sample_from_pair_key(pair_key: str) -> str:
    sample = LANE_SUFFIX_RE.sub("", pair_key)
    if sample.endswith("_rerun"):
        sample = sample.removesuffix("_rerun")
    return sample


def mounted_path_to_s3(path: Path, mount_root: Path, bucket: str) -> str:
    rel = path.resolve().relative_to(mount_root.resolve())
    return f"s3://{bucket}/{rel.as_posix()}"


def discover_pairs(primary_dir: Path, mount_root: Path, bucket: str):
    records = {}
    ignored = []

    for fastq in sorted(primary_dir.glob("*/*/WTS_NebRNA/*.fastq.ora")):
        match = READ_RE.match(fastq.name)
        if not match:
            ignored.append((str(fastq), "filename_does_not_match_expected_R1_R2_pattern"))
            continue

        pair_key = match.group("sample")
        sample = sample_from_pair_key(pair_key)
        read = match.group("read")
        key = (sample, pair_key)
        records.setdefault(key, {"sample": sample, "fastq_1": None, "fastq_2": None})
        column = f"fastq_{read}"

        if records[key][column] is not None:
            ignored.append((str(fastq), f"duplicate_{column}_for_pair_key:{pair_key}"))
            continue

        records[key][column] = mounted_path_to_s3(fastq, mount_root, bucket)

    complete = []
    problems = []
    for (_sample, pair_key), rec in sorted(records.items(), key=lambda item: (item[1]["sample"], item[0][1])):
        if rec["fastq_1"] and rec["fastq_2"]:
            complete.append(rec)
        else:
            problems.append({
                "sample": rec["sample"],
                "pair_key": pair_key,
                "fastq_1": rec["fastq_1"] or "",
                "fastq_2": rec["fastq_2"] or "",
                "problem": "missing_R1_or_R2",
            })

    return complete, problems, ignored


def write_samplesheet(records, output_csv: Path):
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "fastq_1", "fastq_2"])
        writer.writeheader()
        writer.writerows(records)


def write_qc(problems, ignored, output_tsv: Path):
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["kind", "sample", "pair_key", "path", "fastq_1", "fastq_2", "problem"], delimiter="\t")
        writer.writeheader()
        for problem in problems:
            writer.writerow({
                "kind": "pair_problem",
                "sample": problem["sample"],
                "pair_key": problem["pair_key"],
                "path": "",
                "fastq_1": problem["fastq_1"],
                "fastq_2": problem["fastq_2"],
                "problem": problem["problem"],
            })
        for path, reason in ignored:
            writer.writerow({
                "kind": "ignored_file",
                "sample": "",
                "pair_key": "",
                "path": path,
                "fastq_1": "",
                "fastq_2": "",
                "problem": reason,
            })


def parse_args():
    parser = argparse.ArgumentParser(description="Build REdiscoverTE samplesheet from mounted primary WTS ORA files")
    parser.add_argument(
        "--primary-dir",
        type=Path,
        default=Path("/workspace/data/research-data-550435500918-ap-southeast-2/byob-icav2/project-cancer-atlas-research/primary"),
        help="Mounted primary data directory to scan",
    )
    parser.add_argument(
        "--mount-root",
        type=Path,
        default=Path("/workspace/data/research-data-550435500918-ap-southeast-2"),
        help="Mounted root corresponding to the S3 input bucket",
    )
    parser.add_argument(
        "--bucket",
        default="research-data-550435500918-ap-southeast-2",
        help="S3 bucket name corresponding to --mount-root",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/workspace/data/grimmond-research-nextflow-98050479638/refs/REdiscoverTE/samplesheets/cancer_atlas_wts_ora_samplesheet.csv"),
        help="Output samplesheet CSV path",
    )
    parser.add_argument(
        "--qc-output",
        type=Path,
        default=Path("/workspace/data/grimmond-research-nextflow-98050479638/refs/REdiscoverTE/samplesheets/cancer_atlas_wts_ora_samplesheet_qc.tsv"),
        help="Output QC TSV path",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.primary_dir.exists():
        sys.exit(f"Primary directory not found: {args.primary_dir}")
    if not args.mount_root.exists():
        sys.exit(f"Mount root not found: {args.mount_root}")

    records, problems, ignored = discover_pairs(args.primary_dir, args.mount_root, args.bucket)
    if not records:
        sys.exit(f"No complete R1/R2 .fastq.ora pairs found under: {args.primary_dir}")

    write_samplesheet(records, args.output)
    write_qc(problems, ignored, args.qc_output)

    samples = sorted({record["sample"] for record in records})
    print(f"Wrote samplesheet: {args.output}")
    print(f"Wrote QC report:   {args.qc_output}")
    print(f"Complete pairs:    {len(records)}")
    print(f"Unique samples:    {len(samples)}")
    print(f"Pair problems:     {len(problems)}")
    print(f"Ignored files:     {len(ignored)}")


if __name__ == "__main__":
    main()
