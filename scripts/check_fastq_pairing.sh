#!/usr/bin/env bash
set -euo pipefail

# Read-only FASTQ validation. Raw FASTQ files are never modified.
# Default: check all supplied libraries.
# Usage:
#   bash scripts/check_fastq_pairing.sh
#   bash scripts/check_fastq_pairing.sh all
#   bash scripts/check_fastq_pairing.sh groupB_rep1

SCRIPT_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
PROJECT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
INPUT_DIR="$PROJECT_DIR/data"
OUTPUT_BASE="$PROJECT_DIR/results/pairing_check"
TARGET="all"

if [[ $# -ge 1 ]]; then
    TARGET="$1"
fi

python3 - "$INPUT_DIR" "$OUTPUT_BASE" "$TARGET" <<'PY'
import csv
import gzip
import os
import sys
from datetime import datetime
from pathlib import Path

input_dir = Path(sys.argv[1])
output_base = Path(sys.argv[2])
target = sys.argv[3]

run_dir = output_base / (datetime.now().strftime("%Y%m%d_%H%M%S") + f"_{os.getpid()}")
run_dir.mkdir(parents=True, exist_ok=False)


def read_record(handle):
    first = handle.readline()
    if first == "":
        return None
    return [first, handle.readline(), handle.readline(), handle.readline()]


def normalized_id(record):
    if record is None or not record[0]:
        return ""
    fields = record[0].strip().split()
    if not fields:
        return ""
    read_id = fields[0].lstrip("@")
    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]
    return read_id


def sequence_length(record):
    if record is None or not record[1]:
        return -1
    return len(record[1].rstrip("\r\n"))


def structure_issues(record, mate):
    issues = []
    if record is None:
        return [f"missing_{mate}"]
    if any(line == "" for line in record):
        issues.append(f"truncated_{mate}_record")
        return issues
    if not record[0].startswith("@"):
        issues.append(f"invalid_{mate}_header")
    if not record[2].startswith("+"):
        issues.append(f"invalid_{mate}_plus_line")
    sequence = record[1].rstrip("\r\n")
    quality = record[3].rstrip("\r\n")
    if len(sequence) != len(quality):
        issues.append(f"{mate}_sequence_quality_length_mismatch")
    return issues


def write_record(handle, record):
    if record is not None:
        for line in record:
            handle.write(line)


r1_files = sorted(input_dir.glob("*_R1_001.fastq.gz"))
if target != "all":
    r1_files = [
        path for path in r1_files
        if path.name.startswith(target + "_")
    ]

if not r1_files:
    raise SystemExit(f"No R1 FASTQ files matched target: {target}")

summary_path = run_dir / "pairing_summary.tsv"

with summary_path.open("w", newline="", encoding="utf-8") as summary_handle:
    summary_writer = csv.writer(summary_handle, delimiter="\t", lineterminator="\n")
    summary_writer.writerow([
        "sample",
        "positions_checked",
        "matching_ids",
        "id_mismatches",
        "mate_length_differences",
        "malformed_pairs",
        "missing_r1_records",
        "missing_r2_records",
        "status",
    ])

    for r1_path in r1_files:
        sample = r1_path.name.removesuffix("_L001_R1_001.fastq.gz")
        r2_path = input_dir / r1_path.name.replace(
            "_R1_001.fastq.gz", "_R2_001.fastq.gz"
        )

        if not r2_path.exists():
            summary_writer.writerow([
                sample, 0, 0, 0, 0, 0, 0, 1, "FAIL_MISSING_R2_FILE"
            ])
            print(f"{sample}: FAIL - R2 file is missing")
            continue

        details_path = run_dir / f"{sample}.flagged_pairs.tsv"
        flagged_r1_path = run_dir / f"{sample}.flagged_R1.fastq.gz"
        flagged_r2_path = run_dir / f"{sample}.flagged_R2.fastq.gz"

        positions = 0
        matching_ids = 0
        id_mismatches = 0
        length_differences = 0
        malformed_pairs = 0
        missing_r1 = 0
        missing_r2 = 0

        with (
            gzip.open(r1_path, "rt", encoding="utf-8", newline="") as r1_handle,
            gzip.open(r2_path, "rt", encoding="utf-8", newline="") as r2_handle,
            gzip.open(flagged_r1_path, "wt", encoding="utf-8", newline="") as flagged_r1,
            gzip.open(flagged_r2_path, "wt", encoding="utf-8", newline="") as flagged_r2,
            details_path.open("w", newline="", encoding="utf-8") as details_handle,
        ):
            details_writer = csv.writer(
                details_handle, delimiter="\t", lineterminator="\n"
            )
            details_writer.writerow([
                "pair_position",
                "r1_id",
                "r2_id",
                "r1_length",
                "r2_length",
                "issue",
            ])

            while True:
                r1_record = read_record(r1_handle)
                r2_record = read_record(r2_handle)

                if r1_record is None and r2_record is None:
                    break

                positions += 1
                issues = []

                if r1_record is None:
                    missing_r1 += 1
                    issues.append("missing_R1_record")
                if r2_record is None:
                    missing_r2 += 1
                    issues.append("missing_R2_record")

                r1_id = normalized_id(r1_record)
                r2_id = normalized_id(r2_record)
                r1_length = sequence_length(r1_record)
                r2_length = sequence_length(r2_record)

                if r1_record is not None and r2_record is not None:
                    if r1_id == r2_id:
                        matching_ids += 1
                    else:
                        id_mismatches += 1
                        issues.append("read_id_mismatch")

                    if r1_length != r2_length:
                        length_differences += 1
                        issues.append("mate_length_difference")

                format_issues = (
                    structure_issues(r1_record, "R1")
                    + structure_issues(r2_record, "R2")
                )
                if format_issues:
                    malformed_pairs += 1
                    issues.extend(format_issues)

                if issues:
                    details_writer.writerow([
                        positions,
                        r1_id,
                        r2_id,
                        r1_length,
                        r2_length,
                        ",".join(issues),
                    ])
                    write_record(flagged_r1, r1_record)
                    write_record(flagged_r2, r2_record)

        pairing_failures = id_mismatches + missing_r1 + missing_r2 + malformed_pairs
        if pairing_failures:
            status = "FAIL"
        elif length_differences:
            status = "PASS_WITH_LENGTH_DIFFERENCES"
        else:
            status = "PASS"

        summary_writer.writerow([
            sample,
            positions,
            matching_ids,
            id_mismatches,
            length_differences,
            malformed_pairs,
            missing_r1,
            missing_r2,
            status,
        ])

        print(
            f"{sample}: {status}; "
            f"positions={positions:,}; "
            f"ID mismatches={id_mismatches:,}; "
            f"length differences={length_differences:,}; "
            f"malformed={malformed_pairs:,}; "
            f"missing R1={missing_r1:,}; missing R2={missing_r2:,}"
        )

print(f"Results written to: {run_dir}")
print("Flagged FASTQ files contain records listed in each flagged_pairs.tsv file.")
PY
