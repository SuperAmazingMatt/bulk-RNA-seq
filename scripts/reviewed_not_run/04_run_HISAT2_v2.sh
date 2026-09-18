#!/usr/bin/env bash
set -euo pipefail

module load hisat2

project_dir="${BULK_PROJECT_DIR:-.}"
reference_genome="${HISAT2_INDEX:-./reference/hisat2_index/genome}"
input_dir="${project_dir}/results/02_trimmomatic_results"

# Require an explicit new directory so this reviewed script cannot overwrite the
# historical Stage 04 results accidentally.
if [[ -z "${OUTPUT_DIR:-}" ]]; then
    echo "ERROR: Set OUTPUT_DIR to a new, versioned Stage 04 output directory." >&2
    exit 1
fi

output_dir="${OUTPUT_DIR}"
if [[ -e "$output_dir" ]]; then
    echo "ERROR: OUTPUT_DIR already exists; choose a new directory: $output_dir" >&2
    exit 1
fi

log_dir="${output_dir}/logs"
mkdir -p "$log_dir"

shopt -s nullglob
r1_files=("$input_dir"/*_R1.trim_pe.fastq.gz)
if (( ${#r1_files[@]} != 6 )); then
    echo "ERROR: Expected 6 trimmed R1 files, found ${#r1_files[@]} in $input_dir" >&2
    exit 1
fi

{
    printf 'run_started=%s\n' "$(date -Is)"
    printf 'hostname=%s\n' "$(hostname)"
    printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-not_set}"
    printf 'reference_index=%s\n' "$reference_genome"
    printf 'input_dir=%s\n' "$input_dir"
    printf 'output_dir=%s\n' "$output_dir"
    module list 2>&1
    hisat2 --version
} > "${log_dir}/run_manifest.txt" 2>&1

printf 'sample\tR1\tR2\tSAM\tsummary\tstderr_log\n' > "${log_dir}/sample_manifest.tsv"

for file1 in "${r1_files[@]}"; do
    file2="${file1%_R1.trim_pe.fastq.gz}_R2.trim_pe.fastq.gz"
    if [[ ! -f "$file2" ]]; then
        echo "ERROR: Missing R2 mate for $file1" >&2
        exit 1
    fi

    sample_name="${file1##*/}"
    sample_name="${sample_name%_R1.trim_pe.fastq.gz}"
    output_sam="${output_dir}/${sample_name}.hisat2.sam"
    summary_file="${log_dir}/${sample_name}.hisat2_summary.txt"
    stderr_log="${log_dir}/${sample_name}.hisat2.stderr.log"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample_name" "$file1" "$file2" "$output_sam" "$summary_file" "$stderr_log" \
        >> "${log_dir}/sample_manifest.tsv"

    echo "Aligning sample: $sample_name"
    hisat2 -p 4 --new-summary \
        --summary-file "$summary_file" \
        -x "$reference_genome" \
        -1 "$file1" \
        -2 "$file2" \
        -S "$output_sam" \
        2> "$stderr_log"

    if [[ ! -s "$output_sam" ]]; then
        echo "ERROR: HISAT2 produced an empty SAM for $sample_name" >&2
        exit 1
    fi
done

printf 'run_finished=%s\n' "$(date -Is)" >> "${log_dir}/run_manifest.txt"
echo "Stage 04 completed: $output_dir"
