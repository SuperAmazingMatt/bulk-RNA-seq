#!/bin/bash
module load fastqc
module load multiqc

input_dir="${BULK_PROJECT_DIR:-.}"/results/02_trimmomatic_results
output_dir="${BULK_PROJECT_DIR:-.}"/results/03_fastqc_trimmed_results
# Create output directory if it doesn't already exist (-p allows creation of nested directories)
mkdir -p "$output_dir"
for file in "$input_dir"/*_pe.fastq.gz
do
echo "$file"
fastqc "$file" --outdir="$output_dir"
done
multiqc "$output_dir" -o "$output_dir" -n multiqc_report.html --fullnames --force
