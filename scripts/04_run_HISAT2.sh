#!/bin/bash
module load hisat2

reference_genome="${HISAT2_INDEX:-./reference/hisat2_index/genome}"

input_dir="${BULK_PROJECT_DIR:-.}"/results/02_trimmomatic_results
output_dir="${BULK_PROJECT_DIR:-.}"/results/04_hisat2_results
mkdir -p $output_dir
for file1 in "$input_dir"/*_R1.trim_pe.fastq.gz
do
file2="${file1/_R1.trim_pe.fastq.gz/_R2.trim_pe.fastq.gz}"
if [[ ! -f "$file2" ]]; then
echo "Paired file not found for: $file1"
continue
fi
sample_name=$(basename "$file1" _R1.trim_pe.fastq)
output_sam="${output_dir}/${sample_name}.hisat2.sam"
echo "Aligning sample: $sample_name"

hisat2 -x "$reference_genome" -1 "$file1" -2 "$file2" -S "$output_sam"
done
