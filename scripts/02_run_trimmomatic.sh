#!/bin/bash
module load trimmomatic

input_dir="${BULK_PROJECT_DIR:-.}"/data
output_dir="${BULK_PROJECT_DIR:-.}"/results/02_trimmomatic_results
mkdir -p "$output_dir"
for file1 in ${input_dir}/*_R1_001.fastq.gz
do
file2="${file1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

if [[ ! -f "$file2" ]]; then
echo "Warning: Pair file not found for $file1" 
continue
fi
# Extract sample base name (e.g., groupA_rep1)
name=$(basename "$file1" "_L001_R1_001.fastq.gz")
# Define output file names
out1="${output_dir}/${name}_R1.trim_pe.fastq.gz"
out2="${output_dir}/${name}_R2.trim_pe.fastq.gz"
out1se="${output_dir}/${name}_R1.trim_se.fastq.gz"
out2se="${output_dir}/${name}_R2.trim_se.fastq.gz"
echo "Running Trimmomatic for $name..."
java -jar "$TRIMMOMATIC_PATH/trimmomatic-0.38.jar" PE -threads 4 "$file1" "$file2" "$out1" "$out1se" "$out2" "$out2se" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
done
