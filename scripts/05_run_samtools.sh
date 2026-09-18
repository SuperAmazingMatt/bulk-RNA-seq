#!/bin/bash
module load samtools

input_dir="${BULK_PROJECT_DIR:-.}"/results/04_hisat2_results
output_dir="${BULK_PROJECT_DIR:-.}"/results/05_samtools_results
mkdir -p $output_dir
for file in $input_dir/*.sam
do
# Get the basename of the file minus the ".sam" ending
name=$(basename $file .hisat2.sam)
echo $name
# Convert SAM to BAM format
samtools view -@ 4 -b $file > $output_dir/$name.unsorted.bam
# Sort BAM files
samtools sort -@ 4 $output_dir/$name.unsorted.bam > $output_dir/$name.sorted.bam
# Index BAM files
samtools index -@ 4 $output_dir/$name.sorted.bam
# Generate summary stats
samtools flagstat -@ 4 $output_dir/$name.sorted.bam > $output_dir/$name.flagstat
# Delete the unsorted BAM
rm $output_dir/$name.unsorted.bam
done
