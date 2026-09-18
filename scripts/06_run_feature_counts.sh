#!/bin/bash
module load subread
module load multiqc

reference_gtf="${REFERENCE_GTF:-./reference/annotation.gtf}"

bam_dir="${BULK_PROJECT_DIR:-.}"/results/05_samtools_results
bam_ending=".sorted.bam"
output_dir="${BULK_PROJECT_DIR:-.}"/results/06_feature_counts_results
output_file=feature_counts.txt
mkdir -p $output_dir
featureCounts -T 4 -p --countReadPairs -t exon -g gene_id -a $reference_gtf -s 0 -o $output_dir/$output_file $bam_dir/*$bam_ending

# display the summary report 
multiqc "$output_dir" -o "$output_dir" -n feature_counts_summary.html --force 
