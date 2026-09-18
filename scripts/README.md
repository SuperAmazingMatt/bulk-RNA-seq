# Bulk RNA-seq scripts

These 18 scripts demonstrate a complete bulk RNA-seq workflow: stages 01–15, two supplementary paired-read checks, and a separately labelled alignment revision. They are sanitized adaptations of the preserved working scripts. Personal paths, account details, study labels and embedded findings have been removed; input data and generated results are not included.

## Workflow

| Stage | Script | Main tools and purpose |
|---|---|---|
| 01 | [Raw-read QC](01_run_fastqc_raw_data.sh) | FastQC and MultiQC: inspect read and library quality |
| 02 | [Trimming](02_run_trimmomatic.sh) | Trimmomatic: paired-end quality trimming |
| 03 | [Post-trimming QC](03_run_fastqc_trimmed_data.sh) | FastQC and MultiQC: review the trimmed paired reads |
| 04 | [Original alignment script](04_run_HISAT2.sh) | HISAT2: align paired reads; historical naming limitation described below |
| 05 | [BAM processing](05_run_samtools.sh) | SAMtools: convert, sort, index and summarize alignments |
| 06 | [Gene counting](06_run_feature_counts.sh) | featureCounts and MultiQC: count paired fragments against gene annotations |
| 07 | [Normalization and PCA](07_normalization_pca.Rmd) | edgeR and R: filter low-expression genes, inspect PCA and calculate TMM factors |
| 08 | [Differential expression](08_differential_expression.Rmd) | edgeR: fit a quasi-likelihood model, test a contrast and draw a volcano plot |
| 09 | [GO analysis](09_gene_ontology_enrichment.R) | edgeR, AnnotationDbi and org.Mm.eg.db: map selected genes and test GO categories |
| 10 | [KEGG analysis](10_kegg_pathway_enrichment.R) | edgeR and annotation tools: pathway over-representation |
| 11 | [Competitive gene-set testing](11_camera_pathway_enrichment.R) | limma CAMERA: compare gene sets within the expression model |
| 12 | [Pathway visualization](12_visualise_enriched_pathways.R) | ggplot2, GO.db and tidyverse: display the strongest CAMERA terms |
| 13 | [Result figures](13_poster_ready_figures.R) | ggplot2 and ggrepel: assemble PCA, volcano and enrichment displays |
| 14 | [QC figures](14_poster_qc_problem_figures.R) | edgeR and tidyverse: compare complexity, summarize supplied QC evidence and visualize assignment |
| 15 | [Gene-to-pathway context](15_poster_de_pathway_context.R) | edgeR, limma and annotation tools: show selected member genes within supplied gene sets |

Supplementary scripts:

- [Bash/Python pairing check](check_fastq_pairing.sh): validate FASTQ structure, matching read IDs and mate lengths. It defaults to all supplied libraries and writes reports and flagged reads under `results/pairing_check/`.
- [R pairing and length analysis](check_and_plot_pair_lengths.R): scan raw or trimmed FASTQ in chunks, distinguish integrity failures from length differences, and summarize paired lengths. Failure-read export is optional.
- [Reviewed alignment revision](reviewed_not_run/04_run_HISAT2_v2.sh): correct output-name parsing, check mates, capture logs and require a new output directory. This revision was reviewed but was not the version that produced the original results.

## Inputs and local setup

Run from the repository root in an environment with the required tools. The shell scripts expect an environment-module system; module names and the Trimmomatic JAR location may need adjustment. R scripts use edgeR, limma, AnnotationDbi, org.Mm.eg.db, GO.db, tidyverse, ggrepel and ggfortify as indicated in each file. The mouse annotation and reference must describe the same assembly.

- Place your own paired FASTQ inputs in `data/`, using names such as `groupA_rep1_L001_R1_001.fastq.gz` and the matching R2 name. The six example sample slots are `groupA_rep1` through `groupA_rep3`, then `groupB_rep1` through `groupB_rep3`.
- Provide a HISAT2 index through `HISAT2_INDEX` (default `./reference/hisat2_index/genome`) and a GTF through `REFERENCE_GTF` (default `./reference/annotation.gtf`). `BULK_PROJECT_DIR` overrides the numbered shell scripts' default project root `.`; R paths remain relative to the current repository root. The supplementary pairing shell script derives the project root from its own `scripts/` location.
- Check the count-column order against your own sample metadata before analysis. Stages 07–15 assign the six example sample names by column position. Change the design and labels deliberately for a different dataset. Positive log fold change means higher in group A for the group A minus group B contrast.
- Inputs and outputs follow the numbered `results/` directories named in each script. Create each R stage's output directory before running it. Stages 13–15 create their own output directories.
- For Stage 14, copy [plot_settings.example.R](../config/plot_settings.example.R) to `config/plot_settings.R` and complete it using your own QC evidence. The example contains missing values, not substitute findings.
- For Stage 15, copy [pathway_selection.example.R](../config/pathway_selection.example.R) to `config/pathway_selection.R` and supply four GO IDs from your own Stage 11 results. The original study's selected terms are withheld.

Both completed configuration files, `data/`, `reference/` and `results/` are ignored by Git. Review outputs before sharing them because executing the scripts on your own data will generate labels, measurements and reports.

The R Markdown files have analysis evaluation disabled. Inspect the chunks in RStudio with the working directory set to the repository root; enable or execute them only when you intend to run the analysis. When rendering an enabled Rmd, explicitly use the repository root as `knit_root_dir`. These publication copies have not been executed on data.

## Provenance and interpretation

The original Stage 04 script removes a suffix that omits `.gz`, so its output names require correction before downstream reuse. The historical output files were renamed afterward. Prefer the separately labelled reviewed revision for a new run; set `OUTPUT_DIR` to a fresh Stage 04 directory and point Stage 05 at that directory. The revision must not be represented as the historically executed version.

Stage 13 is adapted from a preserved local source with a different fingerprint from the recorded executed server version. Stages 13–15 are display-oriented scripts that reconstruct some calculations, so they are not guaranteed to reproduce earlier results numerically. In particular, Stage 15 recalculates library sizes after filtering, whereas Stage 08 retains its original subsetting behavior.

PCA is calculated before TMM normalization in the preserved workflow. Differential-gene selection uses strict `abs(logFC) > 2` and unadjusted `PValue < 0.05`; this is not an FDR-controlled gene list. The original GO/KEGG code does not specify an explicit tested-gene universe and retains one-to-many identifier mappings. Those limitations remain visible rather than being silently changed in a publication copy.

The repository contains pipeline code and the separately redacted poster portfolio. It excludes raw reads, matrices, measured result tables, unredacted plots, editable posters, private notebooks, document-build helpers, superseded alternatives and original teaching-template files.
