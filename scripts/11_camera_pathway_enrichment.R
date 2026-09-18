#Load Pre-installed Libraries
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tidyverse)
#Fix function conflicts explicitly (Tidyverse vs AnnotationDbi)
#Different R packages can contain functions with the same name. The :: operator tells R
#exactly which package's function to use, preventing conflicts and making the code easier to
#understand and reproduce. Here it will use the function select from AnnotationDbi etc.
select <- AnnotationDbi::select
filter <- dplyr::filter
#Data Loading and Structuring
counts <- read.table("results/06_feature_counts_results/feature_counts.txt", header=TRUE,
                     sep="\t", comment.char="#", check.names=FALSE, stringsAsFactors = FALSE)
counts <- counts[,c("Geneid",colnames(counts)[7:ncol(counts)])]
colnames(counts)[-1] <- c("groupA_rep1", "groupA_rep2",
                          "groupA_rep3", "groupB_rep1",
                          "groupB_rep2", "groupB_rep3")
counts


#rename the first column of the counts data frame to "GeneName".
colnames(counts)[1] <- "GeneName"
#adds _1 to duplicated names
counts$GeneName <- make.unique(as.character(counts$GeneName), sep = "_")
#Use gene names as row names
rownames(counts) <- counts$GeneName
rownames(counts)



#Delete the GeneName column from the counts data frame. Gene names have already been
#moved from the GeneName column into the row names.
counts$GeneName <- NULL
#Brief Note: Rename the gene column → make gene names unique → move gene names to
#row names → remove the now-redundant gene column.
#DGEList Setup & Normalisation
#Create an edgeR object containing the RNA-seq count data and store it as dge.
dge <- DGEList(counts = counts)
dge


#Creates a grouping variable that tells edgeR which experimental condition each sample belongs to
group <- factor(c("groupA", "groupA", "groupA", "groupB", "groupB", "groupB"))
group



dge$samples$group <- group
#Calculates the Counts Per Million (CPM) for every gene in every sample. Raw RNA-seq counts are affected by library size.
cpm_values <- cpm(dge)
cpm_values



keep <- rowSums(cpm_values > 1) >= 2
keep


dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
dge




#Differential Expression Analysis
design <- model.matrix(~0 + group)
colnames(design) <- c("groupA", "groupB")
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
fit



contrast_groupA_vs_groupB <- makeContrasts(groupA - groupB, levels = design)


#Extract gene mapping from database
all_genes_in_db <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
go_mapping <- select(org.Mm.eg.db, keys = all_genes_in_db, keytype = "ENSEMBL",
                     columns = c("GO", "ONTOLOGY"))
go_mapping




#Filter strictly down to Biological Processes (BP)
bp_mapping <- go_mapping %>% filter(ONTOLOGY == "BP" & !is.na(GO) &
                                      !is.na(ENSEMBL))
bp_mapping




#Create our structured gene list mapped by Ensembl IDs to match DGEList matrix directly
gene_sets_list <- split(bp_mapping$ENSEMBL, bp_mapping$GO)
#Index matching positions relative to our active matrix rows
indexed_sets <- ids2indices(gene_sets_list, identifiers = rownames(dge))
#Explanation for the above step : This step specifies which rows of the expression/count matrix belong to each gene set. Essentially, it maps the genes in each biological gene set to their corresponding row positions in the RNA-seq expression matrix. camera() needs to know which rows of your expression/count matrix belong to each gene set.
# ids2indices() converts the gene names in each gene set into the row numbers corresponding to those genes in dge.
#ids2indices() maps the genes in each biological gene set to their corresponding row positions in the RNA-seq expression matrix. This ensures that camera() tests the correct genes against the differential expression model.
#Gene sets → ids2indices() → camera() → enriched pathways → FDR → biological interpretation


#Drop tiny categories to keep calculations statistically valid. Filter Out Very Small Gene Sets
indexed_sets <- indexed_sets[sapply(indexed_sets, length) >= 5]
# Explanation of the above step:
# • indexed_sets contains the gene sets after matching their genes to those present in the RNA-seq dataset.
# • sapply(indexed_sets, length) calculates the number of genes from the RNA-seq dataset present in each gene set.
# • >= 5 keeps only gene sets containing at least five genes from the dataset.
# • Very small gene sets can produce unstable or unreliable enrichment statistics, so filtering them out improves the robustness and interpretability of the analysis.

#Run CAMERA across the target contrast design
camera_results <- camera(dge, index = indexed_sets, design = design, contrast =
                           contrast_groupA_vs_groupB)
camera_results


sig_camera <- camera_results %>% filter(FDR < 0.05)
sig_camera


#Save result data matrix
write.csv(sig_camera, "results/11_camera_pathway_results/CAMERA_Competitive_Enrichment_groupA_vs_groupB.csv", row.names = TRUE)

