#Load Pre-installed Libraries
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tidyverse)
#Fix function conflicts explicitly (Tidyverse vs AnnotationDbi)
select <- AnnotationDbi::select
filter <- dplyr::filter
# Data Loading
counts <- read.table("results/06_feature_counts_results/feature_counts.txt", header=TRUE, sep="\t", comment.char="#",
                     check.names=FALSE, stringsAsFactors = FALSE)
counts <- counts[,c("Geneid",colnames(counts)[7:ncol(counts)])]
colnames(counts)[-1] <- c("groupA_rep1", "groupA_rep2",
                          "groupA_rep3", "groupB_rep1",
                          "groupB_rep2", "groupB_rep3")
colnames(counts)[1] <- "GeneName"
counts$GeneName <- make.unique(as.character(counts$GeneName), sep = "_")
rownames(counts) <- counts$GeneName
counts$GeneName <- NULL
# DGEList Setup & Normalisation
dge <- DGEList(counts = counts)
group <- factor(c("groupA", "groupA", "groupA", "groupB", "groupB", "groupB"))
dge$samples$group <- group
cpm_values <- cpm(dge)
keep <- rowSums(cpm_values > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
#Differential Expression Analysis
design <- model.matrix(~0 + group)
colnames(design) <- c("groupA", "groupB")
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
contrast_groupA_vs_groupB <- makeContrasts(groupA - groupB, levels = design)
qlf_groupA_vs_groupB <- glmQLFTest(fit, contrast = contrast_groupA_vs_groupB)
results_groupA_vs_groupB <- topTags(qlf_groupA_vs_groupB, n = Inf)$table
#Isolate Significant DE Genes
results_groupA_vs_groupB$Significant <- "No"
results_groupA_vs_groupB$Significant[results_groupA_vs_groupB$PValue < 0.05 &
                                       results_groupA_vs_groupB$logFC > 2] <- "Up"
results_groupA_vs_groupB$Significant[results_groupA_vs_groupB$PValue < 0.05 &
                                       results_groupA_vs_groupB$logFC < -2] <- "Down"
sig_groupA_vs_groupB <- results_groupA_vs_groupB[results_groupA_vs_groupB$PValue < 0.05 &
                                                     abs(results_groupA_vs_groupB$logFC) > 2, ]

#Annotation Lookup (Convert Ensembl IDs to Entrez IDs)
test_annotation <- select(org.Mm.eg.db, keys = rownames(sig_groupA_vs_groupB), keytype =
                            "ENSEMBL", columns = c("ENSEMBL", "ENTREZID", "GENENAME"))
entrez_genes <- test_annotation$ENTREZID
entrez_genes <- entrez_genes[!is.na(entrez_genes)]
entrez_genes <- as.character(entrez_genes)
de_genes <- list(`groupA_vs_groupB` = entrez_genes)
#KEGG Pathway Enrichment
#Run KEGG enrichment mapping
kegg_results <- kegga(de_genes, species = "Mm")
#Adjust for multiple testing using Benjamini-Hochberg false discovery rate
kegg_results$FDR <- p.adjust(kegg_results$P.groupA_vs_groupB, method = "BH")
kegg_results <- kegg_results %>% arrange(FDR)
#Filter for statistically significant biological pathways
sig_kegg <- kegg_results %>% filter(FDR < 0.05)
#Save KEGG Pathway Results
write.csv(sig_kegg, "results/10_kegg_pathway_results/KEGG_Pathways_groupA_vs_groupB.csv", row.names = TRUE)

