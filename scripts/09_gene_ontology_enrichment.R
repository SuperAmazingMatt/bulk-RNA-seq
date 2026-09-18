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



#links the Ensembl IDs from the differential expression analysis to standard Entrez identifiers, which can then be used for downstream pathway and gene-set enrichment analysis. Ensembl IDs → Entrez IDs → pathway/gene-set annotation → enrichment analysis
test_annotation



entrez_genes <- test_annotation$ENTREZID
entrez_genes <- entrez_genes[!is.na(entrez_genes)]
entrez_genes <- as.character(entrez_genes)
de_genes <- list(`groupA_vs_groupB` = entrez_genes)
de_genes



#Gene Ontology (GO) Enrichment via goana
go_results <- goana(de_genes, species = "Mm")
go_results


go_results$FDR <- p.adjust(go_results$P.groupA_vs_groupB, method = "BH")
go_results


go_results <- go_results %>% arrange(FDR)
go_results



bp <- go_results %>% filter(Ont == "BP")
mf <- go_results %>% filter(Ont == "MF")
cc <- go_results %>% filter(Ont == "CC")
write.csv(bp %>% filter(FDR < 0.05), "results/09_gene_ontology_results/Biological_Process_groupA_vs_groupB.csv",
          row.names = TRUE)
write.csv(mf %>% filter(FDR < 0.05), "results/09_gene_ontology_results/Molecular_Function_groupA_vs_groupB.csv",
          row.names = TRUE)
write.csv(cc %>% filter(FDR < 0.05), "results/09_gene_ontology_results/Cellular_Component_groupA_vs_groupB.csv",
          row.names = TRUE)

