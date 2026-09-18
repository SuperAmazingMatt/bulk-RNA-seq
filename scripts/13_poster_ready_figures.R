# Publication adaptation of the preserved local Stage 13 source; see scripts/README.md.
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(GO.db)
library(tidyverse)
library(ggrepel)

project_dir <- "."
output_dir <- file.path(project_dir, "results/13_poster_figures_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

groupA_colour <- "#D95F59"
groupB_colour <- "#2F6B9A"
teal_colour <- "#008A8F"
navy_colour <- "#0B2A4A"
grey_colour <- "#C7CED6"

poster_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 19, face = "bold", color = navy_colour),
    plot.subtitle = element_text(size = 12, color = "#52677B", margin = margin(b = 10)),
    axis.title = element_text(size = 13, face = "bold", color = navy_colour),
    axis.text = element_text(size = 11, color = navy_colour),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(14, 18, 14, 18)
  )

save_poster_plot <- function(plot_object, filename, width, height) {
  png_file <- file.path(output_dir, paste0(filename, ".png"))
  svg_file <- file.path(output_dir, paste0(filename, ".svg"))

  ggsave(
    png_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )

  grDevices::svg(
    svg_file,
    width = width,
    height = height,
    bg = "white",
    onefile = TRUE
  )
  print(plot_object)
  grDevices::dev.off()
}

# Reconstruct the verified count import and filtering used in Stages 07 and 08
counts <- read.delim(
  file.path(project_dir, "results/06_feature_counts_results/feature_counts.txt"),
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE
)

counts <- counts[, c("Geneid", colnames(counts)[7:ncol(counts)])]
sample_ids <- c(
  "groupA_rep1",
  "groupA_rep2",
  "groupA_rep3",
  "groupB_rep1",
  "groupB_rep2",
  "groupB_rep3"
)
colnames(counts) <- c("genes", sample_ids)

dge <- DGEList(counts = counts[, -1], genes = counts[, 1])
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, ]

# PCA: preserve the teaching-template calculation before TMM
logCPM <- cpm(dge, log = TRUE, prior.count = 1)
pca <- prcomp(t(logCPM))
pca_variance <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

pca_df <- data.frame(
  sample = c("Group A 1", "Group A 2", "Group A 3", "Group B 1", "Group B 2", "Group B 3"),
  condition = factor(
    c(rep("Group A", 3), rep("Group B", 3)),
    levels = c("Group A", "Group B")
  ),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_hline(yintercept = 0, color = "#E2E7EC", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "#E2E7EC", linewidth = 0.5) +
  geom_point(size = 4.5, stroke = 1) +
  geom_text_repel(
    aes(label = sample),
    size = 4,
    color = navy_colour,
    box.padding = 0.5,
    point.padding = 0.35,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Group A" = groupA_colour,
    "Group B" = groupB_colour
  )) +
  scale_shape_manual(values = c(
    "Group A" = 16,
    "Group B" = 17
  )) +
  labs(
    title = "PCA of sample relationships",
    subtitle = paste0(
      format(sum(keep), big.mark = ","),
      " genes; library-size log2 CPM calculated before TMM"
    ),
    x = sprintf("PC1 (%.2f%%)", pca_variance[1]),
    y = sprintf("PC2 (%.2f%%)", pca_variance[2]),
    color = "Condition",
    shape = "Condition"
  ) +
  poster_theme

save_poster_plot(
  pca_plot,
  "PCA_sample_structure_groupA_vs_groupB",
  width = 10,
  height = 7
)

# Differential expression: preserve the Stage 08 model and teaching thresholds
x <- calcNormFactors(dge)
experiment_groups <- factor(
  sub("_.*", "", colnames(dge$counts)),
  levels = c("groupA", "groupB")
)
x$samples$group <- experiment_groups

design <- model.matrix(~0 + experiment_groups)
colnames(design) <- sub("experiment_groups", "", colnames(design))
rownames(design) <- rownames(x$samples)

contrasts <- makeContrasts(
  groupA_vs_groupB = groupA - groupB,
  levels = design
)

x <- estimateDisp(x, design)
fit <- glmQLFit(x, design)
res <- glmQLFTest(fit, contrast = contrasts[, "groupA_vs_groupB"])
result_table <- topTags(res, n = Inf)$table

gene_ids <- if ("genes" %in% colnames(result_table)) {
  as.character(result_table$genes)
} else {
  rownames(result_table)
}
gene_ids <- sub("\\..*$", "", gene_ids)
gene_symbols <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = unique(gene_ids),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

volcano_df <- result_table %>%
  mutate(
    gene_id = gene_ids,
    gene_symbol = unname(gene_symbols[gene_ids]),
    gene_label = ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol),
    neg_log10_p = -log10(pmax(PValue, .Machine$double.xmin)),
    expression = case_when(
      logFC > 2 & PValue < 0.05 ~ "Higher in groupA",
      logFC < -2 & PValue < 0.05 ~ "Higher in groupB",
      TRUE ~ "Not selected"
    )
  )

volcano_df$expression <- factor(
  volcano_df$expression,
  levels = c("Not selected", "Higher in groupB", "Higher in groupA")
)

label_genes <- bind_rows(
  volcano_df %>%
    filter(expression == "Higher in groupA") %>%
    arrange(PValue) %>%
    slice_head(n = 6),
  volcano_df %>%
    filter(expression == "Higher in groupB") %>%
    arrange(PValue) %>%
    slice_head(n = 6)
) %>%
  distinct(gene_label, .keep_all = TRUE)

volcano_counts <- table(volcano_df$expression)
volcano_x_limit <- ceiling(max(abs(volcano_df$logFC), na.rm = TRUE))

volcano_plot <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_p, color = expression)) +
  geom_point(alpha = 0.58, size = 1.15) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "#52677B") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#52677B") +
  geom_text_repel(
    data = label_genes,
    aes(label = gene_label),
    color = navy_colour,
    size = 3.5,
    box.padding = 0.45,
    point.padding = 0.25,
    segment.color = "#8A98A6",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Not selected" = grey_colour,
    "Higher in groupB" = groupB_colour,
    "Higher in groupA" = groupA_colour
  )) +
  coord_cartesian(xlim = c(-volcano_x_limit, volcano_x_limit)) +
  labs(
    title = "Differential expression: groupA versus groupB",
    subtitle = paste0(
      format(unname(volcano_counts["Higher in groupA"]), big.mark = ","),
      " higher in groupA; ",
      format(unname(volcano_counts["Higher in groupB"]), big.mark = ","),
      " higher in groupB; threshold: |log2FC| > 2 and raw P < 0.05"
    ),
    x = "log2 fold change (groupA - groupB)",
    y = "-log10(raw P-value)",
    color = "Expression"
  ) +
  poster_theme

save_poster_plot(
  volcano_plot,
  "Volcano_differential_expression_groupA_vs_groupB",
  width = 10,
  height = 7
)

# KEGG: top pathways from the existing combined-direction Stage 10 result
kegg_results <- read.csv(
  file.path(project_dir, "results/10_kegg_pathway_results/KEGG_Pathways_groupA_vs_groupB.csv"),
  row.names = 1,
  check.names = FALSE
)

kegg_plot_df <- kegg_results %>%
  mutate(
    KEGG_ID = rownames(kegg_results),
    neg_log10_fdr = -log10(pmax(FDR, .Machine$double.xmin)),
    pathway_label = stringr::str_wrap(Pathway, width = 38),
    de_genes = groupA_vs_groupB
  ) %>%
  arrange(FDR) %>%
  slice_head(n = 10)

kegg_plot_df$pathway_label <- factor(
  kegg_plot_df$pathway_label,
  levels = rev(kegg_plot_df$pathway_label)
)

kegg_plot <- ggplot(kegg_plot_df, aes(x = neg_log10_fdr, y = pathway_label)) +
  geom_segment(
    aes(x = 0, xend = neg_log10_fdr, yend = pathway_label),
    color = "#DDE5EB",
    linewidth = 0.8
  ) +
  geom_point(aes(size = de_genes), color = teal_colour, alpha = 0.9) +
  scale_size_continuous(range = c(4, 10)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.08))) +
  labs(
    title = "KEGG pathway over-representation",
    subtitle = "Top 10 pathways by FDR; up- and downregulated DE genes were pooled",
    x = "-log10(FDR)",
    y = NULL,
    size = "DE genes in pathway\n(combined directions)"
  ) +
  poster_theme +
  theme(axis.text.y = element_text(size = 11, lineheight = 0.95))

save_poster_plot(
  kegg_plot,
  "KEGG_pathway_enrichment_groupA_vs_groupB",
  width = 11,
  height = 8
)

# CAMERA: show the strongest terms in both expression directions
camera_results <- read.csv(
  file.path(project_dir, "results/11_camera_pathway_results/CAMERA_Competitive_Enrichment_groupA_vs_groupB.csv"),
  row.names = 1,
  check.names = FALSE
)

camera_results$GO_ID <- rownames(camera_results)
camera_results$Term <- AnnotationDbi::Term(GO.db::GOTERM[camera_results$GO_ID])
camera_results$Term <- ifelse(
  is.na(camera_results$Term),
  camera_results$GO_ID,
  camera_results$Term
)

camera_plot_df <- bind_rows(
  camera_results %>%
    filter(Direction == "Up") %>%
    arrange(FDR) %>%
    slice_head(n = 6),
  camera_results %>%
    filter(Direction == "Down") %>%
    arrange(FDR) %>%
    slice_head(n = 6)
) %>%
  mutate(
    direction_label = ifelse(Direction == "Up", "Higher in groupA", "Higher in groupB"),
    signed_log10_fdr = ifelse(Direction == "Up", -log10(FDR), log10(FDR)),
    term_label = stringr::str_wrap(Term, width = 43)
  ) %>%
  arrange(signed_log10_fdr)

camera_plot_df$term_label <- factor(
  camera_plot_df$term_label,
  levels = camera_plot_df$term_label
)

camera_plot <- ggplot(camera_plot_df, aes(x = signed_log10_fdr, y = term_label)) +
  geom_vline(xintercept = 0, color = "#52677B", linewidth = 0.7) +
  geom_segment(
    aes(x = 0, xend = signed_log10_fdr, yend = term_label, color = direction_label),
    linewidth = 0.9,
    alpha = 0.5
  ) +
  geom_point(aes(size = NGenes, color = direction_label), alpha = 0.95) +
  scale_color_manual(values = c(
    "Higher in groupA" = groupA_colour,
    "Higher in groupB" = groupB_colour
  )) +
  scale_size_continuous(range = c(4, 10)) +
  scale_x_continuous(
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  labs(
    title = "CAMERA competitive enrichment: GO Biological Process",
    subtitle = "Top 6 terms within each direction by FDR",
    x = "Group B higher  <-  -log10(FDR)  ->  Group A higher",
    y = NULL,
    color = "Expression direction",
    size = "Genes tested in set"
  ) +
  poster_theme +
  theme(axis.text.y = element_text(size = 10.5, lineheight = 0.95))

save_poster_plot(
  camera_plot,
  "CAMERA_GO_BP_enrichment_groupA_vs_groupB",
  width = 12,
  height = 8.5
)

cat("Poster figures written to:", output_dir, "\n")
print(sort(list.files(output_dir)))
