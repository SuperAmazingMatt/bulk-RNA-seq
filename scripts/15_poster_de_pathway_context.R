# Poster companion linking differential genes to enriched biological functions
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(GO.db)
library(tidyverse)

project_dir <- "."
output_dir <- file.path(project_dir, "results/15_poster_de_pathway_context_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

groupA_colour <- "#D95F59"
groupB_colour <- "#2F6B9A"
navy_colour <- "#0B2A4A"

poster_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 19, face = "bold", color = navy_colour),
    plot.subtitle = element_text(size = 12, color = "#52677B", margin = margin(b = 10)),
    plot.caption = element_text(size = 10.5, color = "#52677B", hjust = 0, margin = margin(t = 10)),
    axis.title = element_text(size = 13, face = "bold", color = navy_colour),
    axis.text = element_text(size = 11, color = navy_colour),
    strip.text = element_text(size = 10.8, face = "bold", color = navy_colour, lineheight = 0.95),
    strip.background = element_rect(fill = "#EEF3F6", color = NA),
    plot.margin = margin(14, 18, 14, 18)
  )

save_poster_plot <- function(plot_object, filename, width, height) {
  ggsave(
    file.path(output_dir, paste0(filename, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )

  grDevices::svg(
    file.path(output_dir, paste0(filename, ".svg")),
    width = width,
    height = height,
    bg = "white",
    onefile = TRUE
  )
  print(plot_object)
  grDevices::dev.off()
}

# Reconstruct the verified Stage 08 model and Stage 11 GO-BP membership
counts <- read.table(
  file.path(project_dir, "results/06_feature_counts_results/feature_counts.txt"),
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- counts[, c("Geneid", colnames(counts)[7:ncol(counts)])]
colnames(counts)[-1] <- c(
  "groupA_rep1",
  "groupA_rep2",
  "groupA_rep3",
  "groupB_rep1",
  "groupB_rep2",
  "groupB_rep3"
)
colnames(counts)[1] <- "GeneName"
counts$GeneName <- make.unique(as.character(counts$GeneName), sep = "_")
rownames(counts) <- counts$GeneName
counts$GeneName <- NULL

dge <- DGEList(counts = counts)
group <- factor(c("groupA", "groupA", "groupA", "groupB", "groupB", "groupB"))
dge$samples$group <- group
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~0 + group)
colnames(design) <- c("groupA", "groupB")
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
contrast_groupA_vs_groupB <- makeContrasts(groupA - groupB, levels = design)
de_test <- glmQLFTest(fit, contrast = contrast_groupA_vs_groupB)
de_results <- topTags(de_test, n = Inf)$table
de_results$ENSEMBL <- rownames(de_results)

all_genes_in_db <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
go_mapping <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = all_genes_in_db,
  keytype = "ENSEMBL",
  columns = c("GO", "ONTOLOGY")
)

bp_mapping <- go_mapping %>%
  dplyr::filter(ONTOLOGY == "BP" & !is.na(GO) & !is.na(ENSEMBL))

gene_sets_list <- split(bp_mapping$ENSEMBL, bp_mapping$GO)
indexed_sets <- ids2indices(gene_sets_list, identifiers = rownames(dge))
indexed_sets <- indexed_sets[sapply(indexed_sets, length) >= 5]

camera_results <- read.csv(
  file.path(
    project_dir,
    "results/11_camera_pathway_results/CAMERA_Competitive_Enrichment_groupA_vs_groupB.csv"
  ),
  row.names = 1,
  check.names = FALSE
)

camera_results$GO_ID <- rownames(camera_results)
camera_results$Term <- AnnotationDbi::Term(GO.db::GOTERM[camera_results$GO_ID])
camera_results$Term <- ifelse(is.na(camera_results$Term), camera_results$GO_ID, camera_results$Term)

# Supply study-appropriate terms locally; the original selected terms are withheld.
selection_file <- file.path(project_dir, "config/pathway_selection.R")
if (!file.exists(selection_file)) {
  stop("Copy config/pathway_selection.example.R to config/pathway_selection.R and supply your GO IDs.")
}
selection_env <- new.env(parent = baseenv())
sys.source(selection_file, envir = selection_env)
representative_go_ids <- selection_env$representative_go_ids
if (!is.character(representative_go_ids) || length(representative_go_ids) != 4L ||
    anyNA(representative_go_ids) || anyDuplicated(representative_go_ids) ||
    !all(grepl("^GO:[0-9]{7}$", representative_go_ids))) {
  stop("Provide four distinct GO IDs in the local pathway selection file.")
}
if (!all(representative_go_ids %in% camera_results$GO_ID)) {
  stop("Each selected GO ID must occur in your Stage 11 CAMERA results.")
}

selected_pathways <- camera_results %>%
  dplyr::filter(GO_ID %in% representative_go_ids) %>%
  mutate(selection_order = match(GO_ID, representative_go_ids)) %>%
  arrange(selection_order) %>%
  mutate(
    direction_label = ifelse(Direction == "Up", "Higher in groupA", "Higher in groupB"),
    pathway_label = paste0(
      stringr::str_wrap(Term, width = 39),
      "\nCAMERA FDR = ",
      format(FDR, digits = 2, scientific = TRUE)
    )
  )

selected_members <- purrr::map_dfr(seq_len(nrow(selected_pathways)), function(pathway_row) {
  go_id <- selected_pathways$GO_ID[pathway_row]
  member_ids <- rownames(dge)[indexed_sets[[go_id]]]
  expected_direction <- selected_pathways$Direction[pathway_row]

  de_results %>%
    dplyr::filter(
      ENSEMBL %in% member_ids,
      abs(logFC) > 2,
      PValue < 0.05,
      if (expected_direction == "Up") logFC > 0 else logFC < 0
    ) %>%
    arrange(PValue) %>%
    slice_head(n = 3) %>%
    mutate(
      GO_ID = go_id,
      pathway_label = selected_pathways$pathway_label[pathway_row],
      direction_label = selected_pathways$direction_label[pathway_row]
    )
})

gene_symbols <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = unique(selected_members$ENSEMBL),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

selected_members <- selected_members %>%
  mutate(
    gene_symbol = unname(gene_symbols[ENSEMBL]),
    gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "", ENSEMBL, gene_symbol),
    gene_axis = paste0(gene_symbol, "|", GO_ID),
    direction_label = factor(
      direction_label,
      levels = c("Higher in groupB", "Higher in groupA")
    )
  )

gene_axis_levels <- selected_members %>%
  arrange(pathway_label, logFC) %>%
  pull(gene_axis) %>%
  unique()
selected_members$gene_axis <- factor(selected_members$gene_axis, levels = gene_axis_levels)

pathway_context_plot <- ggplot(
  selected_members,
  aes(x = logFC, y = gene_axis, color = direction_label)
) +
  geom_vline(xintercept = 0, color = "#7A8A99", linewidth = 0.7) +
  geom_segment(aes(x = 0, xend = logFC, yend = gene_axis), linewidth = 1.1, alpha = 0.45) +
  geom_point(size = 4.8) +
  geom_text(
    aes(label = sprintf("%.1f", logFC)),
    hjust = ifelse(selected_members$logFC > 0, -0.55, 1.55),
    color = navy_colour,
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~pathway_label, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c(
    "Higher in groupB" = groupB_colour,
    "Higher in groupA" = groupA_colour
  )) +
  scale_y_discrete(labels = function(x) sub("\\|.*$", "", x)) +
  scale_x_continuous(expand = expansion(mult = c(0.20, 0.20))) +
  labs(
    title = "Differential genes within selected gene sets",
    subtitle = "Up to three concordant members from four representative FDR-significant CAMERA terms",
    x = "Group B higher  ←  log2 fold change (groupA - groupB)  →  Group A higher",
    y = "Selected member genes",
    color = "Expression direction",
    caption = stringr::str_wrap(
      paste0(
        "Four terms were supplied in the local pathway selection file. ",
        "Genes pass the teaching rule |log2FC| > 2 and raw P < 0.05 and are selected by smallest raw P within each set. ",
        "They are examples, not pathway drivers; related GO sets can overlap."
      ),
      width = 145
    )
  ) +
  poster_theme +
  theme(
    panel.spacing = unit(1.1, "lines"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom"
  )

save_poster_plot(
  pathway_context_plot,
  "DE_gene_to_CAMERA_pathway_context",
  width = 13,
  height = 9
)

cat("Poster DE/pathway companion written to:", output_dir, "\n")
print(
  selected_members %>%
    dplyr::select(GO_ID, pathway_label, gene_symbol, logFC, PValue, direction_label)
)
print(sort(list.files(output_dir)))
