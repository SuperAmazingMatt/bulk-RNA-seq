# Poster-ready figures for the Stage 01-06 QC evidence
library(edgeR)
library(tidyverse)

project_dir <- "."
config_path <- file.path(project_dir, "config/plot_settings.R")
if (!file.exists(config_path)) {
  stop("Copy config/plot_settings.example.R to config/plot_settings.R and supply your QC observations and decisions before running Stage 14.")
}
source(config_path, local = TRUE)
if (!exists("plot_settings", inherits = FALSE) || !is.list(plot_settings) ||
    !is.list(plot_settings$qc)) {
  stop("config/plot_settings.R must define plot_settings$qc; see the example file.")
}
qc_settings <- plot_settings$qc

require_text <- function(value, field) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
    stop("Supply a nonempty value for plot_settings$qc$", field, ".")
  }
}
require_count <- function(value, field, positive = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0 || value != floor(value) ||
      (positive && value == 0)) {
    stop("Supply a ", if (positive) "positive" else "nonnegative",
         " whole number for plot_settings$qc$", field, ".")
  }
}

sample_ids <- c("groupA_rep1", "groupA_rep2", "groupA_rep3",
                "groupB_rep1", "groupB_rep2", "groupB_rep3")
sample_labels <- setNames(c("Group A 1", "Group A 2", "Group A 3",
                            "Group B 1", "Group B 2", "Group B 3"), sample_ids)
require_text(qc_settings$focus_sample, "focus_sample")
if (!qc_settings$focus_sample %in% sample_ids) {
  stop("plot_settings$qc$focus_sample must match one of the six generic sample IDs.")
}
require_text(qc_settings$complexity_decision, "complexity_decision")
require_text(qc_settings$assignment_strandedness, "assignment_strandedness")

problem_classes <- c("Poor-quality reads", "Low coverage", "Adaptor contamination",
                     "PCR duplicates", "GC-content bias", "Sample mislabelling")
problem_screen <- qc_settings$problem_screen
problem_fields <- c("problem", "status", "evidence", "decision", "status_group")
if (!is.data.frame(problem_screen) || !all(problem_fields %in% names(problem_screen)) ||
    nrow(problem_screen) != length(problem_classes) ||
    !setequal(problem_screen$problem, problem_classes)) {
  stop("Supply one problem_screen row for each of the six classes in config/plot_settings.example.R.")
}
for (field in problem_fields) {
  if (!is.character(problem_screen[[field]]) || anyNA(problem_screen[[field]]) ||
      any(!nzchar(trimws(problem_screen[[field]])))) {
    stop("Complete every problem_screen$", field, " entry using your own evidence.")
  }
}
if (!all(problem_screen$status_group %in% c("candidate", "context", "not_supported"))) {
  stop("Use candidate, context, or not_supported for each problem_screen$status_group.")
}

read_length <- qc_settings$read_length
if (!is.list(read_length)) stop("Supply the read_length list from config/plot_settings.example.R.")
for (field in c("sample_id", "mate", "downstream_performance", "downstream_context",
                "decision", "decision_reason", "unresolved_notes", "limitations")) {
  require_text(read_length[[field]], paste0("read_length$", field))
}
if (!read_length$sample_id %in% sample_ids || !read_length$mate %in% c("R1", "R2")) {
  stop("read_length must identify a generic sample_id and mate R1 or R2.")
}
for (field in c("short_reads", "long_reads", "paired_matches", "missing_mates",
                "malformed_records", "mismatched_mates")) {
  require_count(read_length[[field]], paste0("read_length$", field))
}
for (field in c("short_length_bp", "long_length_bp", "paired_total")) {
  require_count(read_length[[field]], paste0("read_length$", field), positive = TRUE)
}
if (read_length$short_reads + read_length$long_reads == 0 ||
    read_length$short_length_bp >= read_length$long_length_bp ||
    read_length$paired_matches > read_length$paired_total) {
  stop("Check the read-length counts, short/long lengths, and matched-pair denominator.")
}

output_dir <- file.path(project_dir, "results/14_poster_qc_problem_figures_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

navy_colour <- "#0B2A4A"
teal_colour <- "#008A8F"
amber_colour <- "#D99526"
blue_colour <- "#4C78A8"
grey_colour <- "#C7CED6"

poster_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 19, face = "bold", color = navy_colour),
    plot.subtitle = element_text(size = 12, color = "#52677B", margin = margin(b = 10)),
    plot.caption = element_text(size = 10.5, color = "#52677B", hjust = 0, margin = margin(t = 10)),
    axis.title = element_text(size = 13, face = "bold", color = navy_colour),
    axis.text = element_text(size = 11, color = navy_colour),
    strip.text = element_text(size = 11.5, face = "bold", color = navy_colour),
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

# Recalculate the duplication and expression-complexity evidence from active inputs
raw_qc <- read.delim(
  file.path(
    project_dir,
    "results/01_fastqc_raw_results/multiqc_report_data/multiqc_general_stats.txt"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

duplication_column <- "FastQC_mqc_generalstats_fastqc_percent_duplicates"
raw_qc$sample_id <- sub(
  "_L001_R[12]_001\\.fastq\\.gz$",
  "",
  raw_qc$Sample
)

duplication_df <- raw_qc %>%
  group_by(sample_id) %>%
  summarise(
    value = mean(.data[[duplication_column]]),
    .groups = "drop"
  )

counts <- read.delim(
  file.path(project_dir, "results/06_feature_counts_results/feature_counts.txt"),
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE
)

counts <- counts[, c("Geneid", colnames(counts)[7:ncol(counts)])]
colnames(counts) <- c("genes", sample_ids)
count_matrix <- as.matrix(counts[, -1])
raw_dge <- DGEList(counts = count_matrix)

top_100_share <- apply(count_matrix, 2, function(sample_counts) {
  ordered_counts <- sort(sample_counts, decreasing = TRUE)
  100 * sum(head(ordered_counts, 100)) / sum(ordered_counts)
})

detected_genes <- colSums(cpm(raw_dge) > 1)

complexity_df <- bind_rows(
  duplication_df %>%
    mutate(metric = "FastQC duplication estimate", display = sprintf("%.1f%%", value)),
  tibble(
    sample_id = sample_ids,
    value = unname(top_100_share[sample_ids]),
    metric = "Counts in the top 100 genes",
    display = sprintf("%.1f%%", value)
  ),
  tibble(
    sample_id = sample_ids,
    value = unname(detected_genes[sample_ids]),
    metric = "Genes with raw CPM > 1",
    display = format(round(value), big.mark = ",", scientific = FALSE)
  )
) %>%
  mutate(
    sample = factor(sample_labels[sample_id], levels = rev(unname(sample_labels))),
    metric = factor(
      metric,
      levels = c(
        "FastQC duplication estimate",
        "Counts in the top 100 genes",
        "Genes with raw CPM > 1"
      )
    ),
    focus = ifelse(sample_id == qc_settings$focus_sample, "Selected sample", "Other samples")
  )

duplication_values <- sort(duplication_df$value, decreasing = TRUE)
duplication_gap <- duplication_values[1] - duplication_values[2]

complexity_plot <- ggplot(complexity_df, aes(x = value, y = sample, color = focus)) +
  geom_segment(
    aes(x = 0, xend = value, yend = sample),
    linewidth = 0.9,
    alpha = 0.45
  ) +
  geom_point(size = 4.2) +
  geom_text(
    aes(label = display),
    hjust = -0.18,
    size = 3.8,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("Selected sample" = amber_colour, "Other samples" = navy_colour)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(
    title = "Sequence duplication and expression-complexity diagnostics",
    subtitle = paste0(
      sample_labels[[qc_settings$focus_sample]], " is amber. The two highest FastQC estimates differ by ",
      sprintf("%.1f", duplication_gap),
      " percentage points."
    ),
    x = NULL,
    y = NULL,
    color = NULL,
    caption = stringr::str_wrap(
      paste0(
        "FastQC estimates repeated sequences; it cannot separate PCR copies from genuinely abundant RNA. ",
        "Top-gene share and detected-gene breadth are supporting proxies, not UMI-based proof. ",
        "Decision: ", qc_settings$complexity_decision
      ),
      width = 145
    )
  ) +
  poster_theme +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(1.25, "lines")
  )

save_poster_plot(
  complexity_plot,
  "Duplication_library_complexity_diagnostic",
  width = 14,
  height = 7.8
)

# Display the six QC problem classes using the user's evidence and decisions.
problem_screen <- problem_screen %>%
  mutate(
    row = rev(seq_len(n())),
    row_fill = ifelse(row %% 2 == 0, "#F7F9FB", "#FFFFFF"),
    status_fill = recode(
      status_group,
      "candidate" = amber_colour,
      "context" = blue_colour,
      "not_supported" = teal_colour
    )
  )

problem_screen_plot <- ggplot(problem_screen) +
  geom_rect(
    aes(xmin = 0, xmax = 1, ymin = row - 0.43, ymax = row + 0.43, fill = row_fill),
    color = NA
  ) +
  geom_rect(
    aes(xmin = 0, xmax = 0.012, ymin = row - 0.43, ymax = row + 0.43, fill = status_fill),
    color = NA
  ) +
  geom_text(aes(x = 0.025, y = row, label = problem), hjust = 0, size = 4.15, fontface = "bold", color = navy_colour) +
  geom_text(aes(x = 0.245, y = row, label = status, color = status_fill), hjust = 0, size = 3.9, fontface = "bold") +
  geom_text(aes(x = 0.42, y = row, label = stringr::str_wrap(evidence, 50)), hjust = 0, size = 3.65, lineheight = 0.95, color = navy_colour) +
  geom_text(aes(x = 0.82, y = row, label = stringr::str_wrap(decision, 24)), hjust = 0, size = 3.65, lineheight = 0.95, color = navy_colour) +
  annotate("text", x = 0.025, y = 6.65, label = "Problem class", hjust = 0, size = 3.8, fontface = "bold", color = "#52677B") +
  annotate("text", x = 0.245, y = 6.65, label = "Status", hjust = 0, size = 3.8, fontface = "bold", color = "#52677B") +
  annotate("text", x = 0.42, y = 6.65, label = "What the data show", hjust = 0, size = 3.8, fontface = "bold", color = "#52677B") +
  annotate("text", x = 0.82, y = 6.65, label = "Decision", hjust = 0, size = 3.8, fontface = "bold", color = "#52677B") +
  scale_fill_identity() +
  scale_color_identity() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 6.85), clip = "off") +
  labs(
    title = "QC problem-class evidence screen",
    subtitle = "Six problem classes; text and colour show the supplied status for each class.",
    caption = "Statuses, evidence and decisions are supplied in the local configuration; screening alone does not prove a mechanism."
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(size = 19, face = "bold", color = navy_colour),
    plot.subtitle = element_text(size = 12, color = "#52677B", margin = margin(b = 12)),
    plot.caption = element_text(size = 10.5, color = "#52677B", hjust = 0, margin = margin(t = 10)),
    plot.margin = margin(14, 18, 14, 18)
  )

save_poster_plot(
  problem_screen_plot,
  "Manual_problem_class_screen",
  width = 15,
  height = 8.5
)

# Show the supplied read-length observation, pairing check and downstream decision.
short_reads <- read_length$short_reads
long_reads <- read_length$long_reads
total_reads <- short_reads + long_reads
short_percent <- 100 * short_reads / total_reads

read_length_plot <- ggplot() +
  annotate("rect", xmin = 0.08, xmax = 0.95, ymin = 0.25, ymax = 0.88, fill = "#F7F9FB", color = "#D7E0E7", linewidth = 0.8) +
  annotate("rect", xmin = 1.08, xmax = 1.95, ymin = 0.25, ymax = 0.88, fill = "#F7F9FB", color = "#D7E0E7", linewidth = 0.8) +
  annotate("rect", xmin = 2.08, xmax = 2.95, ymin = 0.25, ymax = 0.88, fill = "#F7F9FB", color = "#D7E0E7", linewidth = 0.8) +
  annotate("rect", xmin = 3.08, xmax = 3.95, ymin = 0.25, ymax = 0.88, fill = "#E8F3F1", color = teal_colour, linewidth = 1.1) +
  annotate("segment", x = 0.96, xend = 1.06, y = 0.565, yend = 0.565, arrow = arrow(length = unit(0.14, "inches")), color = "#7A8A99", linewidth = 0.8) +
  annotate("segment", x = 1.96, xend = 2.06, y = 0.565, yend = 0.565, arrow = arrow(length = unit(0.14, "inches")), color = "#7A8A99", linewidth = 0.8) +
  annotate("segment", x = 2.96, xend = 3.06, y = 0.565, yend = 0.565, arrow = arrow(length = unit(0.14, "inches")), color = "#7A8A99", linewidth = 0.8) +
  annotate("text", x = 0.515, y = 0.80, label = "1  Raw observation", size = 4.4, fontface = "bold", color = navy_colour) +
  annotate("text", x = 0.515, y = 0.61, label = paste(sample_labels[[read_length$sample_id]], read_length$mate), size = 4.2, fontface = "bold", color = navy_colour) +
  annotate("text", x = 0.515, y = 0.49, label = paste0(format(short_reads, big.mark = ","), " / ", format(total_reads, big.mark = ","), " reads\nwere ", read_length$short_length_bp, " bp (", sprintf("%.3f", short_percent), "%)"), size = 3.9, lineheight = 1.05, color = navy_colour) +
  annotate("rect", xmin = 0.20, xmax = 0.83, ymin = 0.34, ymax = 0.39, fill = navy_colour, color = NA) +
  annotate("rect", xmin = 0.20, xmax = 0.20 + 0.63 * short_percent / 100, ymin = 0.34, ymax = 0.39, fill = amber_colour, color = NA) +
  annotate("text", x = 1.515, y = 0.80, label = "2  Pairing check", size = 4.4, fontface = "bold", color = navy_colour) +
  annotate("text", x = 1.515, y = 0.58, label = paste0(format(read_length$paired_matches, big.mark = ","), " / ", format(read_length$paired_total, big.mark = ","), "\nread IDs matched"), size = 4.1, fontface = "bold", lineheight = 1.05, color = navy_colour) +
  annotate("text", x = 1.515, y = 0.41, label = paste0(read_length$missing_mates, " missing, ", read_length$malformed_records, " malformed,\n", read_length$mismatched_mates, " mismatched mates"), size = 3.8, lineheight = 1.05, color = navy_colour) +
  annotate("text", x = 2.515, y = 0.80, label = "3  Downstream performance", size = 4.4, fontface = "bold", color = navy_colour) +
  annotate("text", x = 2.515, y = 0.57, label = read_length$downstream_performance, size = 4.0, fontface = "bold", lineheight = 1.08, color = navy_colour) +
  annotate("text", x = 2.515, y = 0.36, label = read_length$downstream_context, size = 3.8, color = navy_colour) +
  annotate("text", x = 3.515, y = 0.80, label = "4  Decision", size = 4.4, fontface = "bold", color = teal_colour) +
  annotate("text", x = 3.515, y = 0.60, label = read_length$decision, size = 4.7, fontface = "bold", color = navy_colour) +
  annotate("text", x = 3.515, y = 0.44, label = read_length$decision_reason, size = 3.9, lineheight = 1.05, color = navy_colour) +
  annotate("text", x = 3.515, y = 0.31, label = read_length$unresolved_notes, size = 3.5, color = "#52677B") +
  coord_cartesian(xlim = c(0, 4.03), ylim = c(0.18, 0.94), clip = "off") +
  labs(
    title = paste0(sample_labels[[read_length$sample_id]], ": read-length case study (", read_length$short_length_bp, " and ", read_length$long_length_bp, " bp)"),
    subtitle = "Observation → structural check → downstream performance → evidence-based decision",
    caption = read_length$limitations
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(size = 19, face = "bold", color = navy_colour),
    plot.subtitle = element_text(size = 12, color = "#52677B", margin = margin(b = 10)),
    plot.caption = element_text(size = 10.5, color = "#52677B", hjust = 0, margin = margin(t = 10)),
    plot.margin = margin(14, 18, 14, 18)
  )

save_poster_plot(
  read_length_plot,
  "Read_length_impact_chain",
  width = 14,
  height = 6.5
)

# Redesign the Stage 06 assignment output using short sample labels and percentages
assignment_df <- read.delim(
  file.path(
    project_dir,
    "results/06_feature_counts_results/feature_counts_summary_data/featureCounts_assignment_plot.txt"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

assignment_long <- assignment_df %>%
  pivot_longer(-Sample, names_to = "category", values_to = "fragments") %>%
  group_by(Sample) %>%
  mutate(percent = 100 * fragments / sum(fragments)) %>%
  ungroup() %>%
  mutate(
    sample = factor(sample_labels[Sample], levels = unname(sample_labels)),
    category = recode(
      category,
      "Assigned" = "Assigned",
      "Unassigned: Multi Mapping" = "Multi-mapping",
      "Unassigned: No Features" = "No annotated feature",
      "Unassigned: Ambiguity" = "Ambiguous feature",
      "Unassigned: Unmapped" = "Unmapped"
    ),
    category = factor(
      category,
      levels = c("Unmapped", "Ambiguous feature", "No annotated feature", "Multi-mapping", "Assigned")
    )
  )

assigned_labels <- assignment_long %>%
  filter(category == "Assigned") %>%
  mutate(label = sprintf("%.1f%% assigned", percent))

assignment_plot <- ggplot(assignment_long, aes(x = sample, y = percent, fill = category)) +
  geom_col(width = 0.72) +
  geom_text(
    data = assigned_labels,
    mapping = aes(x = sample, y = 35, label = label),
    color = "white",
    size = 3.9,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c(
    "Assigned" = teal_colour,
    "Multi-mapping" = blue_colour,
    "No annotated feature" = "#D8B365",
    "Ambiguous feature" = "#9E9AC8",
    "Unmapped" = grey_colour
  )) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0)), labels = function(x) paste0(x, "%")) +
  labs(
    title = "featureCounts assignment across libraries",
    subtitle = "featureCounts assignment categories shown as a percentage of assessed fragments",
    x = NULL,
    y = "Fragments assessed",
    fill = "Assignment category",
    caption = paste0(
      "Assigned rates span ", sprintf("%.1f", min(assigned_labels$percent)),
      "–", sprintf("%.1f", max(assigned_labels$percent)), "% of assessed fragments. ",
      "Counting configuration: ", qc_settings$assignment_strandedness
    )
  ) +
  poster_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

save_poster_plot(
  assignment_plot,
  "FeatureCounts_assignment_consistency",
  width = 11,
  height = 7.5
)

cat("Poster QC/problem figures written to:", output_dir, "\n")
print(sort(list.files(output_dir)))
