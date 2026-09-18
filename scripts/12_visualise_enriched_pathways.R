#Load Pre-installed Libraries
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(GO.db)
library(AnnotationDbi)
library(tidyverse)


#prepare the CAMERA results
sig_camera <- read.csv("results/11_camera_pathway_results/CAMERA_Competitive_Enrichment_groupA_vs_groupB.csv",
                       row.names = 1, check.names = FALSE)


#Visualise the top pathways
#plots


# CAMERA DOT PLOT
camera_plot <- sig_camera %>%
  mutate(
    GO_ID = rownames(sig_camera),
    negLogFDR =
      
      -log10(FDR)
    
  ) %>%
  arrange(FDR) %>%
  slice_head(n = 20)
# Convert GO IDs to readable GO terms
camera_plot$TERM <- Term(
    GOTERM[camera_plot$GO_ID] )
# If a GO term cannot be found, keep the GO ID
camera_plot$TERM <- ifelse(
    is.na(camera_plot$TERM),
    camera_plot$GO_ID,
    camera_plot$TERM )
# Order pathways by significance
camera_plot$TERM <- factor(
    camera_plot$TERM,
    levels = rev(camera_plot$TERM) )
# Plot
ggplot(
  camera_plot,
  aes(
    x = negLogFDR,
    y = TERM,
    size = NGenes,
    shape = Direction
  )
) +
  geom_point() +
  labs(
    title = "CAMERA Gene Set Enrichment",
    subtitle = "Group A vs Group B",
    x = "
-log10(FDR)",
    y = "GO Biological Process",
    size = "Number of genes",
    shape = "Direction"
  ) +
  

theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )
# Save plot
ggsave(
  "results/12_pathway_visualisation_results/CAMERA_dotplot_groupA_vs_groupB.png",
  width = 10,
  height = 8,
  dpi = 300
)

