setwd("C:/Users/aayudh.das/OneDrive - StratusTX/RNA-seq")
list.files()
library(dplyr)
library(ggplot2)
library(ggbreak)
library(ggrepel)

####VOLCANO
data <- read.csv("GMP1 vs ICB0004-03.csv")
data <- data %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "Not Significant"
  ))

sample_colors <- c(
  "CD7" = "#4BC9FF",
  "CD4" = "#404444",
  "CD8A" = "#D50032",
  "CD8B" = "#D50032"
)
# Define specific HLA genes to highlight
highlight_Tcell_genes <- c("CD7", "CD4", "CD8A","CD8B")

# Classify genes
data <- data %>%
  mutate(GeneGroup = ifelse(name %in% highlight_Tcell_genes, "Tcell_genes_Highlight", Significance))

# Split data
Tcell_genes_Highlighted_data <- subset(data, GeneGroup == "Tcell_genes_Highlight")
non_highlighted_data <- subset(data, GeneGroup != "Tcell_genes_Highlight")

# Plot
png("Tcell_genes.png", width = 8, height = 6, units = 'in', res = 300)
ggplot() +
  geom_point(data = non_highlighted_data, 
             aes(x = log2FoldChange, y = -log10(pvalue), color = GeneGroup),
             alpha = 0.6, size = 1.5) +
  geom_point(data = Tcell_genes_Highlighted_data, 
             aes(x = log2FoldChange, y = -log10(pvalue)),
             color = "#D50032", size = 4) +  # Red dots size 3
  geom_text_repel(data = Tcell_genes_Highlighted_data,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = name),
                  size = 4.5, family = "Arial", max.overlaps = 100) +  # Label size 3
  scale_color_manual(values = c(
    "Up" = "#4BC9FF",
    "Down" = "#4BC9FF",
    "Not Significant" = "#404444"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_x_continuous(breaks = seq(-8, 10, 2), limits = c(-8, 10)) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(x = "Fold Change (Log2)", y = "P-Value (-log10)") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 16),
    axis.title = element_text(family = "Arial", size = 16),
    axis.text = element_text(family = "Arial", size = 14),
    legend.position = "none")
dev.off()