setwd("C:/Users/aayudh.das/OneDrive - Garuda Therapeutics/RNA-seq")
list.files()
library(dplyr)
library(ggplot2)
library(ggbreak)
library(ggplot2)
library(dplyr)
library(ggrepel)

data <- read.csv("GMP1 vs ICB0004-03.csv")
head(data)

data <- data %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 2 ~ "Up",
    padj < 0.05 & log2FoldChange < -2 ~ "Down",
    TRUE ~ "Not Significant"
  ))


# Define specific HLA genes to highlight
highlight_hla_genes <- c("HLA-A", "HLA-DPB1", "HLA-DQB1", "HLA-DRB1", "HLA-B", "HLA-C")

# Classify genes
data <- data %>%
  mutate(GeneGroup = ifelse(name %in% highlight_hla_genes, "HLA_Highlight", Significance))

# Split data
hla_highlighted_data <- subset(data, GeneGroup == "HLA_Highlight")
non_highlighted_data <- subset(data, GeneGroup != "HLA_Highlight")

# Plot
png("HLAgenes.png", width = 6.5, height = 7, units = 'in', res = 300)
ggplot() +
  geom_point(data = non_highlighted_data, 
             aes(x = log2FoldChange, y = -log10(pvalue), color = GeneGroup),
             alpha = 0.6, size = 1.5) +
  geom_point(data = hla_highlighted_data, 
             aes(x = log2FoldChange, y = -log10(pvalue)),
             color = "red", size = 3) +  # Red dots size 3
  geom_text_repel(data = hla_highlighted_data,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = name),
                  size = 4, family = "Arial", max.overlaps = 100) +  # Label size 3
  scale_color_manual(values = c(
    "Up" = "skyblue",
    "Down" = "skyblue",
    "Not Significant" = "grey40"
  )) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
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
