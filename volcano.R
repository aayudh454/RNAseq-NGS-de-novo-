setwd("C:/Users/aayudh.das/OneDrive - Garuda Therapeutics/RNA-seq")
list.files()
library(dplyr)
library(ggplot2)
library(ggbreak)

data <- read.csv("GMP1 vs ICB0004-03.csv")
head(data)

# Add a column to classify significance
data <- data %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 2 ~ "Up",
    padj < 0.05 & log2FoldChange < -2 ~ "Down",
    TRUE ~ "Not Significant"
  ))


# Final volcano plot with Arial font, no legend, and x-axis ticks every 2 units
ggplot(data, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "skyblue", "Down" = "skyblue", "Not Significant" = "grey40")) +
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
    legend.position = "none"
  )

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)


# Define the highlight genes
highlight_genes <- c("RUNX1", "HOXA9", "MLLT3", "MECOM", "HLF", "SPINK2", "PROM1",
                     "HLA-DRA", "ALDH1A1", "KCNK17", "IL33", "CDH5", "ITGA2B",
                     "IL3RA", "CSF1R")

# Classify genes
data <- data %>%
  mutate(GeneGroup = ifelse(name %in% highlight_genes, "Highlighted", Significance))

# Split data
highlighted_data <- subset(data, GeneGroup == "Highlighted")
non_highlighted_data <- subset(data, GeneGroup != "Highlighted")

# Plot
png("Nature_2022_genelist.png", width = 7, height = 7, units = 'in', res = 300)
ggplot() +
  geom_point(data = non_highlighted_data, 
             aes(x = log2FoldChange, y = -log10(pvalue), color = GeneGroup),
             alpha = 0.6, size = 1.5) +
  geom_point(data = highlighted_data, 
             aes(x = log2FoldChange, y = -log10(pvalue)),
             color = "red", size = 3) +  # Slightly bigger red points
  geom_text_repel(data = highlighted_data,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = name),
                  size = 4, family = "Arial", max.overlaps = 100) +
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
