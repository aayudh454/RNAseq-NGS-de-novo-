setwd("C:/Users/aayudh.das/OneDrive - StratusTX/RNA-seq")
list.files()
library(tidyverse)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(ggbreak)

library(tidyverse)

data <- read.csv("piezo_dotplot.csv")
head(data)

# Compute means
data_long <- data %>%
  mutate(
    mean_HSC = rowMeans(select(., T3_R1, T3_R2, T3_R3)),
    mean_iPSC = rowMeans(select(., T2_R1, T2_R2, T2_R3))
  ) %>%
  select(Name, mean_HSC, mean_iPSC) %>%
  pivot_longer(cols = c(mean_iPSC, mean_HSC), names_to = "Group", values_to = "MeanCount") %>%
  mutate(Group = factor(Group, levels = c("mean_iPSC", "mean_HSC"), labels = c("iPSC", "HSC")))

# Define gene-specific colors
sample_colors <- c(
  "TRPV4" = "#4BC9FF",
  "PIEZO2" = "#404444",
  "PIEZO1" = "#D50032"
)

# Generate dot plot
png("peizo_dotplot.png", width = 4, height = 4, units = 'in', res = 300)
ggplot(data_long, aes(x = Group, y = MeanCount, group = Name, color = Name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_line(alpha = 0.3) +
  scale_color_manual(values = sample_colors) +
  scale_y_continuous(breaks = pretty(data_long$MeanCount, n = 5)) +
  labs(
    x = NULL,
    y = "Mean Normalized Gene Counts",
    color = "Gene"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right",
    plot.title = element_blank()
  )
dev.off()

# Optional: order genes by iPSC mean for better visual grouping
gene_order <- data_long %>%
  filter(Group == "iPSC") %>%
  arrange(desc(MeanCount)) %>%
  pull(Name)

data_long$Name <- factor(data_long$Name, levels = gene_order)

# Heatmap
png("peizo_heatmap.png", width = 4, height = 4, units = 'in', res = 300)
ggplot(data_long, aes(x = Group, y = Name, fill = MeanCount)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(MeanCount, 0)), size = 3) +
  scale_fill_gradient(low = "white", high = "#D50032") +
  labs(
    x = NULL,
    y = NULL,
    fill = "Mean Count"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black")
  )
dev.off()


####VOLCANO
data <- read.csv("GMP1 vs ICB0004-03.csv")
data <- data %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "Not Significant"
  ))

sample_colors <- c(
  "TRPV4" = "#4BC9FF",
  "PIEZO2" = "#404444",
  "PIEZO1" = "#D50032"
)
# Define specific HLA genes to highlight
highlight_hla_genes <- c("PIEZO1", "PIEZO2", "TRPV4")

# Classify genes
data <- data %>%
  mutate(GeneGroup = ifelse(name %in% highlight_hla_genes, "Piezo_Highlight", Significance))

# Split data
Piezo_Highlighted_data <- subset(data, GeneGroup == "Piezo_Highlight")
non_highlighted_data <- subset(data, GeneGroup != "Piezo_Highlight")

# Plot
png("Piezo_volcano.png", width = 4, height = 4, units = 'in', res = 300)
ggplot() +
  geom_point(data = non_highlighted_data, 
             aes(x = log2FoldChange, y = -log10(pvalue), color = GeneGroup),
             alpha = 0.6, size = 1.5) +
  geom_point(data = Piezo_Highlighted_data, 
             aes(x = log2FoldChange, y = -log10(pvalue)),
             color = "#D50032", size = 3) +  # Red dots size 3
  geom_text_repel(data = Piezo_Highlighted_data,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = name),
                  size = 4, family = "Arial", max.overlaps = 100) +  # Label size 3
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