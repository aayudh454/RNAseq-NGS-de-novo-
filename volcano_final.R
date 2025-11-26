setwd("C:/Users/aayudh.das/OneDrive - Garuda Therapeutics/RNA-seq")
list.files()
library(dplyr)
library(ggplot2)
library(ggbreak)

data <- read.csv("GMP1 vs ICB0004-03.csv")
head(data)


# Integrated volcano plot1
# Define gene categories
hematopoietic_markers <- c("RUNX1", "HOXA9", "MLLT3", "MECOM", "HLF", "SPINK2", "PROM1",
                           "HLA-DRA", "ALDH1A1", "KCNK17", "IL33", "CDH5", "ITGA2B",
                           "IL3RA", "CSF1R", "THY1", "PTPRC", "SPN", "CXCR4", "DLL4", "CD34")

transcription_factors <- c("RUNX1", "GATA2", "TAL1")
hox_markers <- c("HOXA9", "HOXA10")

data <- data %>%
  mutate(GeneGroup = case_when(
    name %in% hox_markers ~ "HOX Markers",
    name %in% transcription_factors ~ "Transcription Factors",
    name %in% hematopoietic_markers ~ "Hematopoietic Markers",
    padj < 0.05 & log2FoldChange > 2 ~ "Up",
    padj < 0.05 & log2FoldChange < -2 ~ "Down",
    TRUE ~ "Not Significant"
  ))

# Separate highlighted and other points
highlighted_data <- subset(data, GeneGroup %in% c("Hematopoietic Markers", "Transcription Factors", "HOX Markers"))
non_highlighted_data <- subset(data, !(GeneGroup %in% c("Hematopoietic Markers", "Transcription Factors", "HOX Markers")))

# Plot
png("1.volcano_final.png", width = 8, height = 6, units = 'in', res = 300)
ggplot() +
  geom_point(data = non_highlighted_data,
             aes(x = log2FoldChange, y = -log10(pvalue), color = GeneGroup),
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_point(data = highlighted_data,
             aes(x = log2FoldChange, y = -log10(pvalue), color = GeneGroup),
             size = 3) +
 geom_text_repel(data = highlighted_data,
                aes(x = log2FoldChange, y = -log10(pvalue), label = name),
                size = 4, family = "Arial", fontface = "bold", max.overlaps = 100) +
  scale_color_manual(
    values = c(
      "Hematopoietic Markers" = "red",
      "Transcription Factors" = "green3",
      "HOX Markers" = "purple",
      "Up" = "skyblue",
      "Down" = "skyblue",
      "Not Significant" = "grey40"
    ),
    breaks = c("Hematopoietic Markers", "Transcription Factors", "HOX Markers"), 
    name = ""
  ) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed") +
  scale_x_continuous(breaks = seq(-8, 12, 2), limits = c(-8, 12)) +
  scale_y_continuous(limits = c(0, 310)) +
  labs(x = "Fold Change (Log2)", y = "P-Value (-log10)") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 16),
    axis.title = element_text(family = "Arial", size = 16),
    axis.text = element_text(family = "Arial", size = 14),
    legend.position = "bottom"
  ) +
  coord_cartesian(clip = "off")
dev.off()

# Integrated volcano plot2
# Load data
#setwd("C:/Users/sita.patel/OneDrive - Garuda Therapeutics/Documents/Projects/RNA_Seq/R/visualization/")
data <- read.csv("GMP1 vs ICB0004-03.csv")

# Define gene sets
stress_markers <- c("FOS", "JUN", "EGR1")
wnt_pathway <- c("CTNNB1", "AXIN2")
notch_pathway <- c("NOTCH1", "HES1")
inflammatory_pathway <- c("STAT1", "IRF7")

# Define gene groups
data <- data %>%
  mutate(GeneGroup = case_when(
    name %in% stress_markers ~ "Stress Response",
    name %in% wnt_pathway ~ "Wnt/β-catenin",
    name %in% notch_pathway ~ "Notch",
    name %in% inflammatory_pathway ~ "Inflammatory",
    padj < 0.05 & log2FoldChange > 2 ~ "Up",
    padj < 0.05 & log2FoldChange < -2 ~ "Down",
    TRUE ~ "Not Significant"
  ))

# Separate highlighted and non-highlighted
highlighted_data <- subset(data, GeneGroup %in% c("Stress Response", "Wnt/β-catenin", "Notch", "Inflammatory"))
non_highlighted_data <- subset(data, !(GeneGroup %in% c("Stress Response", "Wnt/β-catenin", "Notch", "Inflammatory")))

# Plot
png("volcano_stress_signalling.png", width = 8, height = 6.5, units = 'in', res = 300)
ggplot() +
  geom_point(data = non_highlighted_data,
             aes(x = log2FoldChange, y = -log10(pvalue), color = GeneGroup),
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_point(data = highlighted_data,
             aes(x = log2FoldChange, y = -log10(pvalue), color = GeneGroup),
             size = 3) +
  geom_text_repel(data = highlighted_data,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = name),
                  size = 4, fontface = "bold", family = "Arial", max.overlaps = 100) +
  scale_color_manual(
    values = c(
      "Stress Response" = "gold",
      "Wnt/β-catenin" = "purple",      
      "Notch" = "#32CD32",              
      "Inflammatory" = "red",        
      "Up" = "skyblue",
      "Down" = "skyblue",
      "Not Significant" = "grey40"
    ),
    breaks = c("Stress Response", "Wnt/β-catenin", "Notch", "Inflammatory"),
    name = ""
  ) +
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
    legend.position = "bottom"
  ) +
  coord_cartesian(clip = "off")
dev.off()