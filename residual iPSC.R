setwd("C:/Users/aayudh.das/OneDrive - StratusTX/RNA-seq")
list.files()
library(tidyverse)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(ggbreak)

data <- read.csv("residual iPSC markers.csv")
data$Name[data$Name == "POU5F1P3"] <- "OCT4"
head(data)



# Compute means
data_long <- data %>%
  mutate(
    mean_HSC = rowMeans(select(., T3_R1, T3_R2, T3_R3)),
    mean_iPSC = rowMeans(select(., T2_R1, T2_R2, T2_R3))
  ) %>%
  select(Name, mean_HSC, mean_iPSC) %>%
  pivot_longer(cols = c(mean_iPSC, mean_HSC), names_to = "Group", values_to = "MeanCount") %>%
  mutate(Group = factor(Group, levels = c("mean_iPSC", "mean_HSC"), labels = c("iPSC", "ST-101")))

png("residualpscvsHSC_dotplot.png", width = 6, height = 6, units = 'in', res = 300)
ggplot(data_long, aes(x = Group, y = MeanCount, group = Name, color = Name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_line(alpha = 0.3) +
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
ggplot(data_long, aes(x = Group, y = Name, fill = MeanCount)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(MeanCount, 0)), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Group",
    y = "Gene",
    fill = "Mean Count"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black")
  )


# Compute means and reshape
data_long <- data %>%
  filter(Name != "GAPDH") %>%
  mutate(
    mean_HSC = rowMeans(select(., T3_R1, T3_R2, T3_R3)),
    mean_iPSC = rowMeans(select(., T2_R1, T2_R2, T2_R3))
  ) %>%
  select(Name, mean_HSC, mean_iPSC) %>%
  pivot_longer(cols = c(mean_iPSC, mean_HSC), names_to = "Group", values_to = "MeanCount") %>%
  mutate(
    Group = factor(Group, levels = c("mean_iPSC", "mean_HSC"), labels = c("iPSC", "ST-101")),
    DisplayLabel = ifelse(MeanCount < 50, "0", as.character(round(MeanCount, 0)))
  )

# Order genes by iPSC expression
gene_order <- data_long %>%
  filter(Group == "iPSC") %>%
  arrange(desc(MeanCount)) %>%
  pull(Name)

data_long$Name <- factor(data_long$Name, levels = gene_order)

# Plot heatmap
setwd("C:/Users/aayudh.das/OneDrive - StratusTX/RNA-seq")
png("residualpscvsHSC_heatmap.png", width = 4, height = 4, units = 'in', res = 300)
ggplot(data_long, aes(x = Group, y = Name, fill = MeanCount)) +
  geom_tile(color = "white") +
  geom_text(aes(label = DisplayLabel), size = 3, family = "Arial") +
  scale_fill_gradient(low = "white", high = "#4BC9FF", name = "Mean Gene Counts") +
  labs(
    x = NULL,
    y = "Gene"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black", size = 14, family = "Arial", face = "bold"),
    axis.text.y = element_text(color = "black", size = 14, family = "Arial", face = "bold"),
    axis.title = element_text(color = "black", size = 16, family = "Arial", face = "bold")
  )

dev.off()


setwd("C:/Users/aayudh.das/OneDrive - StratusTX/RNA-seq")
png("residualpscvsHSC_heatmap.png", width = 4, height = 4, units = 'in', res = 300)
ggplot(data_long, aes(x = Group, y = Name, fill = MeanCount)) +
  geom_tile(color = "white") +
  geom_text(aes(label = DisplayLabel), size = 3, family = "Arial") +
  scale_fill_gradient(low = "white", high = "#4BC9FF", name = "Mean Gene Counts") +
  labs(
    x = NULL,
    y = "Gene"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black", size = 14, family = "Arial", face = "bold"),
    axis.text.y = element_text(color = "black", size = 14, family = "Arial", face = "bold"),
    axis.title = element_text(color = "black", size = 16, family = "Arial", face = "bold"),
    legend.position = "none"
  )
dev.off()