# Create a data frame in R with the bold sample names
data <- data.frame(
  Sample = c("ENG1", "PD1", "PD2", "PD3", 
             "PD6", "PD7", "PD8"),
  Viability = c(77.8, 86.1, 88, 93.9, 71, 63.1, 70.3),
  CD38neg = c(82, 87.7, 93, 95, 98.4, 96.6, 90.7),
  CD34pos = c(87.4, 91.2, 89.6, 91, 76.4, 65.9, 90.7),
  CD34_CD90_CD45neg_HSC = c(5.42, 30.6, 7.8, 12, 36, 8.41, 11.5),
  CD34_CD90neg_CD45pos_MLP = c(1.65, 45.4, 75.4, 43.9, 34.3, 2.58, 2.52),
  CD34_CD90neg_CD45neg_MPP = c(12.1, 2.43, 10.6, 40.4, 23, 22.4, 39.7),
  CD34_CD90pos_CD45pos = c(1.65, 9.51, 6.6, 3.4, 19.3, 19, 37),
  CD34_CD90pos_CD45neg = c(9.39, 8.79, 1.87, 3.5, 0.48, 2.24, 0.62),
  CD34_CD90neg_CD45neg = c(0.98, 0.35, 1.59, 1.6, 0.38, 1, 1.84)
)

# View the data
print(data)

# Set row names and remove Sample column for PCA
rownames(data) <- data$Sample
data_matrix <- data[, -1]

# Perform PCA
pca_result <- prcomp(data_matrix, scale. = TRUE)

# Create a data frame for plotting
pca_df <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  PC3 = pca_result$x[,3]
)

# Plot PCA with labels
setwd("C:/Users/aayudh.das/OneDrive - StratusTX/RNA-seq")
png("HSC_PCA.png", width = 5, height = 5, units = 'in', res = 300)
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(color = "#D50032", size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5) +
  ggtitle("PCA of HSC Panel Data") +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text = element_text(size = 16, family = "Arial"),
    axis.title = element_text(size = 16, family = "Arial"),
    plot.title = element_text(size = 18, family = "Arial", hjust = 0.5)
  )
dev.off()