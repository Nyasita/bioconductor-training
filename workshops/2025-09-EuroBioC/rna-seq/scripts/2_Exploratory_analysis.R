##########################################
# Exploratory analysis and quality control
##########################################

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(vsn)
  library(ggplot2)
  library(ComplexHeatmap)
  library(hexbin)
  library(iSEE)
})

library("RColorBrewer")

se <- readRDS("data/GSE96870_se.rds")
nrow(se)

table(rowData(se)$gbkey)

# Keep only the genes (mRNA)
se_mRNA <- se[rowData(se)$gbkey == "mRNA", ]
nrow(se_mRNA)

# Remove rows that do not have > 5 total counts
se <- se[rowSums(assay(se)) > 5, ]

nrow(se)

####### Challenge 1:
# 1. How many of each feature remained?
table(rowData(se)$gbkey)

# 2. Choose difference filtering?
length(which(rowSums(assay(se)) > 10))
length(which(rowSums(assay(se)) > 20))
# 3. Ask pros/cons
# Discussion

# Add in the sum of all counts
se$libSize <- colSums(assay(se))

# Plot the libSize by using the native pipe in R |>
# to extract the colData, turn it to a dataframe
# and pass it to ggplot:

colData(se) |>
  as.data.frame() |>
  ggplot(aes(x = Label, y = libSize / 1e6, fill = Group)) +
    geom_bar(stat = "identity") + theme_bw() +
    labs(x = "Sample", y = "Total count in millions") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Convert the se to a deseq2 object
dds <- DESeq2::DESeqDataSet(se, design = ~sex + time)

# Calculate the size factors for making our libraries comparable
dds <- estimateSizeFactors(dds)

# Plot the size factors against library size
# and look for any patterns by group
ggplot(data.frame(libSize = colSums(assay(dds)),
                  sizeFactor = sizeFactors(dds),
                  Group = dds$Group),
       aes(x = libSize, y = sizeFactor, col = Group)) +
  geom_point(size = 5) + theme_bw() +
  labs(x = "Library size", y = "Size Factor")

# Plot mean expression vs. variance

meanSdPlot(assay(dds), ranks = FALSE)

# Transform the raw data in order to stabilize the variance (mean - variance independent)
vsd <- DESeq2::vst(dds, blind = TRUE)

meanSdPlot(assay(vsd), ranks = FALSE)

# Construct a Heatmap
# Using the Euclidean distance we can see how close or far away are our samples each other. 

# Calculate the distance (Euclidean). We need to transpose the normalised matrix
dst <- dist(t(assay(vsd)))

# Create a colour panel with blue
colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)

ComplexHeatmap::Heatmap(
  as.matrix(dst),
  col = colors,
  name = "Euclidean\ndistance",
  cluster_rows  = hclust(dst),
  cluster_columns = hclust(dst),
  bottom_annotation = columnAnnotation(
    sex = vsd$sex,
    time = vsd$time,
    col = list(sex = c(Female = "red", Male = "blue"),
               time = c(Day0 = "yellow", Day4 = "forestgreen", Day8 = "purple")))
)

# PCA: Reveals the main drivers of variation of the dataset
# It is unsupervised, doesn't take account the labels (e.g., conditions)
# It help us to understand if the variation in our data derives from a true biological event or technicalities

# Create a pca object
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("sex", "time"),
                           returnData = TRUE)

# calc the PCs percentages
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot the PCA with ggplot
ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sex, shape = time), size = 5) +
  theme_minimal() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_manual(values = c(Male = "blue", Female = "red"))
  

######### Challenge 2
# Making the PCA plot of the challenge 2.2
# 1) Samples = Group
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("Group","time"), returnData = TRUE)

# 2) Plot: color by sample, shape by time
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Group, shape = time), size = 5) +
  theme_minimal() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_brewer(palette = "Set2", name = "Group")   # nice distinct colors

###### Challenge 3.
# Compare before and after VST (Variance Stabilizing Transformation)

# After
pcaDataVst <- DESeq2::plotPCA(vsd, intgroup = c("libSize"),
                              returnData = TRUE)

percentVar <- round(100 * attr(pcaDataVst, "percentVar"))

ggplot(pcaDataVst, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = libSize / 1e6), size = 5) +
  theme_minimal() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_continuous("Total count in millions", type = "viridis")

# Before
pcaDataCts <- DESeq2::plotPCA(DESeqTransform(se), intgroup = c("libSize"),
                              returnData = TRUE)

ggplot(pcaDataCts, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = libSize / 1e6), size = 5) +
  theme_minimal() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_continuous("Total count in millions", type = "viridis")

#### EXTRA
# Running PCA directly on the VST data
mat <- assay(vsd)         # expression matrix (genes x samples)
mat <- t(mat)             # transpose → samples x genes

pca <- prcomp(mat, center = TRUE, scale. = FALSE)

# pca$x = PC scores (samples × PCs)
# Each row = a sample, each column = PC1, PC2, …
# Example: use in regression or clustering.
# pca$rotation = PC loadings (genes × PCs)
# Each row = a gene, each column = weight of that gene in that PC.
# Example: see which genes contribute most to PC1.
# pca$sdev = standard deviation of each PC
# Square them and divide by total variance to get % variance explained.

# Top genes driving PC1
head(sort(abs(pca$rotation[,1]), decreasing = TRUE), 10)

# PC3 vs.PC4
# Run PCA on the most variable genes (like plotPCA does)
ntop <- 500
topVarGenes <- head(order(rowVars(mat), decreasing = TRUE), ntop)
pca <- prcomp(mat[, topVarGenes], center = TRUE, scale. = FALSE)

# Inspect results
head(pca$x)         # sample scores (PC1, PC2, PC3, ...)
head(pca$rotation)  # gene loadings

ggplot(as.data.frame(pca$x), aes(x = PC3, y = PC4, color = vsd$sex, shape = vsd$time)) +
  geom_point(size = 5)

library(matrixStats)
mat <- assay(vsd)

ntop <- 500
topVarGenes <- head(order(rowVars(mat), decreasing = TRUE), ntop)

pca <- prcomp(mat[topVarGenes, ], center = TRUE, scale. = FALSE)  
pca$rotation  

loadings <- pca$x

top_pc1 <- names(sort(abs(loadings[,1]), decreasing = TRUE))[1:30]  
top_pc1

scores <- as.data.frame(pca$rotation)


