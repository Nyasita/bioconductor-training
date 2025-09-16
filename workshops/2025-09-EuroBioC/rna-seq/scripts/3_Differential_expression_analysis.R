##################################
# Differential expression analysis
##################################

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(ExploreModelMatrix)
  library(cowplot)
  library(ComplexHeatmap)
  library(apeglm)
})

# Load the data 
se <- readRDS("data/GSE96870_se.rds")

# Filter out genes/features with less than 5 counts in total
se <- se[rowSums(assay(se)) > 5, ]

# Convert the summarised experiment oblext into a DESeq2 oblect
dds <- DESeq2::DESeqDataSet(se, design = ~sex + time)

# Normalization
# Estimate the size factors to make the samples comparable
dds <- estimateSizeFactors(dds)

# Estimate the wise-gene dispersion to calculate the variance withing the groups
dds <- estimateDispersions(dds)

plotDispEsts(dds)

# Testing (WaldTest on a generalised linear model)
dds <- nbinomWaldTest(dds)

# All in one command
dds <- DESeq(dds)

# Explore results
# Check the results names
resultsNames(dds)

# Day8 vs Day0
resTime <- results(dds, contrast = c("time", "Day8", "Day0"))
summary(resTime)

#View(resTime)
head(resTime[order(resTime$pvalue), ])

# Independent filtering
# Plot an MA diagram
plotMA(resTime, alpha = 0.05)

# Shring the genes with low mean and high dispersion. They include low level information

resTimeLfc <- lfcShrink(dds, coef = "time_Day8_vs_Day0", res = resTime)
plotMA(resTimeLfc, alpha = 0.05)

# Heatmap visualization

# Transform counts
vsd <- vst(dds, blind = TRUE)

# Get top DEGs
genes <- resTime[order(resTime$pvalue), ] |>
  head(10) |>
  rownames()

heatmapData <- assay(vsd)[genes,]

# Scale counts for visualization
heatmapData <- t(scale(t(heatmapData)))

# Add annotation
heatmapColAnnot <- data.frame(colData(vsd)[, c("time", "sex")])
heatmapColAnnot <- HeatmapAnnotation(df = heatmapColAnnot)

# Plot as heatmap
ComplexHeatmap::Heatmap(heatmapData,
                        top_annotation = heatmapColAnnot,
                        cluster_rows = TRUE, cluster_columns = FALSE)


####### EXTRA: Volcano plot
library(ggplot2)
library(ggrepel)

df <- as.data.frame(resTime)
df$gene <- rownames(df)
df$neglog10padj <- -log10(df$padj)

# thresholds
lfc_cut  <- 1
padj_cut <- 0.05

# classify genes
df$status <- "ns"
df$status[df$padj < padj_cut & df$log2FoldChange >=  lfc_cut] <- "Up"
df$status[df$padj < padj_cut & df$log2FoldChange <= -lfc_cut] <- "Down"

# find top 10 most significant within Up and Down
top_up <- df[df$status == "Up", ]
top_up <- head(top_up[order(top_up$padj), "gene"], 10)

top_down <- df[df$status == "Down", ]
top_down <- head(top_down[order(top_down$padj), "gene"], 10)

labels <- c(top_up, top_down)

# plot
ggplot(df, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(aes(color = status), alpha = 0.7, size = 1.6, na.rm = TRUE) +
  scale_color_manual(values = c("Down" = "blue", "ns" = "grey", "Up" = "red")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2) +
  geom_hline(yintercept = -log10(padj_cut), linetype = 2) +
  labs(x = "log2 fold change", y = expression(-log[10]~adjusted~p),
       color = NULL, title = "Volcano plot") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right") +
  geom_text_repel(
    data = subset(df, gene %in% labels),
    aes(label = gene),
    size = 3, box.padding = 0.3, point.padding = 0.2, max.overlaps = Inf
  )

# Output results
head(as.data.frame(resTime))
head(as.data.frame(rowRanges(se)))

temp <- cbind(as.data.frame(rowRanges(se)),
              as.data.frame(resTime))

write.csv(temp, file = "data/Day8vsDay0.csv")

