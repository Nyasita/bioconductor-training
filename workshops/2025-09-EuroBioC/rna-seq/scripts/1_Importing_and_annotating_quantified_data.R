#################################################
# Data preparation
# Importing and annotating quantified data into R
#################################################

# Packages that are needed for the workshop

install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("tidyverse", "SummarizedExperiment",
                       "ExploreModelMatrix", "AnnotationDbi", "org.Hs.eg.db", 
                       "org.Mm.eg.db", "csoneson/ConfoundingExplorer",
                       "DESeq2", "vsn", "ComplexHeatmap", "hgu95av2.db",
                       "RColorBrewer", "hexbin", "cowplot", "iSEE",
                       "clusterProfiler", "enrichplot", "kableExtra",
                       "msigdbr", "gplots", "ggplot2", "simplifyEnrichment",
                       "apeglm", "microbenchmark", "Biostrings",
                       "SingleCellExperiment"))

# Avoid the warning messages
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(hgu95av2.db)
  library(SummarizedExperiment)
})

#Load the Data
#Counts
counts <- read.csv("data/GSE96870_counts_cerebellum.csv", 
                   row.names = 1)
dim(counts)

#Sample annotations
coldata <- read.csv("data/GSE96870_coldata_cerebellum.csv",
                    row.names = 1)
dim(coldata)

#Gene annotations
rowranges <- read.delim("data/GSE96870_rowranges.tsv", 
                        sep = "\t", 
                        colClasses = c(ENTREZID = "character"),
                        header = TRUE, 
                        quote = "", 
                        row.names = 5)
dim(rowranges)

table(rowranges$gbkey)

#Assemble SummarizedExperiment
#Check if samples and genes are in the same order
all.equal(colnames(counts), rownames(coldata)) # samples

all.equal(rownames(counts), rownames(rowranges)) # genes

# If the first is not TRUE, you can match up the samples/columns in
# counts with the samples/rows in coldata like this (which is fine
# to run even if the first was TRUE):

tempindex <- match(colnames(counts), rownames(coldata))
coldata <- coldata[tempindex, ]

# Check again:
all.equal(colnames(counts), rownames(coldata)) 

#Create the SummarizedExperiment oblect
# One final check:
stopifnot(rownames(rowranges) == rownames(counts), # features
          rownames(coldata) == colnames(counts)) # samples

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  rowRanges = as(rowranges, "GRanges"),
  colData = coldata
)

# Access the counts
head(assay(se))
dim(assay(se))

# The above works now because we only have one assay, "counts"
# But if there were more than one assay, we would have to specify
# which one like so:
head(assay(se, "counts"))

# Access the sample annotations
colData(se)

dim(colData(se))

# Access the gene annotations
head(rowData(se))

dim(rowData(se))

# Make better sample IDs that show sex, time and mouse ID:

se$Label <- paste(se$sex, se$time, se$mouse, sep = "_")
se$Label

colnames(se) <- se$Label

# Our samples are not in order based on sex and time
se$Group <- paste(se$sex, se$time, sep = "_")
se$Group

# change this to factor data with the levels in order 
# that we want, then rearrange the se object:

se$Group <- factor(se$Group, levels = c("Female_Day0","Male_Day0", 
                                        "Female_Day4","Male_Day4",
                                        "Female_Day8","Male_Day8"))
se <- se[, order(se$Group)]
colData(se)

# Finally, also factor the Label column to keep in order in plots:

se$Label <- factor(se$Label, levels = se$Label)

# Questions:
# 1
table(se$infection)

# 2
se_infected <- se[, se$infection == "InfluenzaA"]
se_noninfected <- se[, se$infection == "NonInfected"]

means_infected <- rowMeans(assay(se_infected)[1:500, ])
means_noninfected <- rowMeans(assay(se_noninfected)[1:500, ])

summary(means_infected)
summary(means_noninfected)

# 3
ncol(se[, se$sex == "Female" & se$infection == "InfluenzaA" & se$time == "Day8"])

#Save SummarizedExperiment

saveRDS(se, "data/GSE96870_se.rds")
rm(se) # remove the object!

#Load the object
se <- readRDS("data/GSE96870_se.rds")

#Subset the se into mRNA
se_mRNA <- se[rowData(se)$gbkey == "mRNA" , ]
dim(se_mRNA)

sessionInfo()


# EXTRA: Gene Annotations

# mapIds(annopkg, keys, column, keytype, ..., multiVals)
mapIds(org.Mm.eg.db, keys = "497097", column = "SYMBOL", keytype = "ENTREZID")

# 1:n relationship
keys <- head(keys(hgu95av2.db, "ENTREZID"))
last <- function(x){x[[length(x)]]}

# Getting the first
mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID")

# Getting the last
mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID", multiVals = last)

# Getting all the matches
mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID", multiVals = "list")
