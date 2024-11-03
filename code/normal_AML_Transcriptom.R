# Library Loadings
library(gplots)
library(ggplot2)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(pheatmap)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(DESeq2)
library(data.table)
library(dplyr)
library(sva)
library(readxl)
library(org.Hs.eg.db)
library("BiocParallel")
register(MulticoreParam(12))


# Read normal and AML data 
normal.counts28 <- read.csv("row.counts.LCSET10313.normal28.csv")
rownames(normal.counts28) <- normal.counts28$X
normal.counts28 <- normal.counts28[,-1]
normal.meta.data28 <- read.csv("LCSET10313.normal28.info.csv")


# AML data 

# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-LAML',data.category = c('Transcriptome Profiling') ,experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts', access = 'open')

# download data - GDCdownload
####### 2 download the object ########
GDCdownload(query_TCGA, directory = "/home/balqees/ola-project")


####### 3 prepare data #########
aml.tcga.data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE,directory = "/home/balqees/ola-project")




######## 4 get the assay ########
aml.raw <- assay(aml.tcga.data, 'stranded_first')
aml.raw <- as.data.frame(aml.raw)

barcodes_101 <- aml.tcga.data@colData@listData$barcode[
  aml.tcga.data@colData@listData[["prior_malignancy"]] == "no" & 
    aml.tcga.data@colData@listData[["prior_treatment"]] == "No"
]
aml.subset <- aml.raw[,barcodes_101]
dim(aml.subset)
######## 5 convert Ensemble Ids to gene symbols ########
# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(aml.subset)

# Remove version numbers for compatibility
ensembl_ids <- sub("\\..*", "", ensembl_ids)

# Map IDs
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Remove NA values and keep only mapped genes
mapped_genes <- !is.na(gene_symbols)
aml.subset <- aml.subset[mapped_genes, ]
gene_symbols <- gene_symbols[mapped_genes]

unique_gene_symbols <- make.unique(gene_symbols, sep = "_")
length(unique_gene_symbols)


# Assign gene symbols to prostate_matrix
rownames(aml.subset) <- unique_gene_symbols
# Write the prepared AML data 151 samples 
write.csv(aml.subset, "aml.subset101.csv")


# Method 1: Using rownames
rownames(aml.subset)[rownames(aml.subset) == "BMAL2"] <- "ARNTL2"
aml.subset["ARNTL2",]
normal.counts28["ARNTL2",]
# Combine TCGA and GTEx data
common_genes <- intersect(rownames(aml.subset), rownames(normal.counts28))
combined_data <- cbind(aml.subset[common_genes,], normal.counts28[common_genes,])

combined_data["ARNTL2",]
# Create sample information dataframe
sample_info <- data.frame(
  condition = c(rep("AML", ncol(aml.subset)), rep("Normal", ncol(normal.counts28))),
  batch = c(rep("TCGA", ncol(aml.subset)), rep("GTEx", ncol(normal.counts28))),
  row.names = colnames(combined_data)
)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(combined_data),
                              colData = sample_info,
                              design = ~ condition )


# set the factor level
dds$condition <- relevel(dds$condition, ref = "Normal")
levels(dds$condition)
# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)

result <- results(dds)
summary(result)
vsd <- vst(dds)


result(dds)
# Filter and sort results
significant_genes <- result[which(result$padj < 0.05 & abs(result$log2FoldChange) > 2),]
significant_genes <- significant_genes[order(significant_genes$padj),]

summary(significant_genes)
# Visualize results
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'AML vs Normal')

# Optional: Save results to a file
write.csv(significant_genes, file = "AML_vs_Normal_DEGs.csv")






