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
#BiocManager::install("EnhancedVolcano")

# Step 2: Load the package
library(EnhancedVolcano)

# Read normal and AML data 
normal.full <- fread("/home/balqees/ola-project/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
normal.full$Name

library(readr)
normal_blood_GTEx_data <- normal.full 


# Read sample annotations
sample_info <- fread("/home/balqees/ola-project/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

# Filter for whole blood samples
blood_samples <- sample_info[SMTSD == "Whole Blood" & SMAFRZE == "RNASEQ" & SMGEBTCHT == "TruSeq.v1" & SMNABTCHT == "RNA isolation_PAXgene Blood RNA (Manual)", SAMPID]

rownames(sample_info) <- sample_info$SAMPID

df.sample.info <- as.data.frame(sample_info)
rownames(df.sample.info) <-  df.sample.info$SAMPID

nblood.data <- df.sample.info[blood_samples,]

LCSET10313.normal.info <- nblood.data[which(nblood.data$SMGEBTCH == "LCSET-10313"),]

normal_blood_GTEx_data <- as.data.frame(normal_blood_GTEx_data)


rownames(normal_blood_GTEx_data) <-  normal_blood_GTEx_data[,1]
normal_blood_GTEx_data <- normal_blood_GTEx_data[,-1]

row.counts.LCSET10313.normal <- normal_blood_GTEx_data[,LCSET10313.normal.info$SAMPID]
write.csv(row.counts.LCSET10313.normal, "row.counts.LCSET10313.normal28.csv")
dim(row.counts.LCSET10313.normal)

row.counts.LCSET10313.normal <- read.csv("row.counts.LCSET10313.normal28.csv")
rownames(row.counts.LCSET10313.normal)     <- row.counts.LCSET10313.normal$X
row.counts.LCSET10313.normal <- row.counts.LCSET10313.normal[,-1]
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

rownames(aml.subset) %in% rownames(row.counts.LCSET10313.normal) %>% sum()


# Using grep to find the matching pattern # "ENSG00000029153.14" 
pattern <- "ENSG00000029153"
rownames(aml.subset)[grep(pattern, rownames(aml.subset))]
rownames(row.counts.LCSET10313.normal)[grep(pattern, rownames(row.counts.LCSET10313.normal))]


# Combine TCGA and GTEx data # 35117 common genes btw normal and aml
common_genes <- intersect(rownames(aml.subset), rownames(row.counts.LCSET10313.normal))
combined_data <- cbind(aml.subset[common_genes,], row.counts.LCSET10313.normal[common_genes,])

combined_data["ENSG00000029153.14",]
# Create sample information dataframe
sample_info <- data.frame(
  condition = c(rep("AML", ncol(aml.subset)), rep("Normal", ncol(row.counts.LCSET10313.normal))),
  batch = c(rep("TCGA", ncol(aml.subset)), rep("GTEx", ncol(row.counts.LCSET10313.normal))),
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


res <- results(dds, 
               alpha = 0.05,        # Sets FDR threshold
               lfcThreshold = 2,    # Sets log2 fold change threshold
               altHypothesis = "greaterAbs")

summary(res)



results(dds)
# Filter and sort results
significant_genes <- result[which(result$padj < 0.05 & abs(result$log2FoldChange) > 2),]
significant_genes <- significant_genes[order(significant_genes$padj),]

summary(significant_genes)
# Visualize results
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'AML vs Normal')

# Optional: Save results to a file
write.csv(significant_genes, file = "AML_vs_Normal_DEGs.csv")

res[which(rownames(res) == "ENSG00000029153.14"),]




