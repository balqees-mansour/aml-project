# Load required libraries
BiocManager::install("recount3")
library(DESeq2)
library(TCGAbiolinks)
library(recount3)
library(dplyr)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(EnhancedVolcano)

# 1- Download and prepare TCGA AML data
query_TCGA <- GDCquery(project = 'TCGA-LAML',data.category = c('Transcriptome Profiling') ,experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts', access = 'open')



# download data - GDCdownload
####### 2 download the object ########
GDCdownload(query_TCGA, directory = "/home/balbio/Differential-Analysis/")


####### 3 prepare data #########
aml.tcga.data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE,directory = "/home/balbio/Differential-Analysis/")

######## 4 get the assay ########
aml.raw <- assay(aml.tcga.data, 'stranded_first')
aml.raw <- as.data.frame(aml.raw)

######## 5 convert Ensemble Ids to gene symbols ########
# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(aml.raw)

# Remove version numbers for compatibility
ensembl_ids <- sub("\\..*", "", ensembl_ids)

# Map IDs
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Remove NA values and keep only mapped genes
mapped_genes <- !is.na(gene_symbols)
aml.raw <- aml.raw[mapped_genes, ]
gene_symbols <- gene_symbols[mapped_genes]

unique_gene_symbols <- make.unique(gene_symbols, sep = "_")
length(unique_gene_symbols)

# Assign gene symbols to prostate_matrix
rownames(aml.raw) <- unique_gene_symbols
# Write the prepared AML data 151 samples 
write.csv(aml.raw, "prepared_AML_Raw151.csv")

normal_gtex <- read.csv("normal_blood_GTEx_data.csv")
rownames(normal_gtex) <- normal_gtex$X
normal_gtex <- normal_gtex[,-1]



# Subset GTEx data to match TCGA sample size
set.seed(123)  # for reproducibility
gtex_subset <- normal_gtex[, sample(ncol(normal_gtex), 150)]

# Combine TCGA and GTEx data
common_genes <- intersect(rownames(aml.raw), rownames(gtex_subset))
aml.raw <- aml.raw[common_genes,]
gtex_subset <- gtex_subset[common_genes,]

combined_data <- cbind(aml.raw[common_genes,], gtex_subset[common_genes,])

# Create sample information dataframe
sample_info <- data.frame(
  condition = c(rep("AML", ncol(aml.raw)), rep("Normal", ncol(gtex_subset))),
  batch = c(rep("TCGA", ncol(aml.raw)), rep("GTEx", ncol(gtex_subset))),
  row.names = colnames(combined_data)
)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(combined_data),
                              colData = sample_info,
                              design = ~ condition )

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
results <- results(dds, contrast = c("condition", "AML", "Normal"))

# Filter and sort results
significant_genes <- results[which(results$padj < 0.05 & abs(results$log2FoldChange) > 1),]
significant_genes <- significant_genes[order(significant_genes$padj),]

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
                
                
#====================================================================================
#====================================================================================
sample_info <- fread("/home/balbio/Differential-Analysis/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
tcga.clinical <- as.data.frame(aml.tcga.data@colData)

write.csv(tcga.clinical, "tcga.clinical.csv")

subset_aml <- aml.tcga.data[,which(aml.tcga.data@colData@listData[["prior_malignancy"]] == "no" & aml.tcga.data@colData@listData[["prior_treatment"]] == "No")]

dim (subset_aml)

dim (aml.tcga.data)


subset_aml$treatments[[10]][5:15]

subset_aml$treatments[[12]][5] == "Radiation Therapy, NOS"

subset_aml$treatments[15:16][1:5]

    