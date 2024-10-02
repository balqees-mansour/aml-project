# Install required packages
# BiocManager::install(c("curatedTCGAData", "MultiAssayExperiment","TCGAutils", "TCGAbiolinks" ))
# install.packages("data.table")
install.packages("gplots")
# load the required libraries
library(gplots)
library(ggplot2)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(pheatmap)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(DESeq2)
#library(readxl)
library(data.table)
library(dplyr)

# links for Normal Data GTEx Prepartions 
#01 read counts data
# https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz

# Sample annotation data 
# https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

dd <-  gtex_data[gtex_data$Description == "ARNT2", ]
# get a list of projects
gdcprojects <- getGDCprojects()
gdcprojects$disease_type
summarypro <- getProjectSummary('TCGA-LAML')

#tcga_ids <- read.csv("/home/admin/ids.csv",stringsAsFactors = FALSE)

######## 1 make an TCGA object #######
# build a query to retrieve gene expression data ------------
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



# Read GTEx normal blood samples 

# Set the path to your file
file_path <- "/home/balbio/Differential-Analysis/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"

# Read the GCT file
gtex_data <- fread(file_path, skip = 2)

# Read sample annotations
sample_info <- fread("/home/balbio/Differential-Analysis/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
gtex_data<- as.data.frame(gtex_data)


# Remove duplicates from the "Description" column
gtex_data_unique <- gtex_data[!duplicated(gtex_data$Description), ]

# Set the "Description" column as row names
rownames(gtex_data_unique) <- gtex_data_unique$Description

GG <- gtex_data_unique["ARNT2",]
gtex_data_unique <- gtex_data_unique[,-c(1,2)]



# Filter for whole blood samples
blood_samples <- sample_info[SMTSD == "Whole Blood" & SMAFRZE == "RNASEQ" & SMGEBTCHT == "TruSeq.v1" & SMNABTCHT == "RNA isolation_PAXgene Blood RNA (Manual)", SAMPID]


available_samples <- colnames(gtex_data_unique)  # Exclude 'Name' and 'Description'
matching_samples <- intersect(blood_samples, available_samples)

# Subset the expression data for whole blood samples
counts.data <- gtex_data_unique
colnames(counts.data)

normal.blood.data <- counts.data[,matching_samples]

normal.blood.data["ARNT2", ]
aml.matrix["ARNT2",]
## EXCLUDE    OMNI     RNASEQ     WES     WGS 
##   2700     450      17382      979     838 

intersected.genes.TN <- intersect(rownames(normal.blood.data), rownames(aml.matrix))

"ARNT2" %in% intersected.genes.TN
# Step 2: Subset the data frames based on the intersected genes
normal_blood_data <- normal.blood.data[intersected.genes.TN, , drop = FALSE]

subset_aml_matrix <- aml.matrix[intersected.genes.TN, , drop = FALSE]

tumor_selected_75 <- subset_aml_matrix[, sample(ncol(subset_aml_matrix), 75)]
# Randomly select 5 columns from a dataframe
random_selected <- normal_blood_data[, sample(ncol(normal_blood_data), 100)]
normal_selected_75 <- normal_blood_data[, sample(ncol(normal_blood_data), 75)]


aml.matrix["ARNT2",]
combined_data <- bind_cols(subset_aml_matrix, random_selected)
combine_150 <- bind_cols(tumor_selected_75, normal_selected_75)

combine_150["ARNT2",]


######## 5 convert Ensemble Ids to gene symbols ########

# Install and load org.Hs.eg.db package
library(org.Hs.eg.db)

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(combined_data)

# Remove version numbers for compatibility
ensembl_ids <- sub("\\..*", "", ensembl_ids)

# Map IDs
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Remove NA values and keep only mapped genes
mapped_genes <- !is.na(gene_symbols)
combined_data <- combined_data[mapped_genes, ]
gene_symbols <- gene_symbols[mapped_genes]

unique_gene_symbols <- make.unique(gene_symbols, sep = "_")
length(unique_gene_symbols)

# Assign gene symbols to prostate_matrix
rownames(combined_data) <- unique_gene_symbols

#=========================================================================
# my phenodata 
phdata_aml <- colnames(tumor_selected_75)
phdata_aml <- as.data.frame(phdata_aml)

# Assuming your data frame is called 'df'
phdata_aml$condition <- "aml"
names(phdata_aml)[1] <- "sample"

# add normal samples information 
ph_norm <- colnames(normal_selected_75)
ph_norm <- as.data.frame(ph_norm)

# Assuming your data frame is called 'df'
ph_norm$condition <- "normal"
names(ph_norm)[1] <- "sample"

samples_info <- bind_rows(phdata_aml,ph_norm)
rownames(samples_info) <- samples_info$sample

write.csv(combine_150, "countmatrixnew.csv")
write.csv(samples_info, "samplesheetnew.csv", row.names = FALSE)




#================================================================

count_data150 <- read.csv("countmatrixnew.csv")
rownames(count_data150) <- count_data150$gene_id
count_data150 <- count_data150[,-1]
colnames(count_data150) <- rownames(coldata150)

coldata150 <- read.csv("samplesheetnew.csv")
rownames(coldata150) <- coldata150$sample

coldata150$batch <- ifelse(coldata150$condition == "aml", "TCGA", "GTEx")
coldata150$batch[1] <- "GTEx"

  
  
  
# making the rownames and column names identical
all(rownames(coldata150) %in% colnames(count_data150))
all(rownames(coldata150) == colnames(count_data150))

coldata150$condition <- as.factor(coldata150$condition)
coldata150$batch <- as.factor(coldata150$batch)

count_data150 <- round(count_data150)
#dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~ batch + condition)
dds <- DESeqDataSetFromMatrix(countData = count_data150,
                              colData = coldata150,
                              design = ~ condition + batch) 

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition

# set the factor level
dds$condition <- relevel(dds$condition, ref = "normal")
levels(dds$condition)
# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)
summary(res)
combined_data["ARNT2",]
# Explore Results ----------------
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

# Filter for log2 fold change > 1 or < -1
res_filtered <- subset(res0.05, log2FoldChange > 1 | log2FoldChange < -1)
summary(res_filtered)

#==============================================================================
#-----------01# BoxPLot for the 20 upregulated genes in tumor vs normal ----------
#------------------------------------------------------------------------------
#-=============================================================================
install.packages("tidyverse")
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(tidyr)
# Assuming you've already run DESeq2 and have the results in 'res'

# Filter for upregulated genes (log2FoldChange > 1 and adjusted p-value < 0.05)
res_upregulated <- res_filtered[which(res_filtered$log2FoldChange > 1 & res_filtered$padj < 0.05),]

# Order by log2FoldChange
res_upregulated <- res_upregulated[order(res_upregulated$log2FoldChange, decreasing = TRUE),]

# Select top 18 upregulated genes
top_upregulated <- head(res_upregulated, 20)

# Extract normalized counts for these genes
normalized_counts <- counts(dds, normalized=TRUE)[rownames(top_upregulated),]

# Prepare data for plotting
plot_data <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(coldata150, by = c("sample" = "sample"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(res_filtered) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

# Create the plot
ggplot(plot_data, aes(x = condition, y = log2(count + 1), fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("normal" = "lightgray", "aml" = "darkgreen")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
            aes(x = condition[1], y = Inf, label = sprintf("q = %.2e", padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)


# Save the plot
ggsave("top_20_upregulated_genes.png", width = 15, height = 10, dpi = 300)

write.csv(res_upregulated, "top_20_upregulated_genes.csv")


#==============================================================================
#----- #02 BoxPLot for the top 20  downregulated genes in tumor vs normal -----
#------------------------------------------------------------------------------
#-=============================================================================
library(DESeq2)
library(tidyverse)
library(ggplot2)

# Assuming you've already run DESeq2 and have the results in 'res'

# Filter for upregulated genes (log2FoldChange > 1 and adjusted p-value < 0.05)
res_downregulated <- res_filtered[which(res_filtered$log2FoldChange < -1 & res_filtered$padj < 0.05),]

# Order by log2FoldChange
res_downregulated <- res_downregulated[order(res_downregulated$log2FoldChange, decreasing = FALSE),]

# Select top 18 upregulated genes
top_downregulated <- head(res_downregulated, 20)
arnt2 <- res_downregulated["ARNT2",]
# Extract normalized counts for these genes
normalized_counts <- counts(dds, normalized=TRUE)[rownames(arnt2),]

# Prepare data for plotting
plot_data <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(coldata150, by = c("sample" = "sample"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(res_filtered) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

# Create the plot
ggplot(plot_data, aes(x = condition, y = log2(count + 1), fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("normal" = "lightpink", "aml" = "purple")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
            aes(x = condition[1], y = Inf, label = sprintf("q = %.2e", padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)


# Save the plot
ggsave("top_20_downregulated_genes.png", width = 15, height = 10, dpi = 300)
write.csv(res_downregulated,"Downregulated_20genes.csv")


#==============================================================================
#----- #02 BoxPLot for ARNT2 gene in downregulated deseq2 result -----
#------------------------------------------------------------------------------
#-=============================================================================
library(DESeq2)
library(tidyverse)
library(ggplot2)

# Assuming you've already run DESeq2 and have the results in 'res'

# Filter for downregulated genes (log2FoldChange < -3 and adjusted p-value < 0.05)
res_downregulated <- res_filtered[which(res_filtered$log2FoldChange < -1 & res_filtered$padj < 0.05),]


normalized_counts <- counts(dds, normalized=TRUE)["ARNT2",]
  
  # Prepare data for plotting
plot_data <- normalized_counts %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
    left_join(coldata150, by = c("sample" = "sample"))
  
# Add log2FoldChange and padj values
plot_data <- plot_data %>%
    left_join(as.data.frame(res_filtered) %>% 
                rownames_to_column("gene") %>%
                dplyr::select(gene, log2FoldChange, padj),
                by = "gene")

ggplot(plot_data, aes(x = condition, y = log2(count + 1), fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("normal" = "lightpink", "aml" = "purple")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
            aes(x = condition[1], y = Inf, label = sprintf("q = %.2e", padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)

# Assuming you've already run DESeq2 and have the results in 'res'

# Filter for downregulated genes (log2FoldChange < -1 and adjusted p-value < 0.05)
res_downregulated <- res_filtered[which(res_filtered$log2FoldChange < -1 & res_filtered$padj < 0.05),]

# Find ARNT2 in the downregulated results
arnt2_result <- res_downregulated["ARNT2",]

# Extract normalized counts for ARNT2
normalized_counts <- counts(dds, normalized=TRUE)["ARNT2",]

# Prepare data for plotting
plot_data <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(coldata150, by = c("sample" = "sample"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(res_filtered) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

# Create the plot
ggplot(plot_data, aes(x = condition, y = log2(count + 1), fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  scale_fill_manual(values = c("normal" = "lightpink", "aml" = "purple")) +
  labs(y = "log2CPM", x = NULL, title = "ARNT2 Expression") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% slice(1),
            aes(x = condition[1], y = Inf, 
                label = sprintf("log2FC = %.2f\nq = %.2e", log2FoldChange, padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)

# Save the plot
ggsave("ARNT2_expression.png", width = 6, height = 4, dpi = 300)




#=============================================================================
#=============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("countmatrixnew.csv", row.names = 1)

# Assuming the first 75 columns are normal samples and the rest are AML samples
# Adjust these numbers if the split is different in your data
normal_samples <- 76:ncol(data)

aml_samples <- 1:75

# Extract ARNT2 data
arnt2_data <- data["ARNT2", ]

# Prepare the data for plotting
plot_data <- data.frame(
  condition = c(rep("normal", length(normal_samples)), rep("aml", length(aml_samples))),
  count = c(arnt2_data[normal_samples], arnt2_data[aml_samples])
)

# Create the plot
ggplot(plot_data, aes(x = condition, y = log2(count + 1), fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  scale_fill_manual(values = c("normal" = "lightpink", "aml" = "purple")) +
  labs(title = "Boxplot of ARNT2 Counts in Normal and AML Conditions",
       y = "log2(count + 1)", 
       x = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1),
                     labels = scales::scientific_format(digits = 1)) +
  coord_cartesian(ylim = c(0, 10))

# To add log2FoldChange and padj values, you would need to calculate these
# from your actual data using a differential expression analysis tool like DESeq2.
# Then you can add this annotation:
# annotate("text", x = 1.5, y = 10, 
#          label = sprintf("log2FC = %.2f\nq = %.2e", log2FoldChange, padj),
#          hjust = 0.5, vjust = 1)



library(ggplot2)
library(dplyr)

# Read the data
data <- read.csv("countmatrixnew.csv", row.names = 1)

# Extract ARNT2 data
arnt2_data <- data["ARNT2", ]

# Prepare the data for plotting
plot_data <- data.frame(
  condition = factor(c(rep("aml", 75), rep("normal", ncol(data) - 75))),
  count = as.numeric(arnt2_data)
)

# Remove any NA values
plot_data <- plot_data[complete.cases(plot_data), ]

# Create the plot
ggplot(plot_data, aes(x = condition, y = log2(count + 1), fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.4, alpha = 0.5) +
  scale_fill_manual(values = c("normal" = "lightpink", "aml" = "purple")) +
  labs(title = "Boxplot of ARNT2 Counts in Normal and AML Conditions",
       y = "log2(count + 1)", 
       x = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1),
                     labels = scales::scientific_format(digits = 1)) +
  coord_cartesian(ylim = c(0, 10))


#






