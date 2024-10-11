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
library("BiocParallel")
register(MulticoreParam(32))
# Read the normal data 
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
gtex_data_unique <- gtex_data_unique[,-c(1,2)]

# Read sample annotations
sample_info <- fread("/home/balbio/Differential-Analysis/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

# Filter for whole blood samples
blood_samples <- sample_info[SMTSD == "Whole Blood" & SMAFRZE == "RNASEQ" & SMGEBTCHT == "TruSeq.v1" & SMNABTCHT == "RNA isolation_PAXgene Blood RNA (Manual)", SAMPID]


available_samples <- colnames(gtex_data_unique)  # Exclude 'Name' and 'Description'
matching_samples <- intersect(blood_samples, available_samples)


# Subset the expression data for whole blood samples
normal_counts <- gtex_data_unique[,matching_samples]

ahmad <-  read_xlsx("/home/balbio/Differential-Analysis/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx")

assad <- fread("/home/balbio/Differential-Analysis/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

aly <-  read_xlsx("/home/balbio/Differential-Analysis/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx")

write.csv(aly$VARDESC , "vardesc.csv")































