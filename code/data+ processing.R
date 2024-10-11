library(SummarizedExperiment)
library(TCGAbiolinks)
library(DESeq2)
library(maftools)
library(limma)
library(dplyr)
################### download methylation data
tcga_barcodes <- c(
  "TCGA-E1-A7YQ-01A-11D-A34K-05", "TCGA-S9-A7QW-01A-11D-A34D-05", "TCGA-HT-7874-01A-11D-2399-05", "TCGA-FG-A4MU-01B-11D-A28N-05", "TCGA-DU-8162-01A-21D-2254-05",
  "TCGA-DB-A4XG-01A-11D-A27L-05", "TCGA-FG-A60L-01A-12D-A31M-05", "TCGA-DU-A7TG-01A-21D-A34K-05", "TCGA-TM-A84O-01A-11D-A368-05", "TCGA-DB-A64Q-01A-11D-A29T-05",
  "TCGA-HT-8108-01A-11D-2399-05", "TCGA-DB-A64S-01A-11D-A29T-05", "TCGA-DB-A75P-01A-11D-A32C-05", "TCGA-DU-7007-01A-11D-2025-05", "TCGA-DH-A7UU-01A-12D-A34D-05",
  "TCGA-S9-A6U6-01A-12D-A33U-05", "TCGA-DU-5847-01A-11D-1706-05", "TCGA-HT-8113-01A-11D-2399-05", "TCGA-E1-A7Z2-01A-21D-A34K-05", "TCGA-S9-A7QY-01A-11D-A34D-05",
  "TCGA-FG-A87N-01A-11D-A368-05", "TCGA-HT-7676-01A-11D-2399-05", "TCGA-DU-A6S7-01A-21D-A32C-05", "TCGA-FG-A4MX-01A-11D-A26N-05", "TCGA-S9-A7IX-01A-12D-A34D-05",
  "TCGA-CS-6667-01A-12D-2025-05", "TCGA-FG-A6J3-01A-11D-A31M-05", "TCGA-S9-A6U8-01A-21D-A33U-05", "TCGA-DB-5275-01A-01D-1467-05", "TCGA-DU-8158-01A-11D-2254-05",
  "TCGA-CS-5397-01A-01D-1894-05", "TCGA-QH-A65Z-01A-11D-A29T-05", "TCGA-TQ-A7RN-01A-11D-A33U-05", "TCGA-DU-8163-01A-11D-2254-05", "TCGA-WY-A858-01A-11D-A368-05",
  "TCGA-DU-7012-01A-11D-2025-05", "TCGA-DH-A66F-01A-11D-A29T-05", "TCGA-FG-A4MT-01A-11D-A26N-05", "TCGA-HT-7473-01A-11D-2025-05", "TCGA-HT-A616-01A-11D-A29T-05",
  "TCGA-DU-A7TD-01A-12D-A34D-05", "TCGA-TQ-A7RJ-01A-11D-A33U-05", "TCGA-HW-7489-01A-11D-2025-05", "TCGA-DU-A5TS-01A-11D-A28N-05", "TCGA-S9-A6UB-01A-21D-A33U-05",
  "TCGA-QH-A6CY-01A-11D-A32C-05", "TCGA-DU-A7TA-01A-11D-A33U-05", "TCGA-HT-A615-01A-11D-A29T-05", "TCGA-HT-7690-01A-11D-2254-05", "TCGA-S9-A6WM-01A-12D-A33U-05",
  "TCGA-HT-7877-01A-11D-2399-05", "TCGA-E1-A7YE-01A-11D-A34D-05", "TCGA-DU-5870-01A-11D-1706-05", "TCGA-DU-8164-01A-11D-2254-05", "TCGA-HT-A5R5-01A-11D-A28N-05",
  "TCGA-TQ-A7RF-01A-11D-A33U-05", "TCGA-P5-A77W-01A-11D-A32C-05", "TCGA-RY-A847-01A-11D-A368-05", "TCGA-S9-A7R8-01A-11D-A34K-05", "TCGA-HT-A61C-01A-11D-A29T-05",
  "TCGA-S9-A7QZ-01A-12D-A34K-05", "TCGA-S9-A6U5-01A-12D-A33U-05", "TCGA-HT-7680-01A-11D-2254-05", "TCGA-DU-6395-01A-12D-1706-05", "TCGA-DH-A7US-01A-11D-A33U-05",
  "TCGA-FG-5964-01A-11D-1706-05", "TCGA-P5-A5EZ-01A-11D-A27L-05", "TCGA-CS-4944-01A-01D-1467-05", "TCGA-DU-A76R-01A-11D-A32C-05", "TCGA-DU-7301-01A-11D-2087-05",
  "TCGA-S9-A6TU-01A-12D-A32C-05", "TCGA-QH-A6X4-01A-51D-A32C-05", "TCGA-S9-A7R2-01A-21D-A34K-05", "TCGA-DB-5279-01A-01D-1467-05", "TCGA-HT-8012-01A-11D-2399-05",
  "TCGA-CS-4943-01A-01D-1467-05", "TCGA-HT-8558-01A-21D-2399-05", "TCGA-E1-A7YL-01A-11D-A34D-05", "TCGA-DU-A6S8-01A-12D-A32C-05", "TCGA-S9-A6WN-01A-12D-A33U-05",
  "TCGA-HT-8019-01A-21D-2399-05", "TCGA-S9-A6WO-01A-21D-A34D-05", "TCGA-S9-A6U0-01A-12D-32C-05", "TCGA-VM-A8C8-01A-11D-A368-05", "TCGA-E1-5307-01A-01D-1894-05",
  "TCGA-HT-7684-01A-11D-2254-05", "TCGA-QH-A86X-01A-11D-A368-05", "TCGA-DU-7014-01A-11D-2025-05", "TCGA-E1-A7YK-01A-11D-A34D-05", "TCGA-CS-5396-01A-02D-1467-05",
  "TCGA-DH-A7UT-01A-12D-A34D-05", "TCGA-FG-7636-01A-11D-2087-05", "TCGA-HW-8319-01A-11D-2399-05", "TCGA-FG-7643-01A-11D-2087-05"
)
methylation <- GDCquery(
  project = "TCGA-LGG",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  access = "open",
  data.type = "Methylation Beta Value",
  sample.type = c("Primary Tumor"),
  barcode = tcga_barcodes
)
methylationoutput = getResults(methylation)
GDCdownload(methylation)
gdcglio_meth= GDCprepare(methylation,summarizedExperiment = TRUE)
methdata= assay(gdcglio_meth)
# Save the data to a CSV file
write.csv(methdata, file = "final methylation data.csv", row.names = TRUE)

################### download mutation data
query_mut <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation")
mutation = getResults(query_mut)
GDCdownload(query_mut) 
mut_data = GDCprepare(query_mut)
mutationdata <- read.maf(maf = mut_data)

# Check and Extract IDH mutation status
idh_mutations <- subsetMaf(mutationdata, genes = c("IDH1", "IDH2"))
idh_status <- idh_mutations@data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")]
#save
write.csv(idh_status, file = "idh status final.csv", row.names = TRUE)

####@@@@@@@@@@@@@@@@@@@@@@@@@@@###################@@@@@@@@@@@@@@@@@@@@@@@@@@@
####loading data
methylationdata= read.csv("final methylation data.csv")
idh_status <- read.csv("idh status final.csv")

########### processing methylation data
# Remove columns with all NA values
methylationdata <- na.omit(methylationdata)

#Update column names to the desired format
colnames(methylationdata) <- sub("TCGA\\.", "TCGA-CS-", colnames(methylationdata))  # Prepend TCGA-CS-
colnames(methylationdata) <- gsub("\\.", "-", colnames(methylationdata))  # Replace dots with dashes

#Check the dimensions of both datasets
dim(methylationdata)
dim(idh_status)

