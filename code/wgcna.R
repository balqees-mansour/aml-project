# script to perform WGCNA
# setwd("~/Desktop/demo/WGCNA")

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()          # allow multi-threading (optional)

# 1. Fetch Data ------------------------------------------------
data=read.csv("combined_data.csv")
rownames(data)<-data[,1]
data<-data[,-1]
phenodata=read.csv("combined_anno.csv")
rownames(phenodata)<-phenodata$X
phenodata<-phenodata[,-1]
# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


### NOTE: If there are batch effects observed, correct for them before moving ahead


# 3. Normalization ----------------------------------------------------------------------


# making the rownames and column names identical
# Step 1: Replace "." with "-" in the column names of normal_blood_GTEx_data
colnames(data) <- gsub("\\.", "-", colnames(data))
all(rownames(phenodata) %in% colnames(data))
all(rownames(phenodata) == colnames(data))


# create dds
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = phenodata,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more than 75% of samples (129*0.75=97)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 97,]
nrow(dds75) # 13284 genes


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

#checking our gene 
"ENSG00000029153.14" %in% norm.counts

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 10
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)




# grey module = all genes that doesn't fall into other modules were assigned to the grey module


# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

# create traits file - binarize categorical variables
traits <- phenodata %>% 
  mutate(disease_state_bin = ifelse(grepl('AML', condition), 1, 0)) %>% 
  select(3)


# binarize categorical variables

phenodata$batch <- factor(phenodata$batch, levels = c("TCGA", "GTEX"))

batch.out <- binarizeCategoricalColumns(phenodata$batch,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)


traits <- cbind(traits, batch.out)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data),
             y = names(heatmap.data)[1:3],
             col = c("blue1", "skyblue", "white", "pink", "red"))



module.gene.mapping <- as.data.frame(bwnet$colors)
view(module.gene.mapping)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()



# 6B. Intramodular analysis: Identifying driver genes ---------------



# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)



# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$data.Severe.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


x=gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(rownames(gene.signf.corr.pvals)) %>% 
  head(25)
rownames(x)
write.csv(x,"hubgenes.csv")

#GENE Sympols
ensembl_ids <- rownames(x)

# Remove version numbers for compatibility
ensembl_ids <- sub("\\..*", "", ensembl_ids)

# Map IDs
library(AnnotationDbi)
library(org.Hs.eg.db)

gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Remove NA values and keep only mapped genes
mapped_genes <- !is.na(gene_symbols)
gene_symbols <- gene_symbols[mapped_genes]
genes=as.data.frame(gene_symbols)
genes$gene_symbols
# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.