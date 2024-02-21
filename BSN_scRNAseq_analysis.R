library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(hdf5r)

# Loading the dataset 2 
#assuming /home/rstudio/BSN_NewData/filtered_h5s contains data from set 2, modify as approppriate

data.dir="/home/rstudio/BSN_NewData/filtered_h5s"
setwd(data.dir)

#set 2 data
wt1.data <- Read10X_h5("SITTH2.filtered_feature_bc_matrix.h5")
wt2.data <- Read10X_h5("SITTD3.filtered_feature_bc_matrix.h5")
wt3.data <- Read10X_h5("SITTH3.filtered_feature_bc_matrix.h5")
het1.data <- Read10X_h5("SITTE2.filtered_feature_bc_matrix.h5")
het2.data <- Read10X_h5("SITTA3.filtered_feature_bc_matrix.h5")
het3.data <- Read10X_h5("SITTE3.filtered_feature_bc_matrix.h5")
het4.data <- Read10X_h5("SITTF2.filtered_feature_bc_matrix.h5")
het5.data <- Read10X_h5("SITTB3.filtered_feature_bc_matrix.h5")
het6.data <- Read10X_h5("SITTF3.filtered_feature_bc_matrix.h5")


# Loading the dataset 1
data.dir="/home/rstudio/BSN_NewData/Batch_1_20221110"
#assuming /home/rstudio/BSN_NewData/Batch_1_20221110 contains data from set 1, modify as approppriate

setwd(data.dir)

wt4.data <- Read10X(data.dir = "SITTF10_count/outs/filtered_feature_bc_matrix/")
wt5.data <- Read10X(data.dir = "SITTA11_count/outs/filtered_feature_bc_matrix/")
wt6.data <- Read10X(data.dir = "SITTB11_count/outs/filtered_feature_bc_matrix/")
het7.data <- Read10X(data.dir = "SITTH1_count/outs/filtered_feature_bc_matrix/")
het8.data <- Read10X(data.dir = "SITTA1_count/outs/filtered_feature_bc_matrix/")
het9.data <- Read10X(data.dir = "SITTB1_count/outs/filtered_feature_bc_matrix/")

# Initialise object

wt1 <- CreateSeuratObject(counts = wt1.data)
wt2 <- CreateSeuratObject(counts = wt2.data)
wt3 <- CreateSeuratObject(counts = wt3.data)
wt4 <- CreateSeuratObject(counts = wt4.data)
wt5 <- CreateSeuratObject(counts = wt5.data)
wt6 <- CreateSeuratObject(counts = wt6.data)
het1 <- CreateSeuratObject(counts = het1.data)
het2 <- CreateSeuratObject(counts = het2.data)
het3 <- CreateSeuratObject(counts = het3.data)
het4 <- CreateSeuratObject(counts = het4.data)
het5 <- CreateSeuratObject(counts = het5.data)
het6 <- CreateSeuratObject(counts = het6.data)
het7 <- CreateSeuratObject(counts = het7.data)
het8 <- CreateSeuratObject(counts = het8.data)
het9 <- CreateSeuratObject(counts = het9.data)


# Calculate percentage of mitochondiral data

wt1[["percent.mt"]] <- PercentageFeatureSet(wt1, pattern = "^MT-")
wt2[["percent.mt"]] <- PercentageFeatureSet(wt2, pattern = "^MT-")
wt3[["percent.mt"]] <- PercentageFeatureSet(wt3, pattern = "^MT-")
wt4[["percent.mt"]] <- PercentageFeatureSet(wt4, pattern = "^MT-")
wt5[["percent.mt"]] <- PercentageFeatureSet(wt5, pattern = "^MT-")
wt6[["percent.mt"]] <- PercentageFeatureSet(wt6, pattern = "^MT-")
het1[["percent.mt"]] <- PercentageFeatureSet(het1, pattern = "^MT-")
het2[["percent.mt"]] <- PercentageFeatureSet(het2, pattern = "^MT-")
het3[["percent.mt"]] <- PercentageFeatureSet(het3, pattern = "^MT-")
het4[["percent.mt"]] <- PercentageFeatureSet(het4, pattern = "^MT-")
het5[["percent.mt"]] <- PercentageFeatureSet(het5, pattern = "^MT-")
het6[["percent.mt"]] <- PercentageFeatureSet(het6, pattern = "^MT-")
het7[["percent.mt"]] <- PercentageFeatureSet(het7, pattern = "^MT-")
het8[["percent.mt"]] <- PercentageFeatureSet(het8, pattern = "^MT-")
het9[["percent.mt"]] <- PercentageFeatureSet(het9, pattern = "^MT-")


# Adding columns for Clone (Wt vs HetB vs Hom) and Replicates

wt1@meta.data['Clone'] <- 'WT'

wt1@meta.data['Replicate'] <- 'wt1'

wt2@meta.data['Clone'] <- 'WT'

wt2@meta.data['Replicate'] <- 'wt2'

wt3@meta.data['Clone'] <- 'WT'

wt3@meta.data['Replicate'] <- 'wt3'

wt4@meta.data['Clone'] <- 'WT'

wt4@meta.data['Replicate'] <- 'wt4'

wt5@meta.data['Clone'] <- 'WT'

wt5@meta.data['Replicate'] <- 'wt5'

wt6@meta.data['Clone'] <- 'WT'

wt6@meta.data['Replicate'] <- 'wt6'

het1@meta.data['Clone'] <- 'HET'

het1@meta.data['Replicate'] <- 'het1'

het2@meta.data['Clone'] <- 'HET'

het2@meta.data['Replicate'] <- 'het2'

het3@meta.data['Clone'] <- 'HET'

het3@meta.data['Replicate'] <- 'het3'

het4@meta.data['Clone'] <- 'HET'

het4@meta.data['Replicate'] <- 'het4'

het5@meta.data['Clone'] <- 'HET'

het5@meta.data['Replicate'] <- 'het5'

het6@meta.data['Clone'] <- 'HET'

het6@meta.data['Replicate'] <- 'het6'

het7@meta.data['Clone'] <- 'HET'

het7@meta.data['Replicate'] <- 'het7'

het8@meta.data['Clone'] <- 'HET'

het8@meta.data['Replicate'] <- 'het8'

het9@meta.data['Clone'] <- 'HET'

het9@meta.data['Replicate'] <- 'het9'


# Adding columns for Batch

wt1@meta.data['Batch'] <- 'Batch_2'

wt2@meta.data['Batch'] <- 'Batch_2'

wt3@meta.data['Batch'] <- 'Batch_2'

het1@meta.data['Batch'] <- 'Batch_2'

het2@meta.data['Batch'] <- 'Batch_2'

het3@meta.data['Batch'] <- 'Batch_2'

het4@meta.data['Batch'] <- 'Batch_2'

het5@meta.data['Batch'] <- 'Batch_2'

het6@meta.data['Batch'] <- 'Batch_2'


wt4@meta.data['Batch'] <- 'Batch_1'

wt5@meta.data['Batch'] <- 'Batch_1'

wt6@meta.data['Batch'] <- 'Batch_1'

het7@meta.data['Batch'] <- 'Batch_1'

het8@meta.data['Batch'] <- 'Batch_1'

het9@meta.data['Batch'] <- 'Batch_1'


# merge datasets

ls()

merged_samples <- merge(wt1, y = c(wt2, wt3, wt4, wt5, wt6,
                                   het1, het2, het3, het4, het5, het6, het7, het8, het9),
                        add.cell.ids = c("wt1", "wt2", "wt3", "wt4", "wt5", "wt6", "het1", "het2", "het3", "het4", "het5", "het6", "het7", "het8", "het9"))

# filtering - merge_merged_seurat_3_samples_filtered_filtered according to merged_seurat_3_samples_filtered vignette filter, merge_merged_seurat_3_samples_filtered_filtered2 according to Bioinformagician

merged_samples <- subset(merged_samples, subset = nCount_RNA > 800 &
                           nFeature_RNA > 500 & nFeature_RNA < 10000 &
                           percent.mt < 5)

# Violin plot of the data

VlnPlot(merged_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Preparing the dataset for SCTransform

merged_samples_split <- SplitObject(merged_samples, split.by = "Replicate")

merged_samples_split <- lapply(X = merged_samples_split, assay = 'RNA', FUN = SCTransform)

# Finding integration features and anchors 

features <- SelectIntegrationFeatures(object.list = merged_samples_split, nfeatures = 3000)

merged_samples_split <- PrepSCTIntegration(object.list = merged_samples_split, anchor.features = features)

object.anchors <- FindIntegrationAnchors(object.list = merged_samples_split, normalization.method = "SCT",
                                         anchor.features = features)
# Integrate object

object.integrated <- IntegrateData(anchorset = object.anchors, normalization.method = "SCT")

# Run PCA

object.integrated <- RunPCA(object.integrated, verbose = FALSE)
object.integrated <- RunUMAP(object.integrated, reduction = "pca", dims = 1:30)

# View UMAP plots

p1 <- DimPlot(object.integrated, reduction = "umap", group.by = "Clone")
p2 <- DimPlot(object.integrated, reduction = "umap", group.by = "Replicate")
p3 <- DimPlot(object.integrated, reduction = "umap", group.by = "Batch")


p1 + p2 + p3

# Save the initial object

saveRDS(object.integrated, file = "integrated_object.rds")

# Find clusters

object.integrated <- FindNeighbors(object.integrated, reduction = "pca", dims = 1:30)
object.integrated <- FindClusters(object.integrated, resolution = 0.4)

# Plot with clusters

p4 <- DimPlot(object.integrated, reduction = "umap", label = TRUE, repel = TRUE)

p4

# Back to RNA assay

DefaultAssay(object.integrated) <- "RNA"

# Looking at markers

FeaturePlot(object.integrated.normalized, features = c("RBFOX3", "SLC1A2", "SLC1A3", "SLC17A6", "SLC17A7", "GAD1", "SLC32A1"), ncol = 3)

# Normalizing the data

object.integrated.normalized <-NormalizeData(object = object.integrated, normalization.method = "LogNormalize", 
                                             scale.factor = 10000)

# Finding markers for all the clusters

object.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
object.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Save the table of all markers

write.table(object.markers,"object.integrated.normalized.markers.csv",quote=F,col.names=NA,sep="/")

# Het vs WT markers quick check

levels(data)

colnames(object.integrated.normalized@meta.data)

Idents(object = data) <- "seurat_clusters"

wt_het.dif.markers <- FindMarkers(object.integrated.normalized, ident.1 = "WT", ident.2 = "HET")

head(wt_het.dif.markers, n = 10)

# Saving a normalized RDS object

saveRDS(object.integrated.normalized, file = "normalized_integrated_object.rds")


# Preparing the data for the deSEQ analysis

counts <- AggregateExpression(data, group.by = c("seurat_clusters", "Sample"),
                              assays = 'RNA',
                              slot = "counts",
                              return.seurat = FALSE)


counts_clusters_merged <- AggregateExpression(data, group.by = c( "Sample"),
                                              assays = 'RNA',
                                              slot = "counts",
                                              return.seurat = FALSE)



# Get the count data from the matrix

counts <- counts$RNA

# Transpose the matrix

counts.tr <- t(counts)

# Convert to a data frame

counts.tr <- as.data.frame(counts.tr)

# Split the matrix according to the cluster ID

splitRows <- gsub('_.*', '', rownames(counts.tr))


# Splitting the data frame according to the Splitrows vector

counts.split <- split.data.frame(counts.tr,
                                 f = factor(splitRows))
# Fixing cloumn names and transpose matrix
# gsub('.*_(.*)', '\\1', '0_HEThet1')

counts.split.modified <- lapply(counts.split, function(x){
  
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})


# Running the deSEQ analysis on the data with batch number as a covariate
# Get counts matrix

########### CLUSTER 0 ############

counts_0 <- counts.split.modified$`0`

# Generate sample level meta data

colData <- data.frame(Sample = colnames(counts_0))

colData <- colData %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

Batch = (factor(c("Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_1", "Batch_1", "Batch_1", 
                  "Batch_2", "Batch_2", "Batch_2", "Batch_1", "Batch_1", "Batch_1")))

colData <- cbind(colData, Batch)

colData$Condition<-factor(colData$Condition)
colData$Condition<-relevel(colData$Condition,ref="WT")

# Create DESeq2 object

deseq2_0 <- DESeqDataSetFromMatrix(countData = counts_0,
                                   colData = colData,
                                   design = ~ Condition+Batch)


# Filtering the DESeq2 dataset

keep_0 <- rowSums(counts(deseq2_0)) >= 10

deseq2_0 <- deseq2_0[keep_0,]

# Run DESeq2

deseq2_0 <- DESeq(deseq2_0)

# Check coefficients for comparison

resultsNames(deseq2_0)

# Generate result

res_0 <- results(deseq2_0, name = "Condition_HET_vs_WT")

res_0

########### CLUSTER 1 ############

counts_1 <- counts.split.modified$`1`

# Generate sample level meta data

colData_1 <- data.frame(Sample = colnames(counts_1))

colData_1 <- colData_1 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_1 <- cbind(colData_1, Batch)

colData_1$Condition<-factor(colData_1$Condition)
colData_1$Condition<-relevel(colData_1$Condition,ref="WT")


# Create DESeq2 object

deseq2_1 <- DESeqDataSetFromMatrix(countData = counts_1,
                                   colData = colData_1,
                                   design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_1 <- rowSums(counts(deseq2_1)) >= 10

deseq2_1 <- deseq2_1[keep_1,]

# Run DESeq2

deseq2_1 <- DESeq(deseq2_1)

# Check coefficients for comparison

resultsNames(deseq2_1)

# Generate result

res_1 <- results(deseq2_1, name = "Condition_HET_vs_WT")

res_1

########### CLUSTER 2 ############

counts_2 <- counts.split.modified$`2`

# Generate sample level meta data

colData_2 <- data.frame(Sample = colnames(counts_2))

colData_2 <- colData_2 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_2 <- cbind(colData_2, Batch)

colData_2$Condition<-factor(colData_2$Condition)
colData_2$Condition<-relevel(colData_2$Condition,ref="WT")


# Create DESeq2 object

deseq2_2 <- DESeqDataSetFromMatrix(countData = counts_2,
                                   colData = colData_2,
                                   design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_2 <- rowSums(counts(deseq2_2)) >= 10

deseq2_2 <- deseq2_2[keep_2,]

# Run DESeq2

deseq2_2 <- DESeq(deseq2_2)

# Check coefficients for comparison

resultsNames(deseq2_2)

# Generate result

res_2 <- results(deseq2_2, name = "Condition_HET_vs_WT")

res_2

########### CLUSTER 2 ############

counts_3 <- counts.split.modified$`3`

# Generate sample level meta data

colData_3 <- data.frame(Sample = colnames(counts_3))

colData_3 <- colData_3 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_3 <- cbind(colData_3, Batch)

colData_3$Condition<-factor(colData_3$Condition)
colData_3$Condition<-relevel(colData_3$Condition,ref="WT")


# Create DESeq2 object

deseq2_3 <- DESeqDataSetFromMatrix(countData = counts_3,
                                   colData = colData_3,
                                   design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_3 <- rowSums(counts(deseq2_3)) >= 10

deseq2_3 <- deseq2_3[keep_3,]

# Run DESeq2

deseq2_3 <- DESeq(deseq2_3)

# Check coefficients for comparison

resultsNames(deseq2_3)

# Generate result

res_3 <- results(deseq2_3, name = "Condition_HET_vs_WT")

res_3

########### CLUSTER 4 ############

counts_4 <- counts.split.modified$`4`

# Generate sample level meta data

colData_4 <- data.frame(Sample = colnames(counts_4))

colData_4 <- colData_4 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_4 <- cbind(colData_4, Batch)

colData_4$Condition<-factor(colData_4$Condition)
colData_4$Condition<-relevel(colData_4$Condition,ref="WT")


# Create DESeq2 object

deseq2_4 <- DESeqDataSetFromMatrix(countData = counts_4,
                                   colData = colData_4,
                                   design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_4 <- rowSums(counts(deseq2_4)) >= 10

deseq2_4 <- deseq2_4[keep_4,]

# Run DESeq2

deseq2_4 <- DESeq(deseq2_4)

# Check coefficients for comparison

resultsNames(deseq2_4)

# Generate result

res_4 <- results(deseq2_4, name = "Condition_HET_vs_WT")

res_4

########### CLUSTER 5 ############

counts_5 <- counts.split.modified$`5`

# Generate sample level meta data

colData_5 <- data.frame(Sample = colnames(counts_5))

colData_5 <- colData_5 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_5 <- cbind(colData_5, Batch)

colData_5$Condition<-factor(colData_5$Condition)
colData_5$Condition<-relevel(colData_5$Condition,ref="WT")


# Create DESeq2 object

deseq2_5 <- DESeqDataSetFromMatrix(countData = counts_5,
                                   colData = colData_5,
                                   design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_5 <- rowSums(counts(deseq2_5)) >= 10

deseq2_5 <- deseq2_5[keep_5,]

# Run DESeq2

deseq2_5 <- DESeq(deseq2_5)

# Check coefficients for comparison

resultsNames(deseq2_5)

# Generate result

res_5 <- results(deseq2_5, name = "Condition_HET_vs_WT")

res_5

########### CLUSTER 6 ############

counts_6 <- counts.split.modified$`6`

# Generate sample level meta data

colData_6 <- data.frame(Sample = colnames(counts_6))

colData_6 <- colData_6 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_6 <- cbind(colData_6, Batch)

colData_6$Condition<-factor(colData_6$Condition)
colData_6$Condition<-relevel(colData_6$Condition,ref="WT")


# Create DESeq2 object

deseq2_6 <- DESeqDataSetFromMatrix(countData = counts_6,
                                   colData = colData_6,
                                   design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_6 <- rowSums(counts(deseq2_6)) >= 10

deseq2_6 <- deseq2_6[keep_6,]

# Run DESeq2

deseq2_6 <- DESeq(deseq2_6)

# Check coefficients for comparison

resultsNames(deseq2_6)

# Generate result

res_6 <- results(deseq2_6, name = "Condition_HET_vs_WT")

res_6

########### CLUSTER 7 ############

counts_7 <- counts.split.modified$`7`

# Generate sample level meta data

colData_7 <- data.frame(Sample = colnames(counts_7))

colData_7 <- colData_7 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_7 <- cbind(colData_7, Batch)

colData_7$Condition<-factor(colData_7$Condition)
colData_7$Condition<-relevel(colData_7$Condition,ref="WT")


# Create DESeq2 object

deseq2_7 <- DESeqDataSetFromMatrix(countData = counts_7,
                                   colData = colData_7,
                                   design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_7 <- rowSums(counts(deseq2_7)) >= 10

deseq2_7 <- deseq2_7[keep_7,]

# Run DESeq2

deseq2_7 <- DESeq(deseq2_7)

# Check coefficients for comparison

resultsNames(deseq2_7)

# Generate result

res_7 <- results(deseq2_7, name = "Condition_HET_vs_WT")

res_7
########### CLUSTER 8 ############

counts_8 <- counts.split.modified$`8`

# Generate sample level meta data

colData_8 <- data.frame(Sample = colnames(counts_8))

colData_8 <- colData_8 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_8 <- cbind(colData_8, Batch)

colData_8$Condition<-factor(colData_8$Condition)
colData_8$Condition<-relevel(colData_8$Condition,ref="WT")


# Create DESeq2 object

deseq2_8 <- DESeqDataSetFromMatrix(countData = counts_8,
                                   colData = colData_8,
                                   design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_8 <- rowSums(counts(deseq2_8)) >= 10

deseq2_8 <- deseq2_8[keep_8,]

# Run DESeq2

deseq2_8 <- DESeq(deseq2_8)

# Check coefficients for comparison

resultsNames(deseq2_8)

# Generate result

res_8 <- results(deseq2_8, name = "Condition_HET_vs_WT")

res_8

########### CLUSTER 9 ############

counts_9 <- counts.split.modified$`9`

# Generate sample level meta data

colData_9 <- data.frame(Sample = colnames(counts_9))

colData_9 <- colData_9 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_9 <- cbind(colData_9, Batch)

colData_9$Condition<-factor(colData_9$Condition)
colData_9$Condition<-relevel(colData_9$Condition,ref="WT")


# Create DESeq2 object

deseq2_9 <- DESeqDataSetFromMatrix(countData = counts_9,
                                   colData = colData_9,
                                   design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_9 <- rowSums(counts(deseq2_9)) >= 10

deseq2_9 <- deseq2_9[keep_9,]

# Run DESeq2

deseq2_9 <- DESeq(deseq2_9)

# Check coefficients for comparison

resultsNames(deseq2_9)

# Generate result

res_9 <- results(deseq2_9, name = "Condition_HET_vs_WT")

res_9

########### CLUSTER 10 ############

counts_10 <- counts.split.modified$`10`

# Generate sample level meta data

colData_10 <- data.frame(Sample = colnames(counts_10))

colData_10 <- colData_10 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_10 <- cbind(colData_10, Batch)

colData_10$Condition<-factor(colData_10$Condition)
colData_10$Condition<-relevel(colData_10$Condition,ref="WT")


# Create DESeq2 object

deseq2_10 <- DESeqDataSetFromMatrix(countData = counts_10,
                                    colData = colData_10,
                                    design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_10 <- rowSums(counts(deseq2_10)) >= 10

deseq2_10 <- deseq2_10[keep_10,]

# Run DESeq2

deseq2_10 <- DESeq(deseq2_10)

# Check coefficients for comparison

resultsNames(deseq2_10)

# Generate result

res_10 <- results(deseq2_10, name = "Condition_HET_vs_WT")

res_10

########### CLUSTER 11 ############

counts_11 <- counts.split.modified$`11`

# Generate sample level meta data

colData_11 <- data.frame(Sample = colnames(counts_11))

colData_11 <- colData_11 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_11 <- cbind(colData_11, Batch)

colData_11$Condition<-factor(colData_11$Condition)
colData_11$Condition<-relevel(colData_11$Condition,ref="WT")


# Create DESeq2 object

deseq2_11 <- DESeqDataSetFromMatrix(countData = counts_11,
                                    colData = colData_11,
                                    design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_11 <- rowSums(counts(deseq2_11)) >= 10

deseq2_11 <- deseq2_11[keep_11,]

# Run DESeq2

deseq2_11 <- DESeq(deseq2_11)

# Check coefficients for comparison

resultsNames(deseq2_11)

# Generate result

res_11 <- results(deseq2_11, name = "Condition_HET_vs_WT")

res_11

########### CLUSTER 12 ############

counts_12 <- counts.split.modified$`12`

# Generate sample level meta data

colData_12 <- data.frame(Sample = colnames(counts_12))

colData_12 <- colData_12 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

Batch_12 = (factor(c("Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_1", "Batch_1", 
                     "Batch_2", "Batch_2", "Batch_2", "Batch_1", "Batch_1", "Batch_1")))

colData_12 <- cbind(colData_12, Batch_12)

colData_12$Condition<-factor(colData_12$Condition)
colData_12$Condition<-relevel(colData_12$Condition,ref="WT")


# Create DESeq2 object

deseq2_12 <- DESeqDataSetFromMatrix(countData = counts_12,
                                    colData = colData_12,
                                    design = ~ Condition+Batch_12)

# Filtering the DESeq2 dataset

keep_12 <- rowSums(counts(deseq2_12)) >= 10

deseq2_12 <- deseq2_12[keep_12,]

# Run DESeq2

deseq2_12 <- DESeq(deseq2_12)

# Check coefficients for comparison

resultsNames(deseq2_12)

# Generate result

res_12 <- results(deseq2_12, name = "Condition_HET_vs_WT")

res_12

########### CLUSTER 13 ############

counts_13 <- counts.split.modified$`13`

# Generate sample level meta data

colData_13 <- data.frame(Sample = colnames(counts_13))

colData_13 <- colData_13 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_13 <- cbind(colData_13, Batch)

colData_13$Condition<-factor(colData_13$Condition)
colData_13$Condition<-relevel(colData_13$Condition,ref="WT")


# Create DESeq2 object

deseq2_13 <- DESeqDataSetFromMatrix(countData = counts_13,
                                    colData = colData_13,
                                    design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_13 <- rowSums(counts(deseq2_13)) >= 10

deseq2_13 <- deseq2_13[keep_13,]

# Run DESeq2

deseq2_13 <- DESeq(deseq2_13)

# Check coefficients for comparison

resultsNames(deseq2_13)

# Generate result

res_13 <- results(deseq2_13, name = "Condition_HET_vs_WT")

res_13

########### CLUSTER 14 ############

counts_14 <- counts.split.modified$`14`

# Generate sample level meta data

colData_14 <- data.frame(Sample = colnames(counts_14))

colData_14 <- colData_14 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_14 <- cbind(colData_14, Batch_12)

colData_14$Condition<-factor(colData_14$Condition)
colData_14$Condition<-relevel(colData_14$Condition,ref="WT")


# Create DESeq2 object

deseq2_14 <- DESeqDataSetFromMatrix(countData = counts_14,
                                    colData = colData_14,
                                    design = ~ Condition+Batch_12)

# Filtering the DESeq2 dataset

keep_14 <- rowSums(counts(deseq2_14)) >= 10

deseq2_14 <- deseq2_14[keep_14,]

# Run DESeq2

deseq2_14 <- DESeq(deseq2_14)

# Check coefficients for comparison

resultsNames(deseq2_14)

# Generate result

res_14 <- results(deseq2_14, name = "Condition_HET_vs_WT")

res_14

########### CLUSTER 15 ############

counts_15 <- counts.split.modified$`15`

# Generate sample level meta data

colData_15 <- data.frame(Sample = colnames(counts_15))

colData_15 <- colData_15 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_15 <- cbind(colData_15, Batch)

colData_15$Condition<-factor(colData_15$Condition)
colData_15$Condition<-relevel(colData_15$Condition,ref="WT")


# Create DESeq2 object

deseq2_15 <- DESeqDataSetFromMatrix(countData = counts_15,
                                    colData = colData_15,
                                    design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep_15 <- rowSums(counts(deseq2_15)) >= 10

deseq2_15 <- deseq2_15[keep_15,]

# Run DESeq2

deseq2_15 <- DESeq(deseq2_15)

# Check coefficients for comparison

resultsNames(deseq2_15)

# Generate result

res_15 <- results(deseq2_15, name = "Condition_HET_vs_WT")

res_15

########### CLUSTER 16 ############

counts_16 <- counts.split.modified$`16`

# Generate sample level meta data

colData_16 <- data.frame(Sample = colnames(counts_16))

colData_16 <- colData_16 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

Batch_16 = (factor(c("Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_1", 
                     "Batch_2", "Batch_2", "Batch_2", "Batch_1", "Batch_1", "Batch_1")))


colData_16 <- cbind(colData_16, Batch_16)

colData_16$Condition<-factor(colData_16$Condition)
colData_16$Condition<-relevel(colData_16$Condition,ref="WT")


# Create DESeq2 object

deseq2_16 <- DESeqDataSetFromMatrix(countData = counts_16,
                                    colData = colData_16,
                                    design = ~ Condition+Batch_16)

# Filtering the DESeq2 dataset

keep_16 <- rowSums(counts(deseq2_16)) >= 10

deseq2_16 <- deseq2_16[keep_16,]

# Run DESeq2

deseq2_16 <- DESeq(deseq2_16)

# Check coefficients for comparison

resultsNames(deseq2_16)

# Generate result

res_16 <- results(deseq2_16, name = "Condition_HET_vs_WT")

res_16

########### CLUSTER 17 ############

counts_17 <- counts.split.modified$`17`

# Generate sample level meta data

colData_17 <- data.frame(Sample = colnames(counts_17))

colData_17 <- colData_17 %>% 
  mutate(Condition =  ifelse(grepl('HET', Sample), 'HET', 
                             ifelse(grepl('HOM', Sample), 'HOM', 'WT'))) %>%
  column_to_rownames(var = 'Sample')

colData_17 <- cbind(colData_17, Batch_12)

colData_17$Condition<-factor(colData_17$Condition)
colData_17$Condition<-relevel(colData_17$Condition,ref="WT")


# Create DESeq2 object

deseq2_17 <- DESeqDataSetFromMatrix(countData = counts_17,
                                    colData = colData_17,
                                    design = ~ Condition+Batch_12)

# Filtering the DESeq2 dataset

keep_17 <- rowSums(counts(deseq2_17)) >= 10

deseq2_17 <- deseq2_17[keep_17,]

# Run DESeq2

deseq2_17 <- DESeq(deseq2_17)

# Check coefficients for comparison

resultsNames(deseq2_17)

# Generate result

res_17 <- results(deseq2_17, name = "Condition_HET_vs_WT")

res_17



# Save results

write.csv(res_0, file="Cluster_0_DE_results_unfiltered.csv")
write.csv(res_1, file="Cluster_1_DE_results_unfiltered.csv")
write.csv(res_2, file="Cluster_2_DE_results_unfiltered.csv")
write.csv(res_3, file="Cluster_3_DE_results_unfiltered.csv")
write.csv(res_4, file="Cluster_4_DE_results_unfiltered.csv")
write.csv(res_5, file="Cluster_5_DE_results_unfiltered.csv")
write.csv(res_6, file="Cluster_6_DE_results_unfiltered.csv")
write.csv(res_7, file="Cluster_7_DE_results_unfiltered.csv")
write.csv(res_8, file="Cluster_8_DE_results_unfiltered.csv")
write.csv(res_9, file="Cluster_9_DE_results_unfiltered.csv")
write.csv(res_10, file="Cluster_10_DE_results_unfiltered.csv")
write.csv(res_11, file="Cluster_11_DE_results_unfiltered.csv")
write.csv(res_12, file="Cluster_12_DE_results_unfiltered.csv")
write.csv(res_13, file="Cluster_13_DE_results_unfiltered.csv")
write.csv(res_14, file="Cluster_14_DE_results_unfiltered.csv")
write.csv(res_15, file="Cluster_15_DE_results_unfiltered.csv")
write.csv(res_16, file="Cluster_16_DE_results_unfiltered.csv")
write.csv(res_17, file="Cluster_17_DE_results_unfiltered.csv")


########### ALL CLUSTERS ############


counts<-counts_clusters_merged$RNA
colData <- data.frame(as.factor(c(rep("HET",9),rep("WT",6))))

#colData <-  colnames(counts)

#colData <- colData %>% 
# mutate(Condition =  ifelse(grepl('HET', Condition), 'HET', 
#                            ifelse(grepl('HOM', Condition), 'HOM', 'WT'))) %>%
# column_to_rownames(var = 'Condition')
#View(colData)


rownames(colData)<-colnames(counts)

colnames(colData)<-"Condition"

#colData<-c(rep("HET",3),rep("WT",3))

Batch = (factor(c("Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_2", "Batch_1", "Batch_1", "Batch_1", 
                  "Batch_2", "Batch_2", "Batch_2", "Batch_1", "Batch_1", "Batch_1")))

colData <- cbind(colData, Batch)

colData$Condition<-factor(colData$Condition)
colData$Condition<-relevel(colData$Condition,ref="WT")



# Create DESeq2 object

deseq2 <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = colData,
                                 design = ~ Condition+Batch)

# Filtering the DESeq2 dataset

keep <- rowSums(counts(deseq2)) >= 10

deseq2 <- deseq2[keep,]

# Run DESeq2

deseq2 <- DESeq(deseq2)

# Check coefficients for comparison

resultsNames(deseq2)

# Generate result

res <- results(deseq2, name = "Condition_HET_vs_WT")

res


write.csv(res, file="Clusters_Neuronal_DE_results_unfiltered.csv")


# Making dot plots for gene marker expression

levels(data)

colnames(data@meta.data)

Idents(data) <- "seurat_clusters"

markers.to.plot <- c("TTC6", "RGS6", "FOXP2", "MDGA2", "NRXN3", "NRG1", "CRABP1", "TAGLN",
                     "PFN1", "MID1", "LARGE1", "LRRC7", "RBFOX1", "NTM", "PDE4B", "CNTNAP5",
                     "CTNNA2", "COPS9", "UQCR11", "COX8A", "DYNLRB1", "MRFAP1",
                     "INPP4B", "ANKS1B", "RGS7")


DotPlot(data, features = markers.to.plot, cols = c("blue", "red", "yellow", "green", "purple", "orange", "pink", "brown", "cyan", "black"), dot.scale = 10, split.by = "seurat_clusters", c(4, 5, 6, 9)) +
  RotatedAxis()

DotPlot(data, features = markers.to.plot, cols = c("blue", "red", "purple", "green", "yellow", "orange", "cyan"), dot.scale = 10, idents =  c(4, 5, 6, 11, 13, 14, 15)) +
  RotatedAxis()

