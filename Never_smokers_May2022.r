#Looking at the never smoking and OneK1K dataset for the first time.
# Data is here:
/directflow/SCCGGroupShare/projects/venchi/Lung/Raw_Data/never_smoking_lungcancer

#I will create a new folder for results here:

cd /directflow/SCCGGroupShare/projects/venchi/Lung/Never_smokers

#I will also be pushing this to Github

fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate r_base

library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratData)
library(cowplot)
library(Nebulosa)

#This is listed as sct integrated.rds"
lung<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/Lung/Raw_Data/never_smoking_lungcancer/never_smoking_onek1k_sct_integrated.rds")

#This is just "sct.rds"
sct<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/Lung/Raw_Data/never_smoking_lungcancer/never_smoking_onek1k_sct.rds")

#Do a UMAP
pdf("Never_smoker_UMAP.pdf", width=11.6, height=8.2)
DimPlot(sct, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf("Never_smoker_UMAP_splitby_ident.pdf", width=11.6, height=8.2)
DimPlot(sct, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident") + NoLegend()
dev.off()

#Check the %MT in both conditions
pdf("QC_prior_to_filtering.pdf", width=11.6, height=8.2)
VlnPlot(sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), split.by = "orig.ident")
dev.off()

#The non-smokers have not been QCed properly
#The science paper used a fairly complex method of QC:
#For each pool of cells captured, the distributions of total number of UMIs, number
#of genes, and percentage of mitochondrial gene expression were normalized using an
#Ordered Quantile Transformation as follows:
#The percentile of each value in the distribution was calculated and then mapped to
#the same percentile in a normal distribution using the inverse of the standard normal  
5
#cumulative density function. The Z-scores corresponding to three negative units were
#mapped back to the original distributions using a binomial generalized linear model with
#a logit link (90). Cells below the Z-score limits were considered as outliers and removed
#from further analysis. For the expression of mitochondrial genes, an additional upper
#threshold corresponding to two Z-scores was used to exclude apoptotic cells. After
#quality control, 1,272,518 cells were kept and 176,928 were filtered out. The average
#read depth of cells retained for further analysis was 34,000.

#For 1k1k
#min n_feature = 252
#Max feature = 6086
#min counts = 600
#Max counts = 59763
#min MT = 0.0731
#Max MT = 7.832512

#For non-smokers
#min feature = 101
#Max feature = 7176
#min counts = 144
#Max counts = 75735
#Min MT = 0
#max = 10

#Subset the matrix so QC is exactly like the 1k1k dataset
#before filtering
Non-smokers      onek1k 
      67092      337979 
sub<-subset(sct, subset = nFeature_RNA > 251 & nFeature_RNA <6087 & percent.mt <7.832513)

#After filtering
Non-smokers      onek1k 
      60311      337979

#Save filtered object
saveRDS(sub, file = "never_smoking_onek1k_sct_QC_filtered.rds")
#I have spoken with Seyhan - this filtering approach seems valid.
Idents(NS)<-"orig.ident"
pdf("QC_post_to_filtering.pdf", width=11.6, height=8.2)
VlnPlot(NS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), split.by = "orig.ident")
dev.off()

#Look at the cell numbers in each condition
#Create an object for non-smokers and for 1k1k
Idents(NS)<-"orig.ident"
never<-subset(NS, idents="Non-smokers")
onek<-subset(NS, idents = "onek1k")

#Cell numbers for never smokers - n=67092
table(never$predicted.celltype.l2)
             ASDC    B intermediate          B memory           B naive 
               22               671               902              3875 
        CD14 Mono         CD16 Mono           CD4 CTL         CD4 Naive 
             8636               933              2053              6428 
CD4 Proliferating           CD4 TCM           CD4 TEM         CD8 Naive 
               49             10222              1878              1272 
CD8 Proliferating           CD8 TCM           CD8 TEM              cDC1 
               19               166              8777                33 
             cDC2               dnT           Doublet             Eryth 
              578               111                44                 6 
              gdT              HSPC               ILC              MAIT 
              886                81                40               745 
               NK  NK Proliferating     NK_CD56bright               pDC 
             7730                55               336               138 
      Plasmablast          Platelet              Treg 
              153              1995              1477 

write.csv((table(never$predicted.celltype.l2)), file = "never_smokers_celltypes.csv")

#Cell numbers for 1k1k, n= 337979 
table(onek$predicted.celltype.l2)
             ASDC    B intermediate          B memory           B naive 
               47              6974              7568             16341 
        CD14 Mono         CD16 Mono           CD4 CTL         CD4 Naive 
             9627              4365              4448             72467 
CD4 Proliferating           CD4 TCM           CD4 TEM         CD8 Naive 
              201             81507              9454             13405 
CD8 Proliferating           CD8 TCM           CD8 TEM              cDC1 
               73              3973             41859                29 
             cDC2               dnT           Doublet             Eryth 
             1338               538                35                64 
              gdT              HSPC               ILC              MAIT 
             5519               465               131              2765 
               NK  NK Proliferating     NK_CD56bright               pDC 
            43359               466              2082               560 
      Plasmablast          Platelet              Treg 
              810               487              7022 



write.csv((table(onek$predicted.celltype.l2)), file = "Onek_celltypes.csv")

#Do another umap now that I've filtered the cells
Idents(NS)<-"predicted.celltype.l2"
pdf("UMAP_post_filtering.pdf", width=11.6, height=8.2)
DimPlot(NS, reduction = "umap", split.by = "orig.ident",  label = TRUE, repel = TRUE) + NoLegend()
dev.off()

#No split
Idents(NS)<-"orig.ident"
pdf("UMAP_post_filtering_compare_conditions.pdf", width=11.6, height=8.2)
DimPlot(NS, reduction = "umap") 
dev.off()

#I will do DEG based on this filtered matrix using Himanshi's code:
CD4_Proliferating <- subset(diff_exp, idents = "CD4 Proliferating")

de_genes1 <- FindMarkers(CD4_Proliferating, ident.1 = "Non-smokers", ident.2 = "onek1k", group.by = "orig.ident", assay = "RNA", test.use = "wilcox", vars.to.regress = "orig.ident")

de_genes1$genes <- 1:nrow(de_genes1)
#Filter results
de_genes1 <- subset(de_genes1, de_genes1$p_val_adj < 0.05/1520)
de_genes1 <- subset(de_genes1, abs(de_genes1$avg_log2FC) > 1.0)
de_genes1 <- de_genes1 %>% arrange(de_genes1$p_val_adj)

write.table(de_genes1, file = "CD4_Proliferating_wilcox.tsv", sep="\t", quote=F)


#Using Muscat:
#http://www.bioconductor.org/packages/devel/bioc/vignettes/muscat/inst/doc/analysis.html
#First load data as a single cell experiment
sce <- as.SingleCellExperiment(diff_exp)

# calculate per-cell quality control (QC) metrics
library(scater)
qc <- perCellQCMetrics(sce)

dim(sce)
[1]  30655 398290

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

[1]  30655 389417

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

[1]   9571 389417

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

#Data preparation
#muscat expects a certain format of the input SCE. Specifically, the following cell metadata (colData) columns have to be provided:

#"sample_id": unique sample identifiers (e.g., PeterPan_ref1, Nautilus_trt3, ???) = "individual"
#"cluster_id": subpopulation (cluster) assignments (e.g., T cells, monocytes, ???) = "pedicted.celltype.l2"
#"group_id": experimental group/condition (e.g., control/treatment, healthy/diseased, ???) = "orig.ident"

sce
class: SingleCellExperiment 
dim: 9571 389417 
metadata(0):
assays(2): counts logcounts
rownames(9571): NOC2L HES4 ... WASH6P RP11-125B2.1
rowData names(0):
colnames(389417): GCCCAGACAGCGAGTA-1_1 CATAAGCGTCAAAGAT-1_1 ...
  TTTGTCAGTATTACCG-9 TTTGTCATCAGATAAG-9
colData names(16): orig.ident nCount_RNA ... nFeature_SCT ident
reducedDimNames(0):
mainExpName: SCT
altExpNames(2): predicted_ADT RNA

#Data preparation
sce$id<-paste0(sce$orig.ident, sce$individual)
sce<-prepSCE(sce,kid = "predicted.celltype.l2", gid = "orig.ident", sid = "id", drop = TRUE)

#For consistency and easy accession throughout this vignette, we will store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

#Data overview
# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

#Differential State Analysis

#Aggregation of single cell to pseudobulk data
pb <- aggregateData(sce,
    assay = "counts", fun = "sum",
    by = c("cluster_id", "sample_id"))

# one sheet per subpopulation
assayNames(pb)