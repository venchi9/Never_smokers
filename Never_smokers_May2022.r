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

#Look at the cell numbers in each condition
#Create an object for non-smokers and for 1k1k
Idents(sct)<-"orig.ident"
NS<-subset(sct, idents="Non-smokers")
onek<-subset(sct, idents = "onek1k")

#Cell numbers for never smokers - n=67092
table(NS$predicted.celltype.l2)
             ASDC    B intermediate          B memory           B naive 
               26               729               968              4246 
        CD14 Mono         CD16 Mono           CD4 CTL         CD4 Naive 
             9152               986              2191              6599 
CD4 Proliferating           CD4 TCM           CD4 TEM         CD8 Naive 
               56             10536              2000              1309 
CD8 Proliferating           CD8 TCM           CD8 TEM              cDC1 
               21               170              9423                33 
             cDC2               dnT           Doublet             Eryth 
              599               122                47                 7 
              gdT              HSPC               ILC              MAIT 
              934                82                41               864 
               NK  NK Proliferating     NK_CD56bright               pDC 
             8303                56               353               141 
      Plasmablast          Platelet              Treg 
              167              5406              1525 

write.csv((table(NS$predicted.celltype.l2)), file = "never_smokers_celltypes.csv")

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