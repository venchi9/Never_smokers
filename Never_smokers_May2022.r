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
