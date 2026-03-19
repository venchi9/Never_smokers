#Create heatmap of the significant eqtls within the never smoking data
#Data is here:  /directflow/SCCGGroupShare/projects/venchi/Lung/Never_smokers/eQTL/final_interaction_data

#Load packages
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ComplexHeatmap)

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/Lung/Never_smokers/eQTL/final_interaction_data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/Lung/Never_smokers/eQTL/figures"

#Load in data
x<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/Lung/Never_smokers/eQTL/final_interaction_data/merged_data.csv")

pdf(str_glue("{fig_dir}/heatmap_eQTLs_by_celltype_statistic.pdf"), width=11.6, height=8.2)
p<-ggplot(x, aes(x = gene_symbol, y = cluster, fill = statistic)) + geom_tile(color = "grey")+ coord_fixed()+ scale_fill_gradient2(low = "blue", mid = "white", high = "red")
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(str_glue("{fig_dir}/heatmap_eQTLs_by_celltype_beta.pdf"), width=11.6, height=8.2)
p<-ggplot(x, aes(x = gene_symbol, y = cluster, fill = beta)) + geom_tile(color = "grey")+ coord_fixed()+ scale_fill_gradient2(low = "blue", mid = "white", high = "red")
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#Ordered by beta
a<-x[order(x$beta),]
x<-x[order(x$beta),]

pdf(str_glue("{fig_dir}/heatmap_eQTLs_by_celltype_beta_test.pdf"), width=11.6, height=8.2)
p<-ggplot(x, aes(x = gene_symbol, y = cluster, fill = beta)) + geom_tile(color = "grey")+ coord_fixed()+ scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-4, 11))
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



#Plot the number of eQTLs by celltype
y<-table(x$cluster)
y<-as.data.frame(y)
y<-y[order(-y$Freq),]
#rename columns
colnames(y)[1]<-"Celltype"
colnames(y)[2]<-"eQTL_Frequency"

#Plot
pdf(str_glue("{fig_dir}/Barplot_eQTLs_by_celltype.pdf"), width=11.6, height=8.2)
barplot(height = y$eQTL_Frequency, names = y$Celltype, col = rainbow(22), las = 2, cex.names=0.7) 
dev.off()

#List of celltypes in order of frequency of eQTLs
cell_level_order<-c("CD14-Mono"  ,    "CD16-Mono"  ,    "CD4-Naive"   ,   "CD4-TCM"  ,     
"cDC"     ,       "B-memory"   ,    "CD8-TCM"   ,     "B-naive"   ,    
"CD8-Naive"  ,    "B-intermediate","CD8-TEM"   ,    "gdT"  ,         
 "NK"      ,       "MAIT"    ,       "NK_CD56bright" , "CD4-TEM"   ,    
 "CD4-CTL"    ,    "Treg"     ,      "dnT"  ,          "HSPC")   


 #Now do the same for the genes
y<-table(x$gene_symbol)
y<-as.data.frame(y)
y<-y[order(-y$Freq),]
colnames(y)[1]<-"Gene"
colnames(y)[2]<-"eQTL_Frequency"

#Plot
pdf(str_glue("{fig_dir}/Barplot_eQTLs_by_gene.pdf"), width=11.6, height=8.2)
barplot(height = y$eQTL_Frequency, names = y$Gene, col = rainbow(22), las = 2, cex.names=0.7, ylim = c(0, 20))
dev.off()

#Get a prioritised list of genes
"BST2"  ,   "IL7R" ,    "EIF2AK2",  "PDE4B"  ,  "IRF1"  ,   "LY6E"   , 
"EMP3" ,   "IRF7"  ,   "TNFSF10" , "BTG2"  ,   "CCL5"     "GPR183"  
"CD55"     "ATP2B1"   "NFKB1"    "NFKBIA"   "TNFRSF1B" "CCR7"    
 "CD69"     "NMI"      "SEMA4D"   "HIF1A"    "IFNAR1"   "IL4R"    
"ABI1"     "AHR"      "CD82"     "FPR1"     "IL1B"     "IL2RB"   
"KLF6"     "LYN"      "MYC"      "NLRP3"    "RASGRP1"  "RGS1"    
 "TAPBP"    "C5AR1"    "EREG"     "FFAR2"    "GCH1"     "HBEGF"   
"IFITM1"   "IL18"     "LCP2"     "MSR1"     "OLR1"     "PLAUR"   
"PTGER2"   "PTGER4"   "PTPRE"    "SELL"     "TLR2"     "ADM"     
 "ATP2A2"   "C3AR1"    "CD40"     "CD48"     "CDKN1A"   "GNAI3"   
 "ICAM1"    "IL10RA"   "IL18RAP"  "LDLR"     "MARCO"    "MEFV"    
 "NOD2"     "PIK3R5"   "PSEN1"    "PTGIR"    "RELA"     "RHOG"    
 "RIPK2"    "RNF144B"  "SGMS2"    "TLR1"     "TNFSF9"


 #Try and cluster the heatmap based on ordered lists
pdf(str_glue("{fig_dir}/heatmap_eQTLs_by_celltype_beta_celltypes_ordered.pdf"), width=11.6, height=8.2)
p<-ggplot(x, aes(x = gene_symbol, y = factor(cluster, level = cell_level_order), fill = beta)) + geom_tile(color = "grey")+ coord_fixed()+ scale_fill_gradient2(low = "blue", mid = "white", high = "red")
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
#Order both
pdf(str_glue("{fig_dir}/heatmap_eQTLs_by_celltype_beta_celltypes_genes_ordered.pdf"), width=11.6, height=8.2)
p<-ggplot(x, aes(x = factor(gene_symbol, level = gene_level_order), y = factor(cluster, level = cell_level_order), fill = beta)) + geom_tile(color = "grey")+ coord_fixed()+ scale_fill_gradient2(low = "blue", mid = "white", high = "red")
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#We actually have some quite high false discovery rates.  If we limit this to a FDR of 5%, let's see how this goes
df <-filter(x, FDR < 0.05)

#Do everything again


#Plot the number of eQTLs by celltype
y<-table(df$cluster)
y<-as.data.frame(y)
y<-y[order(-y$Freq),]
#rename columns
colnames(y)[1]<-"Celltype"
colnames(y)[2]<-"eQTL_Frequency"
cell_level_order<-as.vector(y$Celltype)
#Plot
pdf(str_glue("{fig_dir}/Barplot_eQTLs_by_celltype_filteredFDR.pdf"), width=11.6, height=8.2)
barplot(height = y$eQTL_Frequency, names = y$Celltype, col = rainbow(22), las = 2, cex.names=0.7) 
dev.off()


 #Now do the same for the genes
y<-table(df$gene_symbol)
y<-as.data.frame(y)
y<-y[order(-y$Freq),]
colnames(y)[1]<-"Gene"
colnames(y)[2]<-"eQTL_Frequency"
gene_level_order<-as.vector(y$Gene)
#Plot
pdf(str_glue("{fig_dir}/Barplot_eQTLs_by_gene_filteredFDR.pdf"), width=11.6, height=8.2)
barplot(height = y$eQTL_Frequency, names = y$Gene, col = rainbow(22), las = 2, cex.names=0.7, ylim = c(0, 20))
dev.off()

x<-df
pdf(str_glue("{fig_dir}/heatmap_eQTLs_by_celltype_beta_celltypes_genes_ordered_filteredFDR.pdf"), width=11.6, height=8.2)
p<-ggplot(x, aes(x = factor(gene_symbol, level = gene_level_order), y = factor(cluster, level = cell_level_order), fill = beta)) + geom_tile(color = "grey")+ coord_fixed()+ scale_fill_gradientn(colours = rev(col),  na.value = "grey", breaks=c(-11,-4, 0, 4, 11),labels=c(-11,-4, 0, 4, 11),
limits=c(-11,11))
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()