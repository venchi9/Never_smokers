#Adapted from Himanshi's code to extract the median expression and pvalues.  The loop doesn't work but you can use it step-by-step

library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(tidyr)

# clusters<- c("Plasmablast","HSPC","cDC","CD14-Mono","Platelet","B-memory","CD4-TCM","CD8-Naive","CD8-TEM","dnT","pDC","CD4-CTL","CD16-Mono","CD4-TEM","CD4-Naive","B-intermediate","B-naive","MAIT","NK","gdT","CD8-TCM","Treg","NK_CD56bright")
# clusters <- c("Plasmablast","HSPC","cDC","CD14-Mono","Platelet",)
clusters <- c("cDC","CD14-Mono","CD16-Mono")
chroms <- "2"
condition= "group"

work_dir <- "/directflow/SCCGGroupShare/projects/himaro/projects/never_smoking_lungcancer/eqtl_analysis"
input_dir <-"/directflow/SCCGGroupShare/projects/himaro/projects/never_smoking_lungcancer/eqtl_analysis/eqtl_analysis_output/reanalyse_matrix_eqtl/correct_filter_after_interaction/inflammation_filter/correct_crossmatrixEQTL_inflammation"
setwd(work_dir)

for (cluster in clusters){
for (chrom in chroms){

            exprs_filename <- sprintf("%s/prepGenotype_files/output/CountMatrix/%s/%s/zscores/nvrsmoke_onek1k_%s_%s_Chr%s_ZscoreMatrix.tsv", work_dir, condition, cluster, condition, cluster, chrom)

            genotype_filename <- sprintf("%s/prepGenotype_files/inputs/GeneLoc/%s/%s/genotype_file/%s_%s_Chr%s_SNPs.tsv", work_dir, condition, cluster, condition, cluster, chrom)

            Celltype_interaction_test <- sprintf("%s/%s/%s/nvrsmoke_onek1k_%s_%s_cisEQTLs_combined_Inflammationfiltered.tsv", input_dir, condition, cluster, condition, cluster)
            
            flag <- TRUE
            tryCatch(Celltype_interaction_test <- read_tsv(Celltype_interaction_test), error=function(e) flag<<-FALSE)
            if (!flag) next

            ## Load gene expression data
            gene_mat <- read_tsv(exprs_filename)
            gene_mat <- as.data.frame(gene_mat)
            names(gene_mat)[names(gene_mat) == "gene_id"] <- "ensembl_gene_id"
            rownames(gene_mat) <- gene_mat[, "ensembl_gene_id"]

            Celltype_interaction_expression <- merge(x=Celltype_interaction_test, y=gene_mat, by = 'ensembl_gene_id')
            # Celltype_interaction_expression <- Celltype_interaction_expression[ -c(4,6,7,8) ]
            as_tibble(Celltype_interaction_expression)

            Celltype_interaction_expression <- Celltype_interaction_expression %>% pivot_longer(cols = starts_with("0_"),
                              names_to= c("individuals"),
                              values_to='expression')


            # Load snps data
            ## REF -left one    ALT - right one

            snp_mat <- read_tsv(genotype_filename)
            snp_mat <- as.data.frame(snp_mat)
            rownames(snp_mat) <- snp_mat$id
            snps <- snp_mat %>% separate(col = id, into = c("CHROM", "POS", "REF", "ALT"), sep = ":")
     
            snps_drop <- snps
            col <- colnames(snps_drop)

            snps_drop$snps <- rownames(snps_drop)

            snp_data <- snps_drop %>% pivot_longer(cols = starts_with("0_"),
                              names_to= c("individuals"),
                              values_to='allel')

            Celltype_interaction_expression_snps <- merge(x=Celltype_interaction_expression, y= snp_data, by=c("individuals","snps")) 
            Celltype_interaction_expression_snps <- data.table(Celltype_interaction_expression_snps)
            Celltype_interaction_expression_snps <- Celltype_interaction_expression_snps %>% drop_na()



            Celltype_interaction_expression_snps$group <- ifelse(grepl("_C",Celltype_interaction_expression_snps$individuals), "non_smoker", 
            ifelse(grepl("_LB",Celltype_interaction_expression_snps$individuals), "non_smoker", 
            ifelse(grepl("_SV",Celltype_interaction_expression_snps$individuals), "non_smoker", "onek1k")))

            unique(Celltype_interaction_expression_snps[group=="onek1k"]$individuals)
            unique(Celltype_interaction_expression_snps[group=="non_smoker"]$individuals)
      
            #Celltype_interaction_expression_snps_new[group=="onek1k" & individuals=="0_LB037-C1"]


            ###########################################################################################################################
          print("Plotting the graph ")
     

            dir.create(paste0("/directflow/SCCGGroupShare/projects/venchi/Lung/Never_smokers/box_plots_data/",cluster))
            setwd(paste0("/directflow/SCCGGroupShare/projects/venchi/Lung/Never_smokers/box_plots_data/",cluster))


           


            eQTLbetaFigures <- list()
            for (gene_sym in unique(Celltype_interaction_expression_snps$gene_symbol)){
            ensembl_gene_id <- unique(Celltype_interaction_expression_snps[gene_symbol == gene_sym,]$ensembl_gene_id)
            for (SNP in unique(Celltype_interaction_expression_snps[which(Celltype_interaction_expression_snps$gene_symbol == gene_sym),]$snps)){
                  REF <- unique(Celltype_interaction_expression_snps[which(Celltype_interaction_expression_snps$snps == SNP),]$REF)
                  ALT <- unique(Celltype_interaction_expression_snps[which(Celltype_interaction_expression_snps$snps == SNP),]$ALT)
                  eQTLbetaFigures[[paste0(gene_sym,SNP)]] <- ggplot(Celltype_interaction_expression_snps[snps == SNP & gene_symbol == gene_sym], aes(factor(allel,levels=c("0","1","2")), expression, color = factor(group))) +
                        theme_classic() + 
                        geom_boxplot(aes(group = interaction(factor(allel,levels=c("0","1","2")), factor(group))),color = "black", width = 0.5, position=position_dodge2(preserve = "single"), fill = "white", outlier.shape = NA, lwd=0.5) +
                        geom_point(aes(color = factor(group)), position=position_jitterdodge(dodge.width = 0.7), alpha = 0.80, size = 1) +  #geom_bar(stat = "doge") + 
                    
                        
                              scale_color_manual(values = c(non_smoker = "#e0a307", onek1k = '#6ca3d4')) +
                              theme(legend.position = "none") +
                              labs(x=NULL, y=NULL) +
                               theme(plot.title = element_text(hjust = 0.5), 
                              #text = element_text(size = 6), 
                               plot.subtitle = element_text(hjust = 0.5),
                               plot.caption =  element_text(hjust = 0)) +
                              #legend.key.size = unit(4, "mm"),
                              #legend.key.width = unit(2,"mm")) +
                          ylab(paste0(gene_sym, " Normalized Expression")) +
                          xlab(paste0(SNP, " SNP")) +
                          scale_x_discrete(labels=c("0" = paste0(REF,"/",REF), "1" = paste0(REF,"/",ALT), "2" = paste0(ALT,"/",ALT))) +
                          geom_smooth(position = "identity", aes(group = factor(group)), method = "lm", se = FALSE, size = 0.35, fullrange=T, color="red")



                  ggsave(paste0("eQTLbetaFigures_",gene_sym,"_",SNP,".png"), eQTLbetaFigures[[paste0(gene_sym,SNP)]],width = 6, height = 5)


setwd(paste0("/directflow/SCCGGroupShare/projects/venchi/Lung/Never_smokers/box_plots_data/",cluster))


dt <- as.data.table(Celltype_interaction_expression_snps)

# Median expression per gene × SNP × genotype × group
median_summary <- dt[, .(
  median_expression = median(expression, na.rm = TRUE)
), by = .(gene_symbol, snps, allel, group)]

# Single p-value across genotypes for each gene × SNP × group
pval_summary <- dt[, .(
  p_value = tryCatch({
    expr_list <- split(expression, allel)
    kruskal.test(expr_list)$p.value
  }, error = function(e) NA)
), by = .(gene_symbol, snps, group)]

# Merge medians and p-values
summary_table <- merge(median_summary, pval_summary, by = c("gene_symbol", "snps", "group"))

write.csv(summary_table, "eQTL_expression_summary_by_group.csv", row.names = FALSE)

      }
            }
  

      
      
}
}