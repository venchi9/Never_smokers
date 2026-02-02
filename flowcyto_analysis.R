#!/usr/bin/Rscript


#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
#         File: flowcyto_analysis.R
#               Analysis of flowcytometry results
#       Author: Kavitha KS
#   Created on: 13/11/2024
#    Conda env: /directflow/SCCGGroupShare/projects/kavkri/.conda/envs/seurat5
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

'''
BACKGROUND - the cause of never smoking lung cancer is unknown, however is mediated by IL1B expression in monocytes.  
We have found several eQTLs in our cohort in IL1B when we compared patients with NS lung cancer to NS controls.  
IL1B response to pollution exposure is key here.

METHODS - using some of the 10k10k cohort, we have taken PBMCs, split the sample in half, stimulated half with LPS (simulates bacterial infection) 
and then analysed both the control (untreated) and treated (LPS) cells via flow cytometry.  We measured the expression of IL1B and TNFa
Data - we have done these experiments in batches.  Within each batch there is up to 16 individuals.  For each individual you have an 
MFI (mean fluorescence) for both IL1B and TNFa for both control and LPS.  There are 3 SNPs in IL1B and for each SNP a genotype.  
For each individual, the SNP and genotype are also listed.  Not all individuals have all three SNPs.

HYPOTHESIS - that IL1B genotype correlates with the degree of monocyte stimulation and expression of IL1B and TNFa.  There may be some 
genotypes which correlate with high IL1B expression and therefore may predispose a person to an inflammatory response and cancer (and perhaps other diseases)

WHAT I WANT - I want some box plots grouped by SNP and genotype to compare the level of activation to see if there are genotypes which are more "inflammatory" than others.
Specifics - Can we do boxplots, comparing the different SNP genotypes for 
1.  Baseline IL1B and TNFa expression.  
2.  The treated IL1B and TNFa expression 
3.  Percentage change between control and treated (LPS/control x 100)
Special note - RZ749 didnt go so smoothly.  Can you please do the analyses with and without this experiment and we can see how skewed they are.

data location: https://docs.google.com/spreadsheets/d/1N1g9IlRJG_tCQwiQmMEkFWT85eya9A8J/edit?usp=sharing&ouid=108742766295744171557&rtpof=true&sd=true
'

library(dplyr)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(rstatix) 

my_dir <- "/directflow/SCCGGroupShare/projects/kavkri/venessa/NeverSmokingLungCancer/flowcytometry-IL1B-TNF1a"
setwd(my_dir)

#plots_dir <- "/directflow/SCCGGroupShare/projects/kavkri/venessa/NeverSmokingLungCancer/flowcytometry-IL1B-TNF1a/plots"
plots_dir <- "/directflow/SCCGGroupShare/projects/venchi/Lung/Never_smokers/plots"
data_dir <- "/directflow/SCCGGroupShare/projects/kavkri/venessa/NeverSmokingLungCancer/flowcytometry-IL1B-TNF1a/data"

# table <- "RZ748"
# tsv_file <- sprintf("%s/%s_Results_Table.tsv", data_dir, table)
# table_df <- as.data.frame(read.table(tsv_file, header=TRUE, sep="\t"))

# table_df <- table_df %>%
#     rename(Condition = X)

col_names <- c("PatientID", "Position", "Ref", "Alt", "Allele", "Condition", "Live_cells_Count",	"CD3+ Freq. of Parent", "CD14+ Freq. of Parent", "CD14+ Count", "MFI IL1b Monocytes", "MFI TNFa Monocytes")
all_data <- data.frame(matrix(ncol=12, nrow=0))
colnames(all_data) <- col_names

tables <- c("RZ748", "RZ750", "RZ751", "RZ752", "RZ753", "RZ754", "RZ755")
for (table in tables) {
    print(table)
    tsv_file <- sprintf("%s/%s_Results_Table.tsv", data_dir, table)
    table_df <- as.data.frame(read.table(tsv_file, header=TRUE, sep="\t"))
    print(head(table_df))
    colnames(table_df) <- col_names
    all_data <- rbind(all_data, table_df)   
}

dim(all_data)
head(all_data)
tail(all_data)
all_data$Position <- as.factor(all_data$Position)
all_data$Group <- ifelse(grepl('LPS', all_data$Condition) == TRUE, "Treated (LPS)", "Control")
csv_file <- paste(data_dir, "/All_Results_Table.tsv", sep='')
write.table(all_data, file=csv_file, sep='\t', row.names = FALSE)

my_comparisons <- list(c("0", "1"), c("1", "2"), c("0", "2"))

control_grp <- all_data %>%
    filter(Group == "Control")

control_grp_snp1 <- control_grp %>%
    filter(Position == 113735574)

control_grp_snp2 <- control_grp %>%
    filter(Position == 113797453)   

control_grp_snp3 <- control_grp %>%
    filter(Position == 113217527)    

control_grp_snp1_wo749 <- control_grp_snp1 %>%
    filter(grepl('RZ749', Condition) == FALSE)

control_grp_snp2_wo749 <- control_grp_snp2 %>%
    filter(grepl('RZ749', Condition) == FALSE)  

control_grp_snp3_wo749 <- control_grp_snp3 %>%
    filter(grepl('RZ749', Condition) == FALSE)    

## baseline IL1B expression
pdf_title <- sprintf("%s/control_exp_IL1B.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(control_grp_snp1, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(control_grp_snp2, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(control_grp_snp3, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()

## baseline IL1B expression without RZ749
pdf_title <- sprintf("%s/control_exp_IL1B_without_RZ749.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(control_grp_snp1_wo749, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(control_grp_snp2_wo749, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(control_grp_snp3_wo749, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()

## baseline TNFa expression
pdf_title <- sprintf("%s/control_exp_TNFa.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(control_grp_snp1, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(control_grp_snp2, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(control_grp_snp3, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()

## baseline TNFa expression without RZ749
pdf_title <- sprintf("%s/control_exp_TNFa_without_RZ749.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(control_grp_snp1_wo749, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(control_grp_snp2_wo749, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(control_grp_snp3_wo749, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()



## treated group
treated_grp <- all_data %>%
    filter(Group == "Treated (LPS)")

treated_grp_snp1 <- treated_grp %>%
    filter(Position == 113735574)

treated_grp_snp2 <- treated_grp %>%
    filter(Position == 113797453)   

treated_grp_snp3 <- treated_grp %>%
    filter(Position == 113217527)   

treated_grp_snp1_wo749 <- treated_grp_snp1 %>%
    filter(grepl('RZ749', Condition) == FALSE)

treated_grp_snp2_wo749 <- treated_grp_snp2 %>%
    filter(grepl('RZ749', Condition) == FALSE)  

treated_grp_snp3_wo749 <- treated_grp_snp3 %>%
    filter(grepl('RZ749', Condition) == FALSE)  

## treated IL1B expression
pdf_title <- sprintf("%s/treated_exp_IL1B.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(treated_grp_snp1, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(treated_grp_snp2, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(treated_grp_snp3, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()

## treated IL1b expression without RZ749
pdf_title <- sprintf("%s/treated_exp_IL1B_without_RZ749.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(treated_grp_snp1_wo749, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(treated_grp_snp2_wo749, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(treated_grp_snp3_wo749, x="Allele", y="MFI IL1b Monocytes", color="Allele", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()

## treated TNFa expression
pdf_title <- sprintf("%s/treated_exp_TNFa.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(treated_grp_snp1, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(treated_grp_snp2, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(treated_grp_snp3, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()

## treated TNFa expression without RZ749
pdf_title <- sprintf("%s/treated_exp_TNFa_without_RZ749.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(treated_grp_snp1_wo749, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(treated_grp_snp2_wo749, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(treated_grp_snp3_wo749, x="Allele", y="MFI TNFa Monocytes", color="Allele", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()



## control vs treated

all_data_snp1 <- all_data %>%
    filter(Position == 113735574)

all_data_snp2 <- all_data %>%
    filter(Position == 113797453)   

all_data_snp3 <- all_data %>%
    filter(Position == 113217527) 

pdf_title <- sprintf("%s/cntrl_vs_treated_IL1B.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(all_data_snp1, x="Allele", y="MFI IL1b Monocytes", color="Group", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p2 <- ggboxplot(all_data_snp2, x="Allele", y="MFI IL1b Monocytes", color="Group", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p3 <- ggboxplot(all_data_snp3, x="Allele", y="MFI IL1b Monocytes", color="Group", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p1 + p2 + p3
dev.off()

pdf_title <- sprintf("%s/cntrl_vs_treated_TNFa.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(all_data_snp1, x="Allele", y="MFI TNFa Monocytes", color="Group", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p2 <- ggboxplot(all_data_snp2, x="Allele", y="MFI TNFa Monocytes", color="Group", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p3 <- ggboxplot(all_data_snp3, x="Allele", y="MFI TNFa Monocytes", color="Group", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p1 + p2 + p3
dev.off()

## control vs treated without RZ749
all_data_snp1_wo749 <- all_data_snp1 %>%
    filter(grepl('RZ749', Condition) == FALSE)

all_data_snp2_wo749 <- all_data_snp2 %>%
    filter(grepl('RZ749', Condition) == FALSE)  

all_data_snp3_wo749 <- all_data_snp3 %>%
    filter(grepl('RZ749', Condition) == FALSE) 

pdf_title <- sprintf("%s/cntrl_vs_treated_IL1B_without_RZ749.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(all_data_snp1_wo749, x="Allele", y="MFI IL1b Monocytes", color="Group", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p2 <- ggboxplot(all_data_snp2_wo749, x="Allele", y="MFI IL1b Monocytes", color="Group", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p3 <- ggboxplot(all_data_snp3_wo749, x="Allele", y="MFI IL1b Monocytes", color="Group", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p1 + p2 + p3
dev.off()

pdf_title <- sprintf("%s/cntrl_vs_treated_TNFa_without_RZ749.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(all_data_snp1_wo749, x="Allele", y="MFI TNFa Monocytes", color="Group", add="jitter") + ggtitle("snp 113735574 A>G") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p2 <- ggboxplot(all_data_snp2_wo749, x="Allele", y="MFI TNFa Monocytes", color="Group", add="jitter") + ggtitle("snp 113797453 A>T") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p3 <- ggboxplot(all_data_snp3_wo749, x="Allele", y="MFI TNFa Monocytes", color="Group", add="jitter") + ggtitle("snp 113217527 T>G") + stat_compare_means(aes(group=Group), method = "t.test", label = "p.format")
p1 + p2 + p3
dev.off()


## control vs treated percentage
# all_data_snp1_percent_change <- all_data_snp1 %>%
#     group_by(PatientID) %>%
#     mutate(Gene = ifelse(Group == 'Control', 'TFNa', 'IL1B')) %>%
#     mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI IL1b Monocytes`/`MFI IL1b Monocytes`[Group == 'Control']) * 100, (`MFI TNFa Monocytes`[Group == 'Treated (LPS)']/`MFI TNFa Monocytes`) * 100))

# pdf_title <- sprintf("%s/1-cntrl_vs_treated_percent_change_overall.pdf", plots_dir)
# pdf(pdf_title, width = 5, height = 5)
# ggboxplot(all_data_snp1_percent_change, x="Gene", y="Percent change", color="Gene", add = "jitter")
# dev.off()

## control vs treated TNFa percentage change
all_data_snp1_TNFa_percent_change <- all_data_snp1 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI TNFa Monocytes`/`MFI TNFa Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

all_data_snp2_TNFa_percent_change <- all_data_snp2 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI TNFa Monocytes`/`MFI TNFa Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

all_data_snp3_TNFa_percent_change <- all_data_snp3 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI TNFa Monocytes`/`MFI TNFa Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

pdf_title <- sprintf("%s/cntrl_vs_treated_TNFa_percent_change.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(all_data_snp1_TNFa_percent_change, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113735574 A>G TNFa Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(all_data_snp2_TNFa_percent_change, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113797453 A>T TNFa Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(all_data_snp3_TNFa_percent_change, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113217527 T>G TNFa Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()

## control vs treated TNFa percentage change without RZ749
all_data_snp1_TNFa_percent_change_wo749 <- all_data_snp1_wo749 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI TNFa Monocytes`/`MFI TNFa Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

all_data_snp2_TNFa_percent_change_wo749 <- all_data_snp2_wo749 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI TNFa Monocytes`/`MFI TNFa Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

all_data_snp3_TNFa_percent_change_wo749 <- all_data_snp3_wo749 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI TNFa Monocytes`/`MFI TNFa Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

pdf_title <- sprintf("%s/cntrl_vs_treated_TNFa_percent_change_without_RZ749.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(all_data_snp1_TNFa_percent_change_wo749, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113735574 A>G TNFa Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(all_data_snp2_TNFa_percent_change_wo749, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113797453 A>T TNFa Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(all_data_snp3_TNFa_percent_change_wo749, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113217527 T>G TNFa Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()

## control vs treated IL1B percentage change
all_data_snp1_IL1b_percent_change <- all_data_snp1 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI IL1b Monocytes`/`MFI IL1b Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

all_data_snp2_IL1b_percent_change <- all_data_snp2 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI IL1b Monocytes`/`MFI IL1b Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

all_data_snp3_IL1b_percent_change <- all_data_snp3 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI IL1b Monocytes`/`MFI IL1b Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

pdf_title <- sprintf("%s/cntrl_vs_treated_IL1b_percent_change.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(all_data_snp1_IL1b_percent_change, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113735574 A>G IL1b Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(all_data_snp2_IL1b_percent_change, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113797453 A>T IL1b Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(all_data_snp3_IL1b_percent_change, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113217527 T>G IL1b Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()


## control vs treated IL1B percentage change without RZ749
all_data_snp1_IL1b_percent_change_wo749 <- all_data_snp1_wo749 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI IL1b Monocytes`/`MFI IL1b Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

all_data_snp2_IL1b_percent_change_wo749 <- all_data_snp2_wo749 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI IL1b Monocytes`/`MFI IL1b Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

all_data_snp3_IL1b_percent_change_wo749 <- all_data_snp3_wo749 %>%
    group_by(PatientID) %>%
    mutate(`Percent change` = ifelse(Group == "Treated (LPS)", (`MFI IL1b Monocytes`/`MFI IL1b Monocytes`[Group == 'Control']) * 100, NA)) %>%
    filter(Group == 'Treated (LPS)')

pdf_title <- sprintf("%s/cntrl_vs_treated_IL1b_percent_change_without_RZ749.pdf", plots_dir)
pdf(pdf_title, width = 15, height = 5)
p1 <- ggboxplot(all_data_snp1_IL1b_percent_change_wo749, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113735574 A>G IL1b Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p2 <- ggboxplot(all_data_snp2_IL1b_percent_change_wo749, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113797453 A>T IL1b Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p3 <- ggboxplot(all_data_snp3_IL1b_percent_change_wo749, x="Allele", y="Percent change", color="Allele", add = "jitter") + ggtitle("snp 113217527 T>G IL1b Percent change") + stat_compare_means(comparisons=my_comparisons, method="t.test")
p1 + p2 + p3
dev.off()


#Try log2 and FDR rather than percentage change and t-testing.

## ==========================================
## Log2 fold-change SNP analysis (FDR only)
## ==========================================
library(tidyverse)

# --- 1. DATA PREPARATION ---
analysis_df <- all_data %>%
  select(PatientID, Position, Allele, Group, 
         `MFI IL1b Monocytes`, `MFI TNFa Monocytes`) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(`MFI IL1b Monocytes`, `MFI TNFa Monocytes`)
  ) %>%
  filter(complete.cases(.)) %>%
  # Safe renaming by column index to handle specific naming conventions
  rename(
    IL1b_Ctrl = 4, IL1b_LPS = 5, 
    TNFa_Ctrl = 6, TNFa_LPS = 7
  ) %>%
  mutate(
    IL1b_Response = log2(IL1b_LPS / IL1b_Ctrl),
    TNFa_Response = log2(TNFa_LPS / TNFa_Ctrl),
    Allele = as.factor(Allele)
  )

snps <- unique(analysis_df$Position)
cytokines <- c("IL1b_Response", "TNFa_Response")

# --- 2. STATISTICAL CALCULATION ---
results_master <- data.frame()

for (cyto in cytokines) {
  for (snp in snps) {
    df_sub <- analysis_df %>% filter(Position == snp)
    
    # Linear Model (ANOVA logic)
    if(length(unique(df_sub$Allele)) > 1) {
      fit <- lm(as.formula(paste(cyto, "~ Allele")), data = df_sub)
      f <- summary(fit)$fstatistic
      p_raw <- pf(f[1], f[2], f[3], lower.tail = FALSE)
    } else {
      p_raw <- NA 
    }
    
    results_master <- rbind(results_master, data.frame(
      SNP = snp, 
      Cytokine = cyto, 
      P_Raw = p_raw
    ))
  }
}

# Adjust p-values for the 3 SNPs within each cytokine group
results_master <- results_master %>%
  group_by(Cytokine) %>%
  mutate(P_Adj = p.adjust(P_Raw, method = "BH")) %>%
  ungroup()

# --- 3. PLOTTING TO PDF ---
for (i in 1:nrow(results_master)) {
  current_res <- results_master[i, ]
  
  # Skip if no data
  if(is.na(current_res$P_Adj)) next
  
  df_plot <- analysis_df %>% filter(Position == current_res$SNP)
  
  p <- ggplot(df_plot, aes(x = Allele, y = .data[[current_res$Cytokine]], fill = Allele)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    # Set y-axis limits to -1 to 6
    coord_cartesian(ylim = c(-1, 6)) + 
    theme_classic() +
    labs(
      title = paste("SNP:", current_res$SNP),
      subtitle = paste0("FDR Adjusted p = ", signif(current_res$P_Adj, 3)),
      y = paste(current_res$Cytokine, "(Log2 Fold Change)"),
      x = "Genotype (0/1/2)"
    ) +
    theme(legend.position = "none")

  # Save as PDF
  file_name <- paste0("Final_Plot_", current_res$SNP, "_", current_res$Cytokine, ".pdf")
  ggsave(file.path(plots_dir, file_name), plot = p, device = "pdf", width = 5, height = 5)
}

# Save stats to CSV
write.csv(results_master, file.path(plots_dir, "SNP_Statistical_Results.csv"), row.names = FALSE)

#Methods:
#Statistical Analysis
#To evaluate the association between germline variants and the functional response of immune cell subsets, median fluorescence intensity (MFI) values for IL-1b and TNFa were acquired from monocyte populations using flow cytometry. For each patient, the magnitude of the LPS-induced inflammatory response was quantified by calculating the log2 fold change (LFC) between treated and baseline (control) MFIs. In this scale, a value of 1 represents a doubling of cytokine production relative to control, while a value of 0 indicates no change.
#Individual SNP genotypes were analyzed as categorical predictors (0: Ref/Ref; 1: Ref/Alt; 2: Alt/Alt). The relationship between genotype and the cytokine response was assessed via linear regression models. This methodology identifies variants that significantly modulate the magnitude of cytokine induction (Gene-Environment interaction). To account for multiple comparisons across the three candidate loci, p-values were adjusted using the Benjamini-Hochberg False Discovery Rate (FDR) method.
#All statistical procedures were conducted using R software (version 4.x). Data visualization was performed using the ggplot2 package, with results presented as boxplots showing individual patient data points. Significance thresholds were set at an FDR-adjusted p < 0.05.


#Power calculations
#Summary for your team: "Based on the current pilot data of 91 patients, the best-performing variant (SNP 113217527) would require a total cohort of approximately 433 patients to reach a statistical significance of $p < 0.05$ with 80% power. The other tested variants show negligible effect sizes and are unlikely to reach significance without a significantly larger multi-center cohort."

#I can calculate this:

# You may need to install this package first: install.packages("pwr")
library(pwr)

# Function to calculate required N for a specific model
calculate_required_n <- function(model) {
  # Get the Cohen's f2 effect size from R-squared
  r_sq <- summary(model)$r.squared
  f2 <- r_sq / (1 - r_sq)
  
  if (f2 < 0.001) return(NA) # Effect is too small to calculate reasonably
  
  # Calculate required denominator degrees of freedom (v) 
  # for 80% power at alpha 0.05 with 2 degrees of freedom (3 genotypes)
  pw <- pwr.f2.test(u = 2, f2 = f2, sig.level = 0.05, power = 0.80)
  
  # Total N = v + u + 1
  total_n <- ceiling(pw$v + pw$u + 1)
  return(total_n)
}

# --- RUNNING THE PREDICTION ---
power_results <- data.frame()

for (cyto in cytokines) {
  for (snp in snps) {
    df_sub <- analysis_df %>% filter(Position == snp)
    
    if(length(unique(df_sub$Allele)) > 1) {
      fit <- lm(as.formula(paste(cyto, "~ Allele")), data = df_sub)
      
      n_needed <- calculate_required_n(fit)
      
      power_results <- rbind(power_results, data.frame(
        SNP = snp,
        Cytokine = cyto,
        Current_P = summary(fit)$coefficients[1,4], # Just for reference
        Predicted_N_Needed = n_needed
      ))
    }
  }
}

print(power_results)

##################################################
###Venessa Code
#Repeat plots for the baseline comparisons only.  Correct for multiple testing

library(tidyverse)

# --- 1. COMBINE CONTROL TABLES ---
baseline_df <- bind_rows(
  control_grp_snp1_wo749,
  control_grp_snp2_wo749,
  control_grp_snp3_wo749
) %>%
  select(
    PatientID,
    Position,
    Ref,
    Alt,
    Allele,
    `MFI IL1b Monocytes`,
    `MFI TNFa Monocytes`
  ) %>%
  filter(complete.cases(.)) %>%
  rename(
    IL1b_Baseline = `MFI IL1b Monocytes`,
    TNFa_Baseline = `MFI TNFa Monocytes`
  ) %>%
  # Keep numeric order consistent for 0,1,2
  mutate(Allele = factor(Allele, levels = c(0, 1, 2)))

snps <- unique(baseline_df$Position)
cytokines <- c("IL1b_Baseline", "TNFa_Baseline")


# --- 2. STATISTICAL CALCULATION ---
results_baseline <- data.frame()

for (cyto in cytokines) {
  for (snp in snps) {
    
    df_sub <- baseline_df %>% filter(Position == snp)
    
    if (length(unique(df_sub$Allele)) > 1) {
      fit <- lm(as.formula(paste(cyto, "~ Allele")), data = df_sub)
      f <- summary(fit)$fstatistic
      p_raw <- pf(f[1], f[2], f[3], lower.tail = FALSE)
    } else {
      p_raw <- NA
    }
    
    results_baseline <- rbind(
      results_baseline,
      data.frame(
        SNP = snp,
        Cytokine = cyto,
        P_Raw = p_raw
      )
    )
  }
}

results_baseline <- results_baseline %>%
  group_by(Cytokine) %>%
  mutate(P_Adj = p.adjust(P_Raw, method = "BH")) %>%
  ungroup()

# --- 3. PLOTTING WITH CAPITALISED GENOTYPE LABELS ---
# Mapping numeric Allele -> display label
allele_labels <- c("0" = "Ref/Ref", "1" = "Het", "2" = "Alt/Alt")

for (i in 1:nrow(results_baseline)) {
  
  current_res <- results_baseline[i, ]
  if (is.na(current_res$P_Adj)) next
  
  df_plot <- baseline_df %>% filter(Position == current_res$SNP)
  
  p <- ggplot(
    df_plot,
    aes(x = Allele, y = .data[[current_res$Cytokine]], fill = Allele)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    scale_x_discrete(labels = allele_labels) +           # <- capitalised labels
    theme_classic() +
    labs(
      title = paste("SNP:", current_res$SNP),
      subtitle = paste0("FDR Adjusted p = ", signif(current_res$P_Adj, 3)),
      y = paste(current_res$Cytokine, "(MFI, baseline)"),
      x = "Genotype (0 = Ref/Ref, 1 = Het, 2 = Alt/Alt)"
    ) +
    theme(legend.position = "none")
  
  file_name <- paste0(
    "Baseline_Plot_",
    current_res$SNP, "_",
    current_res$Cytokine, ".pdf"
  )
  
  ggsave(
    file.path(plots_dir, file_name),
    plot = p,
    device = "pdf",
    width = 5,
    height = 5
  )
}

# --- 4. SAVE STATISTICS ---
write.csv(
  results_baseline,
  file.path(plots_dir, "Baseline_SNP_Statistical_Results.csv"),
  row.names = FALSE
)
