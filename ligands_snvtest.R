#setwd("~/project/IntellanceII/Full_Intellance")

# ---- load packages ----
library(dplyr)
library(tidyverse)
library(ggplot2)
library(MASS)

# ---- load data ---- 
source("Scripts/job_gg_theme.R")
source('Scripts/load_full_intellance.R')

#Load VST transformed count data
counts <- read.csv("norm_vst_counts_full_intellance.csv",row.names = c(1))
colnames(counts) <- gsub('\\X', '', colnames(counts))
colnames(counts) <- gsub('\\.', '-', colnames(counts))
ligands <- c("EGFR|chr7:55M","AREG|chr4:74M","EREG|chr4:74M","EPGN|chr4:74M","TGFA|chr2:71M","EGF|chr4:110M","HBEGF|chr5:140M","BTC|chr4:75M","CTGF|chr6:132M")
counts_ligands <- as.data.frame(t(counts[rownames(counts)%in%ligands,])) 
colnames(counts_ligands) <- c("AREG","HBEGF","CTGF","EREG","EGF","EGFR","TGFA","BTC")

#Load vIII data
vIII <- read.delim(file="../../../mnt/neuro-genomic-2/intellanceII/Almac/tables/v3_reads.txt")
vIII$sample <- gsub("/Aligned.*","",vIII$sample)
vIII$vIII_status <- ifelse(vIII$vIII/(vIII$wt+vIII$vIII) > 0.01, "mutated", "normal")
vIII$reads <- vIII$vIII + vIII$wt
vIII[vIII$reads<15,c(4)] <- "normal"

#Merge count and vIII data
counts_ligands$sample <- rownames(counts_ligands)
counts_ligands <- merge(vIII[,c(1,4)],counts_ligands,by="sample")

#Convert to Intellance ID
source("sample_id_conversion.R")
intellance_id <- id_conversion(counts_ligands$sample)
colnames(intellance_id)[1] <- "sample"
counts_ligands <- merge(intellance_id,counts_ligands,by="sample")
rownames(counts_ligands) <- counts_ligands$IntellanceID
counts_ligands$sample <- NULL
colnames(counts_ligands)[1] <- "sample"

#SNVs
A289 <- readRDS(file="SNV_A289_full.RData")
A289$sample <- sub(pattern = "\\.\\d+$", replacement = "",A289$sample) 

R108 <- readRDS(file="SNV_R108_full.RData")
R108$sample <- sub(pattern = "\\.\\d+$", replacement = "",R108$sample) 

G598 <- readRDS(file="SNV_G598_full.RData")
G598$sample <- sub(pattern = "\\.\\d+$", replacement = "",G598$sample) 

#Merge SNVs with count data
snv <- as.data.frame(counts_ligands$sample) 
colnames(snv) <- "sample"
snv$mutation <- "none"

for (i in 1:nrow(snv)) {
  if (snv$sample[i] %in% intersect(A289$sample,R108$sample) == TRUE) {
    snv$mutation[i] <- "A289V/T/+R108K/G" }
  else if (snv$sample[i] %in% intersect(A289$sample,G598$sample) == TRUE){
    snv$mutation[i] <- "A289V/T/D+G598V/A"
  }
  else if (snv$sample[i] %in% intersect(R108$sample,G598$sample) == TRUE){
    snv$mutation[i] <- "R108K/G+G598V/A"
  }  
  else if (snv$sample[i] %in% R108$sample){
    snv$mutation[i] <- "R108K/G"
  }   
  else if (snv$sample[i] %in% G598$sample){
    snv$mutation[i] <- "G598V/A"
  }
  else if (snv$sample[i] %in% A289$sample){
    snv$mutation[i] <- "A289V/T/D"
  } 
}

counts_ligands <- merge(snv,counts_ligands,by="sample")
counts_ligands$amplified <- ifelse(counts_ligands$EGFR>10.71,"amplified","normal")

# ---- significance tests ---
ampli_table <- readRDS("output/CNV_amplified_table.Rds")
ampli_table <- dplyr::left_join(ampli_table,intellance_id,by=c("X"="sample"))

EGFR_ampli <- ampli_table %>% dplyr::filter(EGFR_CNV > 5) %>% dplyr::pull(IntellanceID)

#Within EGFR amplified
counts_ligands_amplified_dna <- counts_ligands %>% 
                                dplyr::filter(sample %in% EGFR_ampli) %>%
                                dplyr::select(-amplified)

counts_ligands_amplified_dna$selected_mut <- ifelse(counts_ligands_amplified_dna$mutation == "none",0,1)
counts_ligands_amplified_dna$A289 <- grepl("A289",counts_ligands_amplified_dna$mutation)
counts_ligands_amplified_dna$G598 <- grepl("G598",counts_ligands_amplified_dna$mutation)
counts_ligands_amplified_dna$R108 <- grepl("R108",counts_ligands_amplified_dna$mutation)

# counts_lig_ampli_vIIIneg <- counts_ligands_amplified_dna %>% dplyr::filter(vIII_status == "normal")
# counts_lig_ampli_vIIIpos <- counts_ligands_amplified_dna %>% dplyr::filter(vIII_status == "mutated")
# sum(counts_ligands_amplified_dna$mutation == "none")
# sum(counts_ligands_amplified_dna$A289)
# sum(counts_ligands_amplified_dna$G598)
# sum(counts_ligands_amplified_dna$R108)
# sum(counts_ligands_amplified_dna$mutation == "none")
# sum(counts_ligands_amplified_dna$A289)
# sum(counts_ligands_amplified_dna$G598)
# sum(counts_ligands_amplified_dna$R108)

fisher.test(counts_ligands_amplified_dna$selected_mut,counts_ligands_amplified$vIII_status) #p=0.095
fisher.test(counts_ligands_amplified_dna$G598,counts_ligands_amplified$vIII_status) #p=0.0.056

counts_ligands_amplified <- counts_ligands %>%
                            dplyr::filter(amplified == "amplified")
counts_ligands_amplified$selected_mut <- ifelse(counts_ligands_amplified$mutation == "none","wildtype","mutated")

wilcox.test(BTC ~ vIII_status, data=counts_ligands_amplified)
wilcox.test(HBEGF ~ vIII_status, data=counts_ligands_amplified)
wilcox.test(TGFA ~ vIII_status, data=counts_ligands_amplified)
wilcox.test(EGF ~ vIII_status, data=counts_ligands_amplified)
wilcox.test(EREG ~ vIII_status, data=counts_ligands_amplified)
wilcox.test(AREG ~ vIII_status, data=counts_ligands_amplified)

wilcox.test(BTC ~ selected_mut, data=counts_ligands_amplified)
wilcox.test(HBEGF ~ selected_mut, data=counts_ligands_amplified)
wilcox.test(TGFA ~ selected_mut, data=counts_ligands_amplified)
wilcox.test(EGF ~ selected_mut, data=counts_ligands_amplified)
wilcox.test(EREG ~ selected_mut, data=counts_ligands_amplified)
wilcox.test(AREG ~ selected_mut, data=counts_ligands_amplified)

#Comparison non-ampli and ampli
counts_lig <- counts_ligands 
counts_lig$selected_mut <- ifelse(counts_lig$mutation == "none",0,1)

wilcox.test(BTC ~ amplified, data=counts_lig) 
wilcox.test(HBEGF ~ amplified, data=counts_lig) 
wilcox.test(TGFA ~ amplified, data=counts_lig) 
wilcox.test(EGF ~ amplified, data=counts_lig) 
wilcox.test(EREG ~ amplified, data=counts_lig) 
wilcox.test(AREG ~ amplified, data=counts_lig) 

fisher.test(counts_lig$selected_mut,counts_lig$ampli)

#Amplified but no SNV or vIII
tmp <- counts_lig %>%
  dplyr::filter(amplified=="amplified") %>%
  dplyr::filter(mutation != "none" | vIII_status != "normal") %>%
  dplyr::pull(sample)

no_mutation_egfr <- counts_lig %>%
                    dplyr::filter(!sample%in%tmp)

wilcox.test(BTC ~ amplified, data=no_mutation_egfr) 
wilcox.test(HBEGF ~ amplified, data=no_mutation_egfr) 
wilcox.test(TGFA ~ amplified, data=no_mutation_egfr) 
wilcox.test(EGF ~ amplified, data=no_mutation_egfr) 
wilcox.test(EREG ~ amplified, data=no_mutation_egfr) 
wilcox.test(AREG ~ amplified, data=no_mutation_egfr) 

# ---- DESeq2 test EGFR SNV ----
mm <- dplyr::left_join(counts_ligands[,c("sample","mutation","amplified")],intellance_id,by=c("sample"="IntellanceID")) %>%
                       dplyr::filter(amplified=="amplified")
mm$selected_mut <- ifelse(mm$mutation == "none","wildtype","mutated")

e <- raw_counts[,colnames(raw_counts)%in%mm$sample.y] %>%
     dplyr::filter(rowMeans(.)>3) %>%
     dplyr::filter(!rownames(.)%in%rownames(.)[grepl("chrM",rownames(.))])

cond <- factor(mm[match(colnames(e),mm$sample.y),]$selected_mut, levels=c('wildtype','mutated') ) 

colnames(e) == mm[match(colnames(e),mm$sample.y),]$sample.y

dds <- DESeq2::DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond) %>%
       DESeq2::estimateSizeFactors(.) %>%
       DESeq2::DESeq(parallel = T) %>%
       DESeq2::results(.)

DESeq2::plotMA(dds, alpha = 0.01, main='DESeq2: Resection 1 vs 2',ylim=c(-3,3))

resOrdered <- dds[order(dds$padj),]
rownames(resOrdered) <- gsub("\\|chr.*","",rownames(resOrdered))
rownames(resOrdered) <- gsub(".*\\|","",rownames(resOrdered))

rownames(resOrdered)[which(rownames(resOrdered)=="AL512306.2")] <- "MDM4"
rownames(resOrdered)[which(rownames(resOrdered)=="AL449423.1 ")] <- "CDKN2B"

dtable <- as.data.frame(resOrdered) %>% dplyr::filter(abs(log2FoldChange) > 0.5 & padj < 0.05)
#write.csv(dtable,"output/DESeq2_results_missense_cond_EGFRampli.csv")

labels <- rownames(dtable)[!rownames(dtable)%in%c("B4GALNT1","BCL11B","C1GALT1P1","METTL1","BCL11A","LINC01305","COL11A1","TSFM")]

EnhancedVolcano(title = NULL,
                subtitle = NULL,
                resOrdered,
                lab = rownames(resOrdered),
                axisLabSize = 12,
                legendPosition = "bottom",
                legendLabSize = 9,
                x = 'log2FoldChange',
                y = 'padj',  
                pCutoff = 0.05,
                FCcutoff = 0.5,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                col = c("black", "forestgreen", "deepskyblue3", "firebrick1"),
                xlim = c(-2,2), 
                selectLab = labels,
                ylim = c(0,7),
                captionLabSize = 9)
#ggsave("Figures/EnhancedVolcano_Missense_Condition_EGFRampli.pdf",width = 8, height = 5.7)

# ---- Ligand plots ----
colors <- c("none"="lightgrey", "R108K/G"="forestgreen", "G598V/A"="red","A289V/T/D"="blue","vIII mutated"="orange","EGFR amplified"="violet")

EGFR <- ggplot(counts_ligands, aes(x=reorder(sample, EGFR), y = EGFR)) +
  geom_point(pch=20) +
  geom_point(data = subset(counts_ligands,mutation=="none"&vIII_status=="normal"), 
             mapping = aes(x = reorder(sample, EGFR), y = 7, col="none"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="R108K/G"), 
             mapping = aes(x = reorder(sample, EGFR), y = 6.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="G598V/A"), 
             mapping = aes(x = reorder(sample, EGFR), y = 6.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D"), 
             mapping = aes(x = reorder(sample, EGFR), y = 6.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x = reorder(sample, EGFR), y = 6.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x =reorder(sample, EGFR), y = 6.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, EGFR), y = 6.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, EGFR), y = 6.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,vIII_status=="mutated"), 
             mapping = aes(x = reorder(sample, EGFR), y = 6.20, col="vIII mutated"),size=1.2, pch=15, alpha=1) +  
  geom_point(data = subset(counts_ligands,vIII_status=="mutated"), 
             mapping = aes(x = reorder(sample, EGFR), y = 6.20, col="vIII mutated"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,EGFR>10.71), 
             mapping = aes(x = reorder(sample, EGFR), y = 6, col="EGFR amplified"),size=1.2, pch=15, alpha=1) +  
  scale_color_manual(values = colors) +
  labs(y="EGFR expression (log2)", colour="EGFR mutation") +
  facet_grid(~amplified , drop=F, scales="free_x", space = "free_x") +
  job_gg_theme
#ggsave("EGFR_sorted_mutations.png")
#ggsave("EGFR_sorted_mutations.pdf")

BTC <- ggplot(counts_ligands, aes(x=reorder(sample, BTC), y = BTC)) +
  geom_point(pch=20) +
  geom_point(data = subset(counts_ligands,mutation=="none"&vIII_status=="normal"), 
             mapping = aes(x = reorder(sample, BTC), y = 2, col="none"),size=1, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="R108K/G"), 
             mapping = aes(x = reorder(sample, BTC), y = 1.90, col="R108K/G"),size=1, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="G598V/A"), 
             mapping = aes(x = reorder(sample, BTC), y = 1.80, col="G598V/A"),size=1, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D"), 
             mapping = aes(x = reorder(sample, BTC), y = 1.70, col="A289V/T/D"),size=1, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x = reorder(sample, BTC), y = 1.70, col="A289V/T/D"),size=1, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x =reorder(sample, BTC), y = 1.90, col="R108K/G"),size=1, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, BTC), y = 1.70, col="A289V/T/D"),size=1, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, BTC), y = 1.80, col="G598V/A"),size=1, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,vIII_status=="mutated"), 
             mapping = aes(x = reorder(sample, BTC), y = 1.60, col="vIII mutated"),size=1, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,EGFR>10.71), 
             mapping = aes(x = reorder(sample, BTC), y = 1.50, col="EGFR amplified"),size=1.2, pch=15, alpha=1) + 
  scale_color_manual(values = colors) +
  labs(y="BTC expression (log2)", colour="EGFR mutation") +
  facet_grid(~amplified , drop=F, scales="free_x", space = "free_x") +
  job_gg_theme
#ggsave("BTC_sorted_mutations.png")
#ggsave("BTC_sorted_mutations.pdf")

TGFA <- ggplot(counts_ligands, aes(x=reorder(sample, TGFA), y = TGFA)) +
  geom_point(pch=20) +
  geom_point(data = subset(counts_ligands,mutation=="none"&vIII_status=="normal"), 
             mapping = aes(x = reorder(sample, TGFA), y = 2, col="none"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="R108K/G"), 
             mapping = aes(x = reorder(sample, TGFA), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="G598V/A"), 
             mapping = aes(x = reorder(sample, TGFA), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D"), 
             mapping = aes(x = reorder(sample, TGFA), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x = reorder(sample, TGFA), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x =reorder(sample, TGFA), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, TGFA), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, TGFA), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,vIII_status=="mutated"), 
             mapping = aes(x = reorder(sample, TGFA), y = 1.20, col="vIII mutated"),size=1.2, pch=15, alpha=1) + 
  geom_point(data = subset(counts_ligands,EGFR>10.71), 
             mapping = aes(x = reorder(sample, TGFA), y = 1, col="EGFR amplified"),size=1.2, pch=15, alpha=1) + 
  scale_color_manual(values = colors) +
  labs(y="TGFA expression (log2)", colour="EGFR mutation") +
  facet_grid(~amplified , drop=F, scales="free_x", space = "free_x") +
  job_gg_theme
#ggsave("TGFA_sorted_mutations.png")
#ggsave("TGFA_sorted_mutations.pdf")

EGF <- ggplot(counts_ligands, aes(x=reorder(sample, EGF), y = EGF)) +
  geom_point(pch=20) +
  geom_point(data = subset(counts_ligands,mutation=="none"&vIII_status=="normal"), 
             mapping = aes(x = reorder(sample, EGF), y = 2, col="none"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="R108K/G"), 
             mapping = aes(x = reorder(sample, EGF), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="G598V/A"), 
             mapping = aes(x = reorder(sample, EGF), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D"), 
             mapping = aes(x = reorder(sample, EGF), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x = reorder(sample, EGF), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x =reorder(sample, EGF), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, EGF), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, EGF), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,vIII_status=="mutated"), 
             mapping = aes(x = reorder(sample, EGF), y = 1.20, col="vIII mutated"),size=1.2, pch=15, alpha=1) +  
  geom_point(data = subset(counts_ligands,EGFR>10.71), 
             mapping = aes(x = reorder(sample, EGF), y = 1, col="EGFR amplified"),size=1.2, pch=15, alpha=1) + 
  scale_color_manual(values = colors) +
  labs(y="EGF expression (log2)", colour="EGFR mutation") +
  facet_grid(~amplified , drop=F, scales="free_x", space = "free_x") +
  job_gg_theme
#ggsave("EGF_sorted_mutations.png")
#ggsave("EGF_sorted_mutations.pdf")

EREG <- ggplot(counts_ligands, aes(x=reorder(sample, EREG), y = EREG)) +
  geom_point(pch=20) +
  geom_point(data = subset(counts_ligands,mutation=="none"&vIII_status=="normal"), 
             mapping = aes(x = reorder(sample, EREG), y = 2, col="none"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="R108K/G"), 
             mapping = aes(x = reorder(sample, EREG), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="G598V/A"), 
             mapping = aes(x = reorder(sample, EREG), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D"), 
             mapping = aes(x = reorder(sample, EREG), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x = reorder(sample, EREG), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x =reorder(sample, EREG), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, EREG), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, EREG), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,vIII_status=="mutated"), 
             mapping = aes(x = reorder(sample, EREG), y = 1.20, col="vIII mutated"),size=1.2, pch=15, alpha=1) +  
  geom_point(data = subset(counts_ligands,EGFR>10.71), 
             mapping = aes(x = reorder(sample, EREG), y = 1, col="EGFR amplified"),size=1.2, pch=15, alpha=1) + 
  scale_color_manual(values = colors) +
  labs(y="EREG expression (log2)", colour="EGFR mutation") +
  facet_grid(~amplified , drop=F, scales="free_x", space = "free_x") +
  job_gg_theme
#ggsave("EREG_sorted_mutations.png")
#ggsave("EREG_sorted_mutations.pdf")

CTGF <- ggplot(counts_ligands, aes(x=reorder(sample, CTGF), y = CTGF)) +
  geom_point(pch=20) +
  geom_point(data = subset(counts_ligands,mutation=="none"&vIII_status=="normal"), 
             mapping = aes(x = reorder(sample, CTGF), y = 2, col="none"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="R108K/G"), 
             mapping = aes(x = reorder(sample, CTGF), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="G598V/A"), 
             mapping = aes(x = reorder(sample, CTGF), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D"), 
             mapping = aes(x = reorder(sample, CTGF), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x = reorder(sample, CTGF), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x =reorder(sample, CTGF), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, CTGF), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, CTGF), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,vIII_status=="mutated"), 
             mapping = aes(x = reorder(sample, CTGF), y = 1.20, col="vIII mutated"),size=1.2, pch=15, alpha=1) +  
  geom_point(data = subset(counts_ligands,EGFR>10.71), 
             mapping = aes(x = reorder(sample, CTGF), y = 1, col="EGFR amplified"),size=1.2, pch=15, alpha=1) + 
  scale_color_manual(values = colors) +
  labs(y="CTGF expression (log2)", colour="EGFR mutation") +
  facet_grid(~amplified , drop=F, scales="free_x", space = "free_x") +
  job_gg_theme


HBEGF <- ggplot(counts_ligands, aes(x=reorder(sample, HBEGF), y = HBEGF)) +
  geom_point(pch=20) +
  geom_point(data = subset(counts_ligands,mutation=="none"&vIII_status=="normal"), 
             mapping = aes(x = reorder(sample, HBEGF), y = 2, col="none"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="R108K/G"), 
             mapping = aes(x = reorder(sample, HBEGF), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="G598V/A"), 
             mapping = aes(x = reorder(sample, HBEGF), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D"), 
             mapping = aes(x = reorder(sample, HBEGF), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x = reorder(sample, HBEGF), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x =reorder(sample, HBEGF), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, HBEGF), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, HBEGF), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,vIII_status=="mutated"), 
             mapping = aes(x = reorder(sample, HBEGF), y = 1.20, col="vIII mutated"),size=1.2, pch=15, alpha=1) +  
  geom_point(data = subset(counts_ligands,EGFR>10.71), 
             mapping = aes(x = reorder(sample, HBEGF), y = 1, col="EGFR amplified"),size=1.2, pch=15, alpha=1) + 
  scale_color_manual(values = colors) +
  labs(y="HBEGF expression (log2)", colour="EGFR mutation") +
  facet_grid(~amplified , drop=F, scales="free_x", space = "free_x") +
  job_gg_theme
#ggsave("HBEGF_sorted_mutations.png")
#ggsave("HBEGF_sorted_mutations.pdf")

AREG <- ggplot(counts_ligands, aes(x=reorder(sample, AREG), y = AREG)) +
  geom_point(pch=20) +
  geom_point(data = subset(counts_ligands,mutation=="none"&vIII_status=="normal"), 
             mapping = aes(x = reorder(sample, AREG), y = 2, col="none"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="R108K/G"), 
             mapping = aes(x = reorder(sample, AREG), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="G598V/A"), 
             mapping = aes(x = reorder(sample, AREG), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D"), 
             mapping = aes(x = reorder(sample, AREG), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x = reorder(sample, AREG), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/+R108K/G"), 
             mapping = aes(x =reorder(sample, AREG), y = 1.80, col="R108K/G"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, AREG), y = 1.40, col="A289V/T/D"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,mutation=="A289V/T/D+G598V/A"), 
             mapping = aes(x = reorder(sample, AREG), y = 1.60, col="G598V/A"),size=1.2, pch=15, alpha=1) +
  geom_point(data = subset(counts_ligands,vIII_status=="mutated"), 
             mapping = aes(x = reorder(sample, AREG), y = 1.20, col="vIII mutated"),size=1.2, pch=15, alpha=1) +  
  geom_point(data = subset(counts_ligands,EGFR>10.71), 
             mapping = aes(x = reorder(sample, AREG), y = 1, col="EGFR amplified"),size=1.2, pch=15, alpha=1) + 
  scale_color_manual(values = colors) +
  labs(y="AREG expression (log2)", colour="EGFR mutation") +
  facet_grid(~amplified , drop=F, scales="free_x", space = "free_x") +
  job_gg_theme
#ggsave("AREG_sorted_mutations.png")
#ggsave("AREG_sorted_mutations.pdf")

