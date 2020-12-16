#setwd("~/project/IntellanceII/Full_Intellance")

# ---- load packages ----
library(ggplot2)
library(tidyverse)
library(dplyr)

# ---- load data ----
source("Scripts/job_gg_theme.R")

counts <- read.csv("norm_vst_counts_full_intellance.csv",row.names = c(1)) #Load VST transformed count data
colnames(counts) <- gsub('\\X', '', colnames(counts))
colnames(counts) <- gsub('\\.', '-', colnames(counts))
counts <- as.data.frame(t(counts[rownames(counts)=="EGFR|chr7:55M",])) 
colnames(counts) <- "EGFR"
counts$X <- rownames(counts)

#Load vIII data
vIII <- read.delim(file="../../../mnt/neuro-genomic-2/intellanceII/Almac/tables/v3_reads.txt")
vIII$sample <- gsub("/Aligned.*","",vIII$sample)
vIII$vIII_status <- ifelse(vIII$vIII/(vIII$wt+vIII$vIII) > 0.01, "mutated", "normal")
vIII$vIII_status.10p <- ifelse(vIII$vIII/(vIII$wt+vIII$vIII) > 0.1, "mutated", "normal")
vIII$reads <- vIII$vIII + vIII$wt
vIII[vIII$reads<15,c(4)] <- "normal"

#Merge count and vIII data
counts <- dplyr::left_join(counts,vIII[,c(1,4,5)],by=c("X"="sample"))

# ---- amplification status ----
ampli <- read.csv("../PanelBased/CNV_S1495.csv") %>% 
         dplyr::filter(X=="EGFR") 
ampli <- ampli[,-which(colnames(ampli)=="X")]
ampli <- as.data.frame(t(ampli))
ampli$X <- rownames(ampli)
ampli$X <- gsub("\\.","-",ampli$X)

metadata <- read.delim(file="intellance_identifiers_updates.txt")
ampli <- dplyr::inner_join(metadata[,c(4,7)],ampli,by=c("Almac.ID"="X"))
counts <-  dplyr::inner_join(counts,ampli,by=c("X"="GenomeScan.RNA.ID"))

#EGFR high-expr cut-off: 10.71
ggplot(data=counts,aes(x=EGFR,y=V1,col=vIII_status)) + 
  geom_point() +
  labs(x="EGFR VST counts",y="Copy Number Variation") +
  geom_hline(yintercept=5,linetype="dotted") +
  geom_vline(xintercept=10.71,linetype="dotted") +
  annotate("text", x=10.3, y=6.5, label= "(10.71;5)") +
  job_gg_theme

#ggsave("Figures/counts_cnv_egfr_status.pdf")
counts_amplified <- counts[counts$V1 > 5,] 
counts$amplified <- ifelse(counts$V1 > 5,"ampli","nonampli")
colnames(counts_amplified)[6] <- "EGFR_CNV"

wilcox.test(EGFR ~ vIII_status, data = counts_amplified)
wilcox.test(EGFR_CNV ~ vIII_status, data = counts_amplified)

#saveRDS(counts_amplified,"output/CNV_amplified_table.Rds")
