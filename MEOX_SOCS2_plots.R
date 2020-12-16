#setwd("~/project/IntellanceII/Full_Intellance")

# ---- load packages ----
library(ggpubr)
library(ggplot2)
library(tidyverse)

# ---- load data ----
source('R/job_gg_theme.R')

#Load VST transformed count data
counts <- read.csv("norm_vst_counts_full_intellance.csv",row.names = c(1))
colnames(counts) <- gsub('\\X', '', colnames(counts))
colnames(counts) <- gsub('\\.', '-', colnames(counts))

#Load vIII data
vIII <- read.delim(file="~/project/data/intellanceII/Almac/tables/v3_reads.txt")
vIII$sample <- gsub("/Aligned.*","",vIII$sample)
vIII$vIII_status <- ifelse(vIII$vIII/(vIII$wt+vIII$vIII) > 0.01, "mutated", "normal")
vIII$reads <- vIII$vIII + vIII$wt
vIII[vIII$reads<15,c(4)] <- "normal"

#---- MEOX2 ----
counts_meox <- as.data.frame(t(counts[rownames(counts)%in%c("MEOX2|chr7:16M","EGFR|chr7:55M"),])) 

#Merge count and vIII data
counts_meox$sample <- rownames(counts_meox)
counts_meox <- merge(vIII[,c(1,4)],counts_meox,by="sample")

job_gg_theme <- theme(
  text = element_text(family = 'Helvetica'),
  legend.position = 'bottom',
  plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
  panel.background = element_rect(fill = 'white', colour = 'white'),
  axis.title = element_text(face = "bold",size = rel(1)),
  axis.title.y = element_text(angle=90,vjust =2),
  axis.text = element_text(),
  axis.ticks.x=element_blank(),
  axis.line = element_line(colour="black"),
)

ggplot(counts_meox, aes(x=`MEOX2|chr7:16M`, y = `EGFR|chr7:55M`, col=vIII_status)) +
  geom_point(pch=20) +
  labs(x="MEOX2",y="EGFR") +
  job_gg_theme

# ---- SOCS2 ----
counts_socs <- as.data.frame(t(counts[rownames(counts)%in%c("SOCS2|chr12:94M","EGFR|chr7:55M"),])) 

counts_socs$sample <- rownames(counts_socs)
counts_socs <- merge(vIII[,c(1,4)],counts_socs,by="sample")

job_gg_theme <- theme(
  text = element_text(family = 'Helvetica'),
  legend.position = 'bottom',
  plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
  panel.background = element_rect(fill = 'white', colour = 'white'),
  axis.title = element_text(face = "bold",size = rel(1)),
  axis.title.y = element_text(angle=90,vjust =2),
  axis.text = element_text(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.line = element_line(colour="black"),
)

socs <- ggplot(subset(counts_socs,`EGFR|chr7:55M`>10.71), aes(x=reorder(sample, `SOCS2|chr12:94M`), y = `SOCS2|chr12:94M`, col=vIII_status)) +
  geom_point(pch=20) +
  labs(x="Index",y="SOCS2") +
  job_gg_theme #only amplified samples

job_gg_theme <- theme(
  text = element_text(family = 'Helvetica'),
  legend.position = 'bottom',
  plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
  panel.background = element_rect(fill = 'white', colour = 'white'),
  axis.title = element_text(face = "bold",size = rel(1)),
  axis.title.y = element_text(angle=90,vjust =2),
  axis.text = element_text(),
  axis.ticks.x=element_blank(),
  axis.line = element_line(colour="black"),
)

ggplot(counts_socs, aes(x=`EGFR|chr7:55M`, y = `SOCS2|chr12:94M`, col=vIII_status)) +
  geom_point(pch=20) +
  labs(x="EGFR VST counts",y="SOCS2 VST counts") +
  stat_cor(method="spearman",color = "red") +
  job_gg_theme
#ggsave("output/SOCS2_vs_EGFR_plot.pdf")
