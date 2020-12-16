#setwd("~/project/IntellanceII/Full_Intellance/RF_script")

# ---- load packages ----
library(VennDiagram)
library(corrplot)
library(RColorBrewer)
library(tidyverse)
library(dplyr)

# ---- Load data ----
#Load counts
counts_intellance <- read.csv("../norm_vst_counts_full_intellance.csv",row.names = c(1))
colnames(counts_intellance) <- gsub('\\X', '', colnames(counts_intellance))
names(counts_intellance) <- gsub('\\.', '-', colnames(counts_intellance))

counts_tcga <- read.csv("~/project/IntellanceII/TCGA/vst_normalized_tcga_gbm_counts.csv",row.names = c(1))
names(counts_tcga) <- gsub('\\.', '-', colnames(counts_tcga))

#Variable importance TCGA and Intellance II
varimp_intellance <- read.csv("varimp_intellance.csv")
counts_intellance <- counts_intellance[which(rownames(counts_intellance)%in%varimp_intellance$gene),]
rownames(counts_intellance) <- gsub("\\|.*$","",rownames(counts_intellance))
varimp_intellance$gene <- gsub("\\|.*$","",varimp_intellance$gene)

varimp_tcga <- read.csv("~/project/IntellanceII/TCGA/RF_script/varimp_tcga.csv")
counts_tcga <- counts_tcga[which(rownames(counts_tcga)%in%varimp_tcga$gene),]

# ---- Venn Diagram ----
v <- venn.diagram(
  x = list(varimp_intellance$gene, varimp_tcga$gene),
  category.names = c("IntellanceII" , "TCGA"),
  filename = NULL,
  output=TRUE,
  col=c("chocolate3", '#21908dff'),
  fill = c(alpha("chocolate3",0.3), alpha('#21908dff',0.3)),
  lwd = 2,
  cex = c(3,3,3),
  cat.cex = 3,
  cat.default.pos = "outer",
  cat.fontface = "bold",
  cat.pos = c(-5,5),
  ext.pos = 20
)
#v[[7]]$label <- paste(intersect(varimp_intellance$gene, varimp_tcga$gene), collapse="\n")

#pdf("../Figures/venndiagram_tcga_intellance.pdf",width=3*7.67 / 1.45, height = 6.83 / 1.45 * 2)
#grid.draw(v)
#dev.off()

intersect <- intersect(varimp_intellance$gene,varimp_tcga$gene)
length(intersect)

# ---- Correlation plot ----
counts_tcga_intersect <- counts_tcga[which(rownames(counts_tcga)%in%intersect),]
counts_intellance_intersect <- counts_intellance[which(rownames(counts_intellance)%in%intersect),]

top30_intellance <- varimp_intellance[order(varimp_intellance$mean_scores, decreasing=TRUE),] %>%
                    dplyr::slice(1:30)
counts_top30_intellance <- counts_intellance[which(rownames(counts_intellance)%in%top30_intellance$gene),]

M <- cor(t(counts_top30_intellance),method="spearman")
corrplot(M, method="ellipse", type="upper",order="hclust", tl.col="black",col=brewer.pal(n=8, name="PuOr"),tl.cex=0.7)

#pdf("../Figures/corrplot_informative_genes_RF.pdf",width=3*7.67 / 1.45, height = 6.83 / 1.45 * 2)
M2 <- cor(t(counts_intellance_intersect),method="spearman")
corrplot(M2, method="ellipse", type="upper",order="hclust", tl.col="black",col=brewer.pal(n=8, name="PuOr"),tl.cex=1.5)
#dev.off()
