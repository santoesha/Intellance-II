#setwd("~/project/IntellanceII/Full_Intellance")

# ---- Load packages ----
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(DESeq2)

#Source files
source('load_full_intellance.R')

#Remove all genes with mean gene count < 3
raw_counts <- raw_counts[rowMeans(raw_counts) > 3,] #left with 31113 genes 

# ---- DESeq2 with random conditions ----
cond <- as.factor(runif(ncol(raw_counts)) > 0.5)
dds <- DESeqDataSetFromMatrix(as.matrix(raw_counts), data.frame(row.names=colnames(raw_counts), cond), design=~cond) %>%
       DESeq(.)

plot(sort(colSums(raw_counts))) #sorted counts
boxplot(raw_counts[,c(1:50)],xlab="Sample ID",ylab="Counts") #before normalization
boxplot(counts(dds, normalized=T)[,c(1:50)], xlab="Sample ID",ylab="Counts") #after normalization

#---- VST ----
vst <- varianceStabilizingTransformation(dds,blind=TRUE)
counts <- as.data.frame(assay(vst)) 

plot(colMeans(counts),colSds(as.matrix(counts)),xlab="Mean sample count",ylab="Standard Deviation")
plot(colMeans(raw_counts),colSds(as.matrix(raw_counts)),xlab="Mean sample count",ylab="Standard Deviation")

#PCA
plotPCA(vst,intgroup="cond")
rv <- rowVars(assay(vst)) 
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pc <- prcomp(t(assay(vst)[select,]))
loadings <- as.data.frame(pc$rotation)
aload <- abs(loadings)
sweep(aload, 2, colSums(aload), "/")
pca <- pc[["x"]]
pca <- pca[,c(1,2)]
pca <- as.data.frame(pca)

#Remove outlier samples 
outliers <- rownames(pca[pca$PC1>50,])
counts <- counts[,-c(which(colnames(counts) %in% outliers))]


#---- plots ----
#Plot sorted EGFR counts
EGFR <- t(counts["EGFR|chr7:55M",])
plot(sort(EGFR,decreasing = FALSE),ylab="Counts",main="EGFR counts") #counts<12 -> not amplified

#Ligands
ligands <- c("EGFR|chr7:55M","AREG|chr4:74M","EREG|chr4:74M","EPGN|chr4:74M","TGFA|chr2:71M","EGF|chr4:110M","HBEGF|chr5:140M","BTC|chr4:75M")
counts_ligands <- as.data.frame(t(counts[rownames(counts)%in%ligands,])) 
colnames(counts_ligands) <- c("AREG","HBEGF","EREG","EGF","EGFR","TGFA","BTC")
counts_ligands_gathered <- counts_ligands %>%
  gather(key = "variable", value = "Ligand",
         -EGFR)
ggplot(counts_ligands_gathered, aes(x = EGFR, y = Ligand)) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  labs(x="EGFR expression (log2)",y="Ligand expression (log2)") +
  geom_point() +
  facet_wrap(~variable)+
  scale_color_viridis_d()+
  stat_cor(method="spearman",color = "red")
#ggsave("output/ligand_correlation_egfr_intellance.pdf")

#EGFR vs FGFR3
FGFR3_vs_EGFR <- as.data.frame(t(counts[rownames(counts)%in%c("EGFR|chr7:55M","FGFR3|chr4:2M"),])) 
colnames(FGFR3_vs_EGFR) <- c("FGFR3", "EGFR")
ggplot(FGFR3_vs_EGFR,aes(x=EGFR,y=FGFR3)) +
  geom_point() +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  labs(x="EGFR expression (log2)",y="FGFR3 expression (log2)") +
  stat_cor(method="spearman",color = "red")

# ---- DESeq2 EGFR as condition -----
raw_counts <- raw_counts[,-which(colnames(raw_counts)%in%outliers)]
EGFR_condition_factor <- ifelse(EGFR<12,"low","high")
condition <- factor(EGFR_condition_factor)
condition <- relevel(condition, ref="low")
coldata_EGFR <- data.frame(row.names=colnames(raw_counts), condition)
countdata_EGFR <- as.matrix(raw_counts)
dds_EGFR <- DESeqDataSetFromMatrix(countdata_EGFR, coldata_EGFR, design=~condition)
dds_EGFR <- DESeq(dds_EGFR)

result_EGFR <- results(dds_EGFR)
resOrdered <- result_EGFR[order(result_EGFR$padj),] %>% data.frame()
sigresult_EGFR <- resOrdered %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble() %>% 
  filter(padj < 0.01 & abs(log2FoldChange) > 0.5)
significant_genes_EGFR_condition <- sigresult_EGFR$gene
write.csv(significant_genes_EGFR_condition, file="significant_genes_EGFR_condition.csv")

plotMA(result_EGFR, alpha = 0.01,main='DESeq2: D. EGFR low vs. high', ylim=c(-2,2))

#DESeq2 EGFR as condition but not in expression data 
data_noegfr <- raw_counts[-which(rownames(raw_counts)=="EGFR|chr7:55M"),]
coldata_noegfr <- data.frame(row.names=colnames(data_noegfr), condition)
countdata_noegfr <- as.matrix(data_noegfr)
dds_noEGFR <- DESeqDataSetFromMatrix(countdata_noegfr, coldata_noegfr, design=~condition)
dds_noEGFR <- DESeq(dds_noEGFR)
result_noEGFR <- results(dds_noEGFR)

vst_noegfr <- varianceStabilizingTransformation(dds_noEGFR)
counts_noegfr <- assay(vst_noegfr)
plotPCA(vst_noegfr)

resOrdered <- result_noEGFR[order(result_noEGFR$padj),] %>% as.data.frame()

EnhancedVolcano(title = NULL,
                subtitle = NULL,
                result_noEGFR,
                lab = rownames(result_noEGFR),
                axisLabSize = 12,
                legendPosition = "right",
                legendLabSize = 9,
                x = 'log2FoldChange',
                y = 'pvalue',  
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                col = c("black", "forestgreen", "blue", "red"),
                xlim = c(-8,8), 
                captionLabSize = 5)