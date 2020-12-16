#setwd("~/project/IntellanceII/Full_Intellance/RF_script")

# ---- Load packages ----
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)

# ---- Model: Fold 1 ----
#Performance
performance_fold1 <- read.delim("performance.txt",header = FALSE, sep=" ")
training_perf_fold1 <- mean(performance_fold1$V1)
validation_perf_fold1 <- mean(performance_fold1$V2)

#Boruta
boruta_fold1 <- read.delim("boruta.txt") %>% separate(col = genes, into = c("number", "gene"), sep = " ") %>% drop_na() %>% dplyr::select(-c(1)) %>% dplyr::group_by(gene) %>% summarise(count = n()) 
number_of_genes_fold1 <- read.delim("number_of_genes.txt",header=FALSE)
number_of_genes_fold1 <- mean(number_of_genes_fold1$V1)

#Varimp
varimp_fold1 <- read.delim("varimp.txt") %>% separate(col = Overall, into = c("gene", "score"), sep = " ") %>% drop_na() 
varimp_fold1$score <- as.numeric(varimp_fold1$score) 
varimp_fold1 <- varimp_fold1 %>% dplyr::group_by(gene) %>% summarise(count_fold1=n(),mean = mean(score))
varimp_fold1$weighted_fold1 <- varimp_fold1$count_fold1 * varimp_fold1$mean

# ---- Model: Fold 2 ----
#Performance
performance_fold2 <- read.delim("performance_fold2.txt",header = FALSE, sep=" ")
training_perf_fold2 <- mean(performance_fold2$V1)
validation_perf_fold2 <- mean(performance_fold2$V2)

#Boruta
boruta_fold2 <- read.delim("boruta_fold2.txt") %>% separate(col = genes, into = c("number", "gene"), sep = " ") %>% drop_na() %>% dplyr::select(-c(1)) %>% dplyr::group_by(gene) %>% summarise(count = n()) 
number_of_genes_fold2 <- read.delim("number_of_genes_fold2.txt",header=FALSE)
number_of_genes_fold2 <- mean(number_of_genes_fold2$V1)

#Varimp
varimp_fold2 <- read.delim("varimp_fold2.txt") %>% separate(col = Overall, into = c("gene", "score"), sep = " ") %>% drop_na() 
varimp_fold2$score <- as.numeric(varimp_fold2$score) 
varimp_fold2 <- varimp_fold2 %>% dplyr::group_by(gene) %>% summarise(count_fold2=n(),mean = mean(score))
varimp_fold2$weighted_fold2 <- varimp_fold2$count_fold2 * varimp_fold2$mean

# ---- Model: Fold 3 ----
#Performance
performance_fold3 <- read.delim("performance_fold3.txt",header = FALSE, sep=" ")
training_perf_fold3 <- mean(performance_fold3$V1)
validation_perf_fold3 <- mean(performance_fold3$V2)

#Boruta
boruta_fold3 <- read.delim("boruta_fold3.txt") %>% separate(col = genes, into = c("number", "gene"), sep = " ") %>% drop_na() %>% dplyr::select(-c(1)) %>% dplyr::group_by(gene) %>% summarise(count = n()) 
number_of_genes_fold3 <- read.delim("number_of_genes_fold3.txt",header=FALSE)
number_of_genes_fold3 <- mean(number_of_genes_fold3$V1)

#Varimp
varimp_fold3 <- read.delim("varimp_fold3.txt") %>% separate(col = Overall, into = c("gene", "score"), sep = " ") %>% drop_na() 
varimp_fold3$score <- as.numeric(varimp_fold3$score) 
varimp_fold3 <- varimp_fold3 %>% dplyr::group_by(gene) %>% summarise(count_fold3=n(),mean = mean(score))
varimp_fold3$weighted_fold3 <- varimp_fold3$count_fold3 * varimp_fold3$mean

# ---- Average Model ----
train_performance <- mean(c(performance_fold1[1:30,"V1"],performance_fold2[1:30,"V1"],performance_fold3[1:30,"V1"])) #trainingperf
train_sd <- sd(c(performance_fold1[1:30,"V1"],performance_fold2[1:30,"V1"],performance_fold3[1:30,"V1"])) 

validation_performance <- mean(c(performance_fold1[1:30,"V2"],performance_fold2[1:30,"V2"],performance_fold3[1:30,"V2"])) #trainingperf
validation_sd <- sd(c(performance_fold1[1:30,"V2"],performance_fold2[1:30,"V2"],performance_fold3[1:30,"V2"])) 

selected_genes <- data.frame(gene=c(varimp_fold1$gene,varimp_fold2$gene,varimp_fold3$gene) %>% unique()) 
selected_genes <- dplyr::inner_join(selected_genes,varimp_fold1[,c(1,4)],by="gene")
selected_genes <- dplyr::inner_join(selected_genes,varimp_fold2[,c(1,4)],by="gene")
selected_genes <- dplyr::inner_join(selected_genes,varimp_fold3[,c(1,4)],by="gene")
selected_genes$mean_scores <- rowMeans(selected_genes[,-1])
selected_genes$mean_scores <- selected_genes$mean_scores/max(selected_genes$mean_scores)

top30 <- selected_genes[order(selected_genes$mean_scores, decreasing=TRUE),] %>% dplyr::slice(1:30)

ggplot(data=top30, aes(x=reorder(gene,mean_scores), y=mean_scores)) +
  geom_bar(stat="identity") + 
  coord_flip() +
  labs(x="Gene",y="Variable importance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=20),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#ggsave("../Figures/intellanceII_featureimportance.pdf",width=3*7.67 / 1.45, height = 6.83 / 1.45 * 2)
W