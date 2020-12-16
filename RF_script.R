#setwd("~/project/IntellanceII/Full_Intellance/RF_script")

# ---- Load packages ----
library(ranger)
library(caret)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(Boruta)

# ---- Load data ----
counts <- read.csv("../norm_vst_counts_full_intellance.csv",row.names = c(1))
colnames(counts) <- gsub('\\X', '', colnames(counts))
names(counts) <- gsub('\\.', '-', colnames(counts))

# ---- Data partition: train and test set ----
set.seed(10)
counts_RF <- counts
rownames(counts_RF)[which(rownames(counts_RF)=="EGFR|chr7:55M")] <- "EGFR"
counts_RF <- as.data.frame(t(counts_RF))
intrain <- createDataPartition(y = counts_RF$EGFR, p=0.85, list=FALSE)
 
train_set <- counts_RF[intrain,]
test_set <- counts_RF[-intrain,]

#write(rownames(test_set),"test_set_samples.txt")
train_set <- as.data.frame(t(train_set))

# ---- 3 times repeated 10FCV ----
for (i in 1:3) {
  flds <- createFolds(train_set, k = 10, list = TRUE, returnTrain = FALSE)
  for(i in 1:10) { 
    train_set_fold <- train_set[,-flds[[i]]] 
    validation_set <- train_set[,flds[[i]]] 
    train_set_fold <- as.data.frame(t(train_set_fold))
                 
    #Filtering genes by: correlation 
    cor_mat <- cor(scale(train_set_fold),method="spearman")
    cor_EGFR <- cor_mat[,which(colnames(cor_mat)=="EGFR")]
    cor_EGFR <- as.data.frame(cor_EGFR)
    cor_EGFR <- cor_EGFR %>% dplyr::filter(abs(cor_EGFR)>0.25) 
    train_set_fold <- train_set_fold[,which(colnames(train_set_fold)%in%rownames(cor_EGFR))]
  
    #Filtering genes by: Median Absolute Deviation
    train_set_fold <- as.data.frame(t(train_set_fold))
    train_set_fold$mad <- apply(train_set_fold, 1, mad)
    train_set_fold <- train_set_fold %>% dplyr::filter(train_set_fold$mad > 0.5) 
    train_set_fold$mad <- NULL
    train_set_fold <- as.data.frame(t(train_set_fold))
  
    EGFR <-  subset(train_set_fold,select = c("EGFR"))
  
    #Feature selection: Boruta
    boruta <- Boruta(EGFR~.,data=train_set_fold)
    print(boruta)
    boruta <- as.data.frame(boruta$finalDecision)  
    boruta <- boruta %>% dplyr::filter(boruta$`boruta$finalDecision` %in% c("Confirmed","Tentative"))
    genes <- as.data.frame(rownames(boruta))
    genes <- as.data.frame(gsub("`","",genes$`rownames(boruta)`))
    colnames(genes)[1] <- "genes"
    
    #Filter neighboring genes
    genes <- genes %>% dplyr::filter(!genes %in% c("EGFR-AS1|chr7:55M"))#,"ELDR|chr7:55M","AC006971.1|chr7:55M","AC074351.1|chr7:55M","AC073324.1|chr7:55M","SEC61G|chr7:55M","LANCL2|chr7:55M","AC073347.1|chr7:55M"))
    train_set_fold <- train_set_fold[,which(colnames(train_set_fold)%in%genes$genes)]
    
    fit_ct <- caret::train(
              x = train_set_fold,
              y = EGFR$EGFR,
              method = "ranger",
              importance = "permutation")
    print(fit_ct)
       
    var_imp <- varImp(fit_ct)
    
    #Train set predictions
    actual_train <- EGFR$EGFR
    predicted_train <- fit_ct$finalModel$predictions
    R2_train <- 1 - (sum((actual_train-predicted_train)^2)/sum((actual_train-mean(actual_train))^2))

    #Validation set prediction
    validation_set <- as.data.frame(t(validation_set))
    actual_validation <- validation_set$EGFR
    predicted_validation <- unname(predict(fit_ct, validation_set))
    R2_val <- 1 - (sum((actual_validation-predicted_validation)^2)/sum((actual_validation-mean(actual_validation))^2))
    # 
    # #Write results
    # write(t(c(R2_train,R2_val)), 'performance.txt', append=TRUE)
    # write.table(genes, 'boruta.txt', append=TRUE)
    # write.table(var_imp$importance, 'varimp.txt', append=TRUE)
    # write(nrow(genes),"number_of_genes.txt",append=TRUE)
    }
}

# ---- Combine results  ----
#Performance
performance <- read.delim("performance.txt",header = FALSE, sep=" ")
training_perf <- mean(performance$V1)
validation_perf <- mean(performance$V2)

c(mean(performance[1:10,"V1"]),mean(performance[11:20,"V1"]),mean(performance[21:30,"V1"])) #standard deviation
c(sd(performance[1:10,"V1"]),sd(performance[11:20,"V1"]),sd(performance[21:30,"V1"])) #standard deviation

#Boruta
boruta <- read.delim("boruta.txt") %>% separate(col = genes, into = c("number", "gene"), sep = " ") %>% drop_na() %>% dplyr::select(-c(1)) %>% dplyr::group_by(gene) %>% summarise(count = n()) 
number_of_genes <- read.delim("number_of_genes.txt",header=FALSE)

c(mean(number_of_genes[1:10,"V1"]),mean(number_of_genes[11:20,"V1"]),mean(number_of_genes[21:30,"V1"])) #standard deviation
c(sd(number_of_genes[1:10,"V1"]),sd(number_of_genes[11:20,"V1"]),sd(number_of_genes[21:30,"V1"])) #standard deviation

number_of_genes <- mean(number_of_genes[1:30,1])

#Varimp
varimp <- read.delim("varimp.txt") %>% separate(col = Overall, into = c("gene", "score"), sep = " ") %>% drop_na() 
varimp$score <- as.numeric(varimp$score) 
varimp <- varimp %>% dplyr::group_by(gene) %>% summarise(count=n(),mean = mean(score))
varimp$weighted <- varimp$count * varimp$mean

# ----Final predictor ----
train_set <- counts_RF[intrain,]

varimp <- varimp[order(varimp$weighted, decreasing=TRUE),] %>% dplyr::slice(1:round(number_of_genes))
EGFR <- as.data.frame(train_set[,which(colnames(train_set)=="EGFR")]) 
colnames(EGFR)[1] <- "EGFR"
final_train <- as.data.frame(train_set[,which(colnames(train_set)%in%varimp$gene)]) 

final_fit <- caret::train(
  x = final_train,
  y = EGFR$EGFR,
  method = "ranger",
  importance = "permutation")
print(final_fit)

final_var_imp <- varImp(final_fit)
plot(final_var_imp,top=20)

# ---- test set performance ----
actual <- test_set$EGFR
test_set <- test_set[,which(colnames(test_set)%in%varimp$gene)]

predicted <- unname(predict(final_fit, test_set))
R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))