#setwd("~/project/IntellanceII/Full_Intellance")

# ---- load packages ----
library(tidyverse)
library(VariantAnnotation)
library(stringr)

# ---- Load data ----
path <- "../../../mnt/neuro-genomic-1/intellanceII/Almac/ADX17020_Client_Folder-48243198/ADX17020_Client_Folder-48243198"
list <- list.files(path=path,recursive=T,pattern = "_SmallVariants.genome.vcf$",full.names = T) 
list <- list[-c(3,162,65,86)] #replicates, and sample with only 700 reads
tmp_vcf_data <- lapply(list,read.table,stringsAsFactors=FALSE)

#Filter EGFR and PASS
tmp_vcf_data_EGFR <- lapply(tmp_vcf_data, function(x) filter(x, V1 == "chr7" & V2 %in% c(55086714:55324313) & V5 != "." & V7 == "PASS"))

#readRDS(file="tmp_vcf_data_EGFR.RData")

# ---- Convert Sample IDs ----
sample_names <- gsub(".*Folder-48243198/S1495","S1495",list)  
sample_names <- gsub("/S1495.*","",sample_names)
names(tmp_vcf_data_EGFR) <- sample_names

source("sample_id_conversion.R")
intellance_id <- names(tmp_vcf_data_EGFR)
intellance_id <- id_conversion(intellance_id)
names(tmp_vcf_data_EGFR) <- intellance_id$IntellanceID

mutations <- lapply(tmp_vcf_data_EGFR, function (x) x[c('V8','V10')])
index <- sapply(mutations, length)
res <- as.data.frame(do.call(rbind,lapply(mutations, `length<-`,
                                          max(index))))
res$sample <- rownames(res)

# ---- Extract mutations and VAF ----
res$V8 <- gsub("EGFR.*","",res[,"V8"])

missense <- res[grep(pattern="missense_variant",res$V8),]
frameshift <- res[grep(pattern="frameshift_variant",res$V8),]

missense$V8 <- gsub("EGFR.*","",missense[,"V8"])
missense$V8 <- gsub(".*:p","",missense[,"V8"])
missense$V8 <- str_replace(missense$V8, ".", "")
missense$V10 <- gsub(":20:.*","",missense[,"V10"])
missense$V10 <- gsub(".*:0.","0.",missense[,"V10"])
names(missense) <- c("missense_mutation","VAF","sample")

frameshift$V8 <- gsub("EGFR.*","",frameshift[,"V8"])
frameshift$V8 <- gsub(".*:p","",frameshift[,"V8"])
frameshift$V8 <- str_replace(frameshift$V8, ".", "")
frameshift$V10 <- gsub(":20:.*","",frameshift[,"V10"])
frameshift$V10 <- gsub(".*:0.","0.",frameshift[,"V10"])
names(frameshift) <- c("frameshift","VAF","sample")

#Interesting point mutations 
R108 <- missense %>% filter(grepl("Arg108",missense_mutation))
A289 <- missense %>% filter(grepl("Ala289",missense_mutation))
G598 <- missense %>% filter(grepl("Gly598",missense_mutation))

#readRDS(file="SNV_R108_full.RData")
#readRDS(file="SNV_A289_full.RData")
#readRDS(file="SNV_G598_full.RData")
