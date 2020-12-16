# ---- Load packages ----
library(dplyr)
library(tidyverse)
library(ggplot2)

# ---- Load data ----

#G-SAM
setwd("~/project/G-SAM")
source("~/project/G-SAM/scripts/R/gsam_metadata.R")
source("scripts/R/snv.R")

selected_mut <- c("G598V","G598A","A289V","A289T","A289D","R108K","R108G")
vaf$selected_mut <-  ifelse(rownames(vaf)%in%grep(paste(selected_mut,collapse="|"),vaf$egfr_snv),"mutated","none")

#IntellanceII
setwd("~/project/IntellanceII/Full_Intellance")
egfr_snv <- read.delim("output/snv_fullintellance.txt", sep = " ")

source("sample_id_conversion.R")
metadata <- read.delim(file="intellance_identifiers_updates.txt")

source("Scripts/job_gg_theme.R")

# ---- MDM2, CDK4 and CDNK2A plot ----
ampli_genes <- read.csv("../PanelBased/CNV_S1495.csv") %>% 
               dplyr::filter(X=="CDK4" | X =="MDM2" | X == "CDKN2A" | X == "EGFR") %>%
               t(.) 

colnames(ampli_genes) <- ampli_genes[1,] 
ampli_genes <- as.data.frame(ampli_genes[-1,]) 

ampli_genes$X <- rownames(ampli_genes) %>%
                 gsub("\\.","-",.)

ampli_genes <- dplyr::left_join(ampli_genes,metadata[,c(1,4,7)],by=c("X"="Almac.ID")) %>%
               dplyr::left_join(.,egfr_snv,by=c("IntellanceID"="sample")) %>%
               dplyr::mutate(CDK4=as.numeric(as.character(CDK4))) %>%
               dplyr::mutate(MDM2=as.numeric(as.character(MDM2))) %>%
               dplyr::mutate(CDKN2A=as.numeric(as.character(CDKN2A))) %>%
               dplyr::mutate(EGFR=as.numeric(as.character(EGFR))) %>%
               dplyr::mutate(EGFR_ampli = ifelse(EGFR >= 5,T,F)) %>%
               dplyr::mutate(CDK4_ampli = ifelse(CDK4 >= 5,T,F)) %>%
               dplyr::mutate(MDM2_ampli = ifelse(MDM2 >= 5,T,F)) %>%
               dplyr::mutate(CDKN2A_deletion = ifelse(CDKN2A < 1,T,F)) %>%
               dplyr::filter(EGFR_ampli == T) %>%
               dplyr::mutate(selected_mut = ifelse(mutation == "none","normal","mutated")) %>%
               dplyr::mutate(CDK4orMDM2 = ifelse(CDK4_ampli|MDM2_ampli,"CDK4/MDM2","none")) 

ampli_plot <-  ampli_genes %>%
               dplyr::group_by(selected_mut,CDK4orMDM2) %>%
               dplyr::summarise(count=n()) %>%
               dplyr::filter(!is.na(selected_mut))
fisher <- fisher.test(ampli_genes$CDK4orMDM2,ampli_genes$selected_mut)

ggplot(ampli_plot, aes(fill=CDK4orMDM2, y=count, x=selected_mut)) + 
       geom_bar(position="stack", stat="identity") +
       xlab("EGFR extracellular missense mutatation") +
       ylab("Number of samples") +
       scale_fill_grey() +
       job_gg_theme +
       theme(legend.title = element_blank()) +
       annotate("text", x = 1 , y = 120, label = paste0("p = ",round(fisher$p.value,3)), col="black" , size = 5)
#ggsave("output/CDK4_MDM2_plot.pdf",width=6,height=5)
