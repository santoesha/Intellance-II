#Load data
intellance.expression <- read.delim('~/project/data/intellanceII/Almac/output/tables/gencode-readcounts-103372-001.txt', stringsAsFactors = F, comment.char = '#')
rownames(intellance.expression) <- intellance.expression$Geneid
intellance.expression <- intellance.expression[,-c(2:6)]
colnames(intellance.expression) <- gsub(".Aligned.sortedByCoord.out.bam","",colnames(intellance.expression))
rownames(intellance.expression) <- gsub("\\..*", "", rownames(intellance.expression))
colnames(intellance.expression) <- gsub(".*hg38.","",colnames(intellance.expression))
colnames(intellance.expression) <- gsub(".", "-", colnames(intellance.expression), fixed=TRUE)

gene_ids <- read.csv(file="gene_ids.csv",row.names = c(1))
colnames(intellance.expression)[colnames(intellance.expression) == "Geneid"] <- "ens.id"
intellance.expression <- merge(gene_ids,intellance.expression)
duplicates <- intellance.expression$gene.id[duplicated(intellance.expression$gene.id)]
intellance.expression <- intellance.expression[-which(intellance.expression$gene.id %in% duplicates),]
rownames(intellance.expression) <- intellance.expression$gene.id
raw_counts <- intellance.expression[,-c(1,2)]

rm(gene_ids,intellance.expression,duplicates)
