#Input: Almac DNA/RNA IDs
#Output: Intellance IDs

id_conversion <- function(a) {
  metadata <- read.delim(file="intellance_identifiers_updates.txt")
  a <- as.vector(a)
  
  if (sum(is.na(match(a,metadata$Almac.DNA))) != length(a)) { 
    a <- metadata[c(match(a,metadata$Almac.DNA)),c("Almac.DNA","IntellanceID")]  
  } else if (sum(is.na(match(a,metadata$Almac.RNA))) != length(a)) {
    a <- metadata[c(match(a,metadata$Almac.RNA)),c("Almac.RNA","IntellanceID")]
  } 
  else {
    a <- metadata[c(match(a,metadata$GenomeScan.RNA.ID)),c("GenomeScan.RNA.ID","IntellanceID")]
  }
  return(a)
}