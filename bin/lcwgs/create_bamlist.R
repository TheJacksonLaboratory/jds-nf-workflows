#!/usr/bin/env Rscript
string <- read.table("bamlist.txt")
bamlist <- as.matrix(t(string))
rownames(bamlist) <- NULL
for(i in 1:nrow(bamlist)){
  input <- bamlist[i,]
  no_bracket <- gsub("\\[", "", input)
  no_bracket <- gsub("\\]", "", no_bracket)
  no_comma <- gsub(".bam,", ".bam", no_bracket)
  bamlist[i,] <- no_comma
}
bamlist[,1] <- bamlist[,1][grep(bamlist[,1], pattern = ".bam")]
bamlist <- bamlist[!duplicated(bamlist[,1]),]
colnames(bamlist) <- NULL
write.table(bamlist, file = "bamlist.txt", quote = F, row.names = F, col.names = F)
