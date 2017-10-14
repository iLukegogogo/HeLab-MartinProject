require(dplyr)
require(foreach)
mm10.refseq.protein.coding.promoter              <- read.delim("~/Project/Martin/peak.overlapping/mm10.refseq.protein.coding.promoter.bed", header=FALSE, stringsAsFactors=FALSE)
colnames(mm10.refseq.protein.coding.promoter)[4] <- 'gene.name'
the.2c.gene.name                                 <- read.table("~/Project/Martin/2c.gene.name.full.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
flag            <- mm10.refseq.protein.coding.promoter$gene.name %in% the.2c.gene.name
region          <- mm10.refseq.protein.coding.promoter[flag,]
convert.to.gtf.line <- function(df){
    gtf.line <- sprintf("%s\tpseduo.RT\texon\t%d\t%d\t.\t+\t.\trepeat \"pseduo.RT.2C\"; gene_id \"%s\"",df$V1,df$V2,df$V3,df$gene.name)  
}

rs <- foreach(i= 1:nrow(region),.combine='c') %do% {
   convert.to.gtf.line(region[i,])  
}
write.table(x=rs,file="pseduo.RT.gtf",quote=FALSE,row.names=FALSE,col.names=FALSE)
