require(plyr)
require(dplyr)
require(ggplot2)
require(foreach)
require(edgeR)


###########################
# Mouse chip-seq data analysis. Given a TF, rank all RT according to the fraction of reads coming from the RT
###########################
load('RData//mouse.chip.seq.csv.chip.seq.results.RData')
#load('RData//mouse.batch1.csv.chip.seq.results.RData') the previous results were genearted by this RData file

exp.meta             <- ddply(exp.meta,.(ip),function(x) x[1,])
exp.meta$exp.name <- gsub(x=exp.meta$exp.name,pattern="\\.\\d+",replacement = "",perl=TRUE)
ip.srr.id            <- exp.meta$ip %>% as.character 
ip.srr.bam           <- paste(ip.srr.id,".all.rmdup.bam",sep='')
rt.read.count.matrix <- rt.read.count.matrix[,ip.srr.bam]

rt.fraction.matrix.v1 <- foreach(j = 1:ncol(rt.read.count.matrix),.combine='cbind') %do% {
    f <- rt.read.count.matrix[,j]  / lib.size[colnames(rt.read.count.matrix)[j]]
    f
}
colnames(rt.fraction.matrix.v1) <- colnames(rt.read.count.matrix)

rt.fraction.matrix.v2 <- foreach(j = 1:ncol(rt.read.count.matrix),.combine='cbind') %do% {
    f <- rt.read.count.matrix[,j]  / sum(rt.read.count.matrix[,j]) 
    f
}
colnames(rt.fraction.matrix.v2) <- colnames(rt.read.count.matrix)


rt.TF.ranking.v1 <- foreach(i = 1:nrow(rt.fraction.matrix.v1),.combine='rbind') %do% {
    row        <- rt.fraction.matrix.v1[i,]
    row        <- scale(row,center = TRUE,scale = TRUE)
    names(row) <- paste(exp.meta$exp.name,row,exp.meta$ip,sep="_")
    row        <- sort(row,decreasing = TRUE)
    names(row)
}
rownames(rt.TF.ranking.v1) <- rownames(rt.fraction.matrix.v1)

rt.TF.ranking.v2 <- foreach(i = 1:nrow(rt.fraction.matrix.v2),.combine='rbind') %do% {
  row        <- rt.fraction.matrix.v2[i,]
  row        <- scale(row,center = TRUE,scale = TRUE)
  names(row) <- paste(exp.meta$exp.name,row,exp.meta$ip,sep="_")
  row        <- sort(row,decreasing = TRUE)
  names(row)
}
rownames(rt.TF.ranking.v2) <- rownames(rt.fraction.matrix.v2)

write.csv(x=rt.TF.ranking.v1 %>% t,file="Results/mouse.rt.TF.ranking.matrix.v1.csv",row.names=FALSE,quote=FALSE)
write.csv(x=rt.TF.ranking.v2 %>% t,file="Results/mouse.rt.TF.ranking.matrix.v2.csv",row.names=FALSE,quote=FALSE)

########################
#####
#######################

df <- foreach(i = 1:nrow(rt.fraction.matrix.v1),.combine='rbind') %do% {
  row        <- rt.fraction.matrix.v1[i,]
  row        <- scale(row,center = TRUE,scale = TRUE)
  TF.name    <- exp.meta$exp.name[row >=3]
  rt.name    <- rownames(rt.fraction.matrix.v1)[i]
  data.frame(rt.name=rt.name,TF.name=TF.name)
}
df <- unique(df)
df <- df[is.na(df$TF.name)==FALSE,]

TF.statistics <- ddply(df,.(TF.name),nrow) %>% arrange(desc(V1))
rt.statistics <- ddply(df,.(rt.name),nrow) %>% arrange(desc(V1))
rownames(TF.statistics) <- TF.statistics$TF.name
rownames(rt.statistics) <- rt.statistics$rt.name

MT2.mm.TF <- df$TF.name[df$rt.name == 'MT2_Mm']
TF.statistics$is.MT2.Mm <- (TF.statistics$TF.name %in% MT2.mm.TF) %>% as.character
ggplot(TF.statistics,aes(x=TF.name,y=V1,fill = factor(is.MT2.Mm))) + geom_bar(stat="identity")+coord_flip()+theme(axis.text=element_text(size=10))


MERVL.int.TF <- df$TF.name[df$rt.name == 'MERVL-int']
TF.statistics$is.MERVL.int.TF <- (TF.statistics$TF.name %in% MERVL.int.TF) %>% as.character
ggplot(TF.statistics,aes(x=TF.name,y=V1,fill = factor(is.MERVL.int.TF))) + geom_bar(stat="identity")+coord_flip()+theme(axis.text=element_text(size=10))


################################
# Human chip-seq data analysis. Given a TF, rank all RT according to the fraction of reads coming from the RT
#################################
load('RData//human.chip.seq.csv.chip.seq.results.RData')
#load('RData//human.batch1.csv.chip.seq.results.RData') previous version 
exp.meta             <- ddply(exp.meta,.(ip),function(x) x[1,])
exp.meta$exp.name <- gsub(x=exp.meta$exp.name,pattern="\\.\\d+",replacement = "",perl=TRUE)
ip.srr.id            <- exp.meta$ip %>% as.character
ip.srr.bam           <- paste(ip.srr.id,".all.rmdup.bam",sep='')
rt.read.count.matrix <- rt.read.count.matrix[,ip.srr.bam]

rt.fraction.matrix.v1 <- foreach(j = 1:ncol(rt.read.count.matrix),.combine='cbind') %do% {
  f <- rt.read.count.matrix[,j]  / lib.size[colnames(rt.read.count.matrix)[j]]
  f
}

rt.fraction.matrix.v2 <- foreach(j = 1:ncol(rt.read.count.matrix),.combine='cbind') %do% {
  f <- rt.read.count.matrix[,j]  / sum(rt.read.count.matrix[,j]) 
  f
}


rt.TF.ranking.v1 <- foreach(i = 1:nrow(rt.fraction.matrix.v1),.combine='rbind') %do% {
  row        <- rt.fraction.matrix.v1[i,]
  row        <- scale(row,center = TRUE,scale = TRUE)
  names(row) <- paste(exp.meta$exp.name,row,exp.meta$ip,sep="_")
  row        <- sort(row,decreasing = TRUE)
  names(row)
}
rownames(rt.TF.ranking.v1) <- rownames(rt.fraction.matrix.v1)

rt.TF.ranking.v2 <- foreach(i = 1:nrow(rt.fraction.matrix.v2),.combine='rbind') %do% {
  row        <- rt.fraction.matrix.v2[i,]
  row        <- scale(row,center = TRUE,scale = TRUE)
  names(row) <- paste(exp.meta$exp.name,row,exp.meta$ip,sep="_")  
  row        <- sort(row,decreasing = TRUE)
  names(row)
}
rownames(rt.TF.ranking.v2) <- rownames(rt.fraction.matrix.v2)

write.csv(x=rt.TF.ranking.v1 %>% t,file="Results/human.rt.TF.ranking.matrix.v1.csv",row.names=FALSE,quote=FALSE)
write.csv(x=rt.TF.ranking.v2 %>% t,file="Results/human.rt.TF.ranking.matrix.v2.csv",row.names=FALSE,quote=FALSE)


##############################################################
# Dux RNA-Seq data to show dosage-depdent global RT induction
###############################################################
load('RData//Dux.dox.csv.rna.seq.results.RData')
Dux.id <- 'ENSMUSG00000075046.6'

#GFP pos, high dose Dux over-expression
gene.data <- gene.read.count.matrix[,c('SRR4032348.bam','SRR4032349.bam','SRR4032350.bam','SRR4032351.bam')]
rt.data   <- rt.read.count.matrix[,c('SRR4032348.bam','SRR4032349.bam','SRR4032350.bam','SRR4032351.bam')]  
de.data   <- rbind(gene.data,rt.data)
flag      <- apply(de.data,1,function(x) sum(x>=10) >=4)
de.data   <- de.data[flag,]

edgeR.res <- DGEList(counts=round(de.data), 
                     group=factor(c('Wild.type','Wild.type','Dux.OE.GFP.pos','Dux.OE.GFP.pos'),levels = c('Wild.type','Dux.OE.GFP.pos'))
)
edgeR.res <- calcNormFactors(edgeR.res, method="RLE")

edgeR.res <- estimateCommonDisp(edgeR.res)
edgeR.res <- estimateTagwiseDisp(edgeR.res)
et        <- exactTest(edgeR.res)
de.res    <- topTags(et, n=Inf)$table
rt.table  <- de.res[!grepl(x = rownames(de.res),pattern = 'ENS'),]

plot(x=rt.table$logCPM,y=rt.table$logFC,xlab='Averaged logCPM',ylab='logFC',cex.lab=1.4,cex.axis=1.5,ylim=c(-10,10),pch=19,font=2,cex=1.3,col='grey')
flag <- rt.table$FDR<0.01 & rt.table$logFC>1 
points(rt.table[flag,'logCPM'],rt.table[flag,'logFC'],col='red',pch=19)
abline(h=0)


#GFP neg, high dose Dux over-expression
gene.data <- gene.read.count.matrix[,c('SRR4032348.bam','SRR4032349.bam','SRR4032352.bam','SRR4032353.bam')]
rt.data   <- rt.read.count.matrix[,c('SRR4032348.bam','SRR4032349.bam','SRR4032352.bam','SRR4032353.bam')]  
de.data   <- rbind(gene.data,rt.data)
flag      <- apply(de.data,1,function(x) sum(x>=10) >=4)
de.data   <- de.data[flag,]

edgeR.res <- DGEList(counts=round(de.data), 
                     group=factor(c('Wild.type','Wild.type','Dux.OE.GFP.neg','Dux.OE.GFP.neg'),levels = c('Wild.type','Dux.OE.GFP.neg'))
)
edgeR.res <- calcNormFactors(edgeR.res, method="RLE")

edgeR.res <- estimateCommonDisp(edgeR.res)
edgeR.res <- estimateTagwiseDisp(edgeR.res)
et        <- exactTest(edgeR.res)
de.res    <- topTags(et, n=Inf)$table
rt.table  <- de.res[!grepl(x = rownames(de.res),pattern = 'ENS'),]

plot(x=rt.table$logCPM,y=rt.table$logFC,xlab='Averaged logCPM',ylab='logFC',cex.lab=1.4,cex.axis=1.5,ylim=c(-10,10),pch=19,font=2,cex=1.3,col='grey')
flag <- rt.table$FDR<0.01 & rt.table$logFC>1 
points(rt.table[flag,'logCPM'],rt.table[flag,'logFC'],col='red',pch=19)
abline(h=0)

gene.table  <- de.res[grepl(x = rownames(de.res),pattern = 'ENS'),]
plot(x=gene.table$logCPM,y=gene.table$logFC,xlab='Averaged logCPM',ylab='logFC',cex.lab=1.4,cex.axis=1.5,ylim=c(-15,15),pch=19,font=2,cex=1.3)
flag <- gene.table$FDR<0.01 & gene.table$logFC>1 
points(gene.table[flag,'logCPM'],gene.table[flag,'logFC'],col='red',pch=19)
abline(h=0)

##############################################################
# Otx2 chip-seq data to show  MT2_Mm enrichment
###############################################################
load('RData//otx2.csv.chip.seq.results.RData')
df   <- MA.list[['OTX2.3']]
plot(x=df$A,y=df$M,xlab='A',ylab='M',cex.lab=2,cex.axis=2,ylim=c(-5,5),pch=19,font=2,cex=1.5,col='grey')
flag <- df$fdr<=0.01 & df$M>=1
points(df$A[flag],df$M[flag],col='red',pch=19,cex=1.5)
abline(h=0)

##############################################################
# Otx2 RNA-Seq data to show  Dux induction and MERVL-int induction
###############################################################
load('RData//otx2.csv.rna.seq.results.RData')
Dux.id                                <- 'ENSMUSG00000075046.6'
Otx2.id                               <- 'ENSMUSG00000021848.13'


gene.data <- gene.read.count.matrix[,c('SRR1201006.bam','SRR1200997.bam','SRR1200998.bam')]
rt.data   <- rt.read.count.matrix[,c('SRR1201006.bam','SRR1200997.bam','SRR1200998.bam')]
de.data   <- rbind(gene.data,rt.data)
flag      <- apply(de.data,1,function(x) sum(x>=10) >=2)
de.data   <- de.data[flag,]

edgeR.res <- DGEList(counts=round(de.data), 
                     group=factor(c('Otx2.OE','Wild.type','Wild.type'),levels = c('Wild.type','Otx2.OE'))
                     )
edgeR.res <- calcNormFactors(edgeR.res, method="RLE")

edgeR.res <- estimateCommonDisp(edgeR.res)
edgeR.res <- estimateTagwiseDisp(edgeR.res)
et        <- exactTest(edgeR.res)
de.res    <- topTags(et, n=Inf)$table

plotSmear(et,cex.lab=1.5,cex.axis=1.5,de.tags = rownames(de.res)[de.res$FDR<0.01],font=2,cex=1,col='grey')
points(de.res['MERVL-int','logCPM'],de.res['MERVL-int','logFC'],col='red',pch=19,cex=1.5)
points(de.res['MT2_Mm','logCPM'],de.res['MT2_Mm','logFC'],col='red',pch=19,cex=1.5)
points(de.res[Dux.id,'logCPM'],de.res[Dux.id,'logFC'],col='red',pch=19,cex=1.5)
points(de.res[Otx2.id,'logCPM'],de.res[Otx2.id,'logFC'],col='red',pch=19,cex=1.5)
abline(h=0)

##############################################################
# Arid3a chip-seq data to show  MT2_Mm enrichment
###############################################################
load('RData//Arid3a.csv.chip.seq.results.RData')
df   <- MA.list[['Arid3a.2']]
df   <- df[df$A >= 4,]
plot(x=df$A,y=df$M,xlab='A',ylab='M',cex.lab=2,cex.axis=2,ylim=c(-5,5),pch=19,font=2,cex=1.5,col='grey',xlim=c(3,16))
flag <- df$fdr<=0.01 & df$M>=1
points(df$A[flag],df$M[flag],col='red',pch=19,cex=1.5)
abline(h=0)


##############################################################
# Oct4 chip-seq data to show  MT2_Mm enrichment
###############################################################
load('RData//oct4.csv.chip.seq.results.RData')
df   <- MA.list[['Pou5f1']]
plot(x=df$A,y=df$M,xlab='A',ylab='M',cex.lab=2,cex.axis=2,ylim=c(-5,5),pch=19,font=2,cex=1.5,col='grey')
flag <- df$fdr<=0.01 & df$M>=1
points(df$A[flag],df$M[flag],col='red',pch=19,cex=1.5)
abline(h=0)



##############################################################
# Nelfa chip-seq data to show   enrichment on 2C gene promoters
###############################################################
load('RData//Nelfa.csv.chip.seq.results.RData')
df   <- MA.list[['Nelfa.1']]
plot(x=df$A,y=df$M,xlab='A',ylab='M',cex.lab=2,cex.axis=2,ylim=c(-5,5),pch=19,font=2,cex=1.5,col='grey')
flag <- df$fdr<=0.01 & df$M>=1
points(df$A[flag],df$M[flag],col='red',pch=19,cex=1.5)
abline(h=0)




######################Trash##############
# flag      <- apply(gene.data,1,function(x) sum(x>=10) >=3)
# gene.data <- gene.data[flag,]
# lib.size  <- calcNormFactors(gene.data, method='RLE')
# 
# normalized.gene.data <- foreach(j = 1:ncol(gene.data),.combine='cbind') %do% {
#     gene.data[,j] / lib.size[j]  
#   
# }
# colnames(normalized.gene.data) <- colnames(gene.data)
# 
# normalized.rt.data <- foreach(j = 1:ncol(rt.data),.combine='cbind') %do% {
#   rt.data[,j] / lib.size[j]  
#   
# }
# colnames(normalized.rt.data) <- colnames(rt.data)
# 
# par(mfrow=c(2,2))
# foreach(i = 1:nrow(exp.meta))  %do% {
#     ip.name    <- sprintf("%s.bam",exp.meta$ip[i] )    %>% as.character
#     input.name <- sprintf("%s.bam",exp.meta$input[i] ) %>% as.character
#     M          <- log2(normalized.rt.data[,ip.name] + 1) - log2(normalized.rt.data[,input.name] + 1)
#     A          <- (log2(normalized.rt.data[,ip.name] + 1) + log2(normalized.rt.data[,input.name] + 1))/2
#     plot(x=A,y=M,main=exp.meta$exp.name[i],xlab='A',ylab='M',cex.lab=1.5,cex.axis=1.5,ylim=c(-4,10),pch=19)
#     points(A['MERVL-int'],M['MERVL-int'],col='red',pch=19)
#     abline(h=0)
# }
# 
# par(mfrow=c(2,2))
# foreach(i = 1:nrow(exp.meta))  %do% {
#   ip.name    <- sprintf("%s.bam",exp.meta$ip[i] )    %>% as.character
#   input.name <- sprintf("%s.bam",exp.meta$input[i] ) %>% as.character
#   M          <- log2(normalized.gene.data[,ip.name] + 1) - log2(normalized.gene.data[,input.name] + 1)
#   A          <- (log2(normalized.gene.data[,ip.name] + 1) + log2(normalized.gene.data[,input.name] + 1))/2
#   plot(x=A,y=M,main=exp.meta$exp.name[i],xlab='A',ylab='M',cex.lab=1.5,cex.axis=1.5,ylim=c(-10,14),pch=19)
#   abline(h=0)
# }




# Dux.dox.neg.vs.no.dox.1 <- gene.MA.list[[1]] %>% arrange(desc(A))
# Dux.dox.neg.vs.no.dox.2 <- gene.MA.list[[2]] %>% arrange(desc(A))
# Dux.dox.pos.vs.no.dox.1 <- gene.MA.list[[3]] %>% arrange(desc(A))
# Dux.dox.pos.vs.no.dox.2 <- gene.MA.list[[4]] %>% arrange(desc(A))
# rownames(Dux.dox.neg.vs.no.dox.1) <- Dux.dox.neg.vs.no.dox.1$rt.id
# rownames(Dux.dox.neg.vs.no.dox.2) <- Dux.dox.neg.vs.no.dox.2$rt.id
# rownames(Dux.dox.pos.vs.no.dox.1) <- Dux.dox.pos.vs.no.dox.1$rt.id
# rownames(Dux.dox.pos.vs.no.dox.2) <- Dux.dox.pos.vs.no.dox.2$rt.id
# 
# #Draw MA plot
# par(mfcol=c(1,2))
# plot(x=Dux.dox.neg.vs.no.dox.1$A,y=Dux.dox.neg.vs.no.dox.1$M,xlab='A',ylab='M')
# points(Dux.dox.neg.vs.no.dox.1[Dux.id,'A'],Dux.dox.neg.vs.no.dox.1[Dux.id,'M'],col='red',pch=19)
# abline(h=0)
# plot(x=Dux.dox.neg.vs.no.dox.2$A,y=Dux.dox.neg.vs.no.dox.2$M,xlab='A',ylab='M')
# points(Dux.dox.neg.vs.no.dox.2[Dux.id,'A'],Dux.dox.neg.vs.no.dox.2[Dux.id,'M'],col='red',pch=19)
# abline(h=0)
# 
# plot(x=Dux.dox.pos.vs.no.dox.1$A,y=Dux.dox.pos.vs.no.dox.1$M,xlab='A',ylab='M')
# points(Dux.dox.pos.vs.no.dox.1[Dux.id,'A'],Dux.dox.pos.vs.no.dox.1[Dux.id,'M'],col='red',pch=19)
# abline(h=0)
# plot(x=Dux.dox.pos.vs.no.dox.2$A,y=Dux.dox.pos.vs.no.dox.2$M,xlab='A',ylab='M')
# points(Dux.dox.pos.vs.no.dox.2[Dux.id,'A'],Dux.dox.pos.vs.no.dox.2[Dux.id,'M'],col='red',pch=19)
# abline(h=0)
# 
# #Draw distribution of M
# top.cnt <- 4000
# par(mfcol=c(1,2))
# hist(Dux.dox.neg.vs.no.dox.1[1:top.cnt,'M'],col=rgb(0,0,0,0.5),breaks=250,xlab='M',main='')
# hist(Dux.dox.pos.vs.no.dox.1[1:top.cnt,'M'],col=rgb(0,1,0,0.5),breaks=250,add=TRUE)
# hist(Dux.dox.neg.vs.no.dox.2[1:top.cnt,'M'],col=rgb(0,0,0,0.5),breaks=250,xlab='M',main='')
# hist(Dux.dox.pos.vs.no.dox.2[1:top.cnt,'M'],col=rgb(0,1,0,0.5),breaks=250,add=TRUE)
# 
# 
# draw.df.1        <- rbind(Dux.dox.neg.vs.no.dox.1[1:top.cnt,],Dux.dox.pos.vs.no.dox.1[1:top.cnt,])
# draw.df.1$status <- c(rep('GFP.negative',top.cnt),rep('GFP.positive',top.cnt))
# draw.df.2        <- rbind(Dux.dox.neg.vs.no.dox.2[1:top.cnt,],Dux.dox.pos.vs.no.dox.2[1:top.cnt,])
# draw.df.2$status <- c(rep('GFP.negative',top.cnt),rep('GFP.positive',top.cnt))
# ggplot(draw.df.1) + geom_density(aes(x=M,color=status))
# ggplot(draw.df.2) + geom_density(aes(x=M,color=status))
# 
# 
# 
# Dux.dox.neg.vs.no.dox.1 <- rt.MA.list[[1]] %>% arrange(desc(A))
# Dux.dox.neg.vs.no.dox.2 <- rt.MA.list[[2]] %>% arrange(desc(A))
# Dux.dox.pos.vs.no.dox.1 <- rt.MA.list[[3]] %>% arrange(desc(A))
# Dux.dox.pos.vs.no.dox.2 <- rt.MA.list[[4]] %>% arrange(desc(A))
# rownames(Dux.dox.neg.vs.no.dox.1) <- Dux.dox.neg.vs.no.dox.1$rt.id
# rownames(Dux.dox.neg.vs.no.dox.2) <- Dux.dox.neg.vs.no.dox.2$rt.id
# rownames(Dux.dox.pos.vs.no.dox.1) <- Dux.dox.pos.vs.no.dox.1$rt.id
# rownames(Dux.dox.pos.vs.no.dox.2) <- Dux.dox.pos.vs.no.dox.2$rt.id
# 
# rt.id <- 'MERVL-int'
# par(mfcol=c(1,2))
# plot(x=Dux.dox.neg.vs.no.dox.1$A,y=Dux.dox.neg.vs.no.dox.1$M,xlab='A',ylab='M',pch=19)
# points(Dux.dox.neg.vs.no.dox.1[rt.id,'A'],Dux.dox.neg.vs.no.dox.1[rt.id,'M'],col='red',pch=19)
# abline(h=0)
# plot(x=Dux.dox.neg.vs.no.dox.2$A,y=Dux.dox.neg.vs.no.dox.2$M,xlab='A',ylab='M',pch=19)
# points(Dux.dox.neg.vs.no.dox.2[rt.id,'A'],Dux.dox.neg.vs.no.dox.2[rt.id,'M'],col='red',pch=19)
# abline(h=0)
# 
# plot(x=Dux.dox.pos.vs.no.dox.1$A,y=Dux.dox.pos.vs.no.dox.1$M,xlab='A',ylab='M',pch=19)
# points(Dux.dox.pos.vs.no.dox.1[rt.id,'A'],Dux.dox.pos.vs.no.dox.1[rt.id,'M'],col='red',pch=19)
# abline(h=0)
# plot(x=Dux.dox.pos.vs.no.dox.2$A,y=Dux.dox.pos.vs.no.dox.2$M,xlab='A',ylab='M',pch=19)
# points(Dux.dox.pos.vs.no.dox.2[rt.id,'A'],Dux.dox.pos.vs.no.dox.2[rt.id,'M'],col='red',pch=19)
# abline(h=0)
# 
# top.cnt <- 200
# par(mfcol=c(1,2))
# hist(Dux.dox.neg.vs.no.dox.1[1:top.cnt,'M'],col=rgb(0,0,0,0.5),breaks=50,xlab='M',main='',xlim=c(-2,10))
# hist(Dux.dox.pos.vs.no.dox.1[1:top.cnt,'M'],col=rgb(0,1,0,0.5),breaks=50,add=TRUE)
# hist(Dux.dox.neg.vs.no.dox.2[1:top.cnt,'M'],col=rgb(0,0,0,0.5),breaks=50,xlab='M',main='',xlim=c(-2,10))
# hist(Dux.dox.pos.vs.no.dox.2[1:top.cnt,'M'],col=rgb(0,1,0,0.5),breaks=50,add=TRUE)
# 
# draw.df.1        <- rbind(Dux.dox.neg.vs.no.dox.1[1:top.cnt,],Dux.dox.pos.vs.no.dox.1[1:top.cnt,])
# draw.df.1$status <- c(rep('GFP.negative',top.cnt),rep('GFP.positive',top.cnt))
# draw.df.2        <- rbind(Dux.dox.neg.vs.no.dox.2[1:top.cnt,],Dux.dox.pos.vs.no.dox.2[1:top.cnt,])
# draw.df.2$status <- c(rep('GFP.negative',top.cnt),rep('GFP.positive',top.cnt))
# ggplot(draw.df.1) + geom_density(aes(x=M,color=status))
# ggplot(draw.df.2) + geom_density(aes(x=M,color=status))
