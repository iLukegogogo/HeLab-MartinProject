#! /usr/bin/Rscript

require(plyr)
require(dplyr)
require(data.table)
require(foreach)


get.zscore <- function(m,a,n1,n2){ #n1 for IP, n2 for input
  condition.mean <- 0
  tmp            <- 2^(a-log2(sqrt(n1)) - log2(sqrt(n2)))    
  condition.var  <- (4*(1-tmp)* (1/log(2)) * (1/log(2))) / ((n1+n2)*tmp)
  z.score        <- (m-condition.mean) / sqrt(condition.var)
  z.score
}



options  <- commandArgs(trailingOnly = TRUE)
exp.file <- options[1]
exp.meta <- read.csv(file=exp.file,header=TRUE,stringsAsFactors=FALSE)





data           <- fread(input = 'RData/rt.read.count.matrix',header = TRUE,stringsAsFactors = FALSE) %>% as.data.frame
data$Chr       <- NULL
data$Start     <- NULL
data$Strand    <- NULL
data$End       <- NULL
data$Length    <- NULL
rownames(data) <- data$Geneid
data$Geneid    <- NULL
rt.read.count.matrix <- as.matrix(data)



data            <- read.table(file = 'RData/rt.read.count.matrix.summary',header = TRUE,stringsAsFactors = FALSE) %>% as.data.frame
data$Status     <- NULL
m               <- as.matrix(data)
lib.size        <- apply(m,2,sum)
names(lib.size) <- colnames(m)


MA.list <- foreach(i = 1:nrow(exp.meta)) %do% {
    row            <- exp.meta[i,]
    input          <- row$input
    ip             <- row$ip
    exp.name       <- row$exp.name
    input.bam.file <- paste(strsplit(x=input,split=":") %>% unlist,".all.rmdup.bam",sep="") %>% as.character
    ip.bam.file    <- paste(strsplit(x=ip,   split=":") %>% unlist,".all.rmdup.bam",sep="") %>% as.character
    
    input.rt.read.count.matrix <- rt.read.count.matrix[,input.bam.file]
    ip.rt.read.count.matrix    <- rt.read.count.matrix[,ip.bam.file]
    
    if(length(input.bam.file) == 1){
        input.rt.read.count <- input.rt.read.count.matrix
    }else{
        input.rt.read.count <- apply(input.rt.read.count.matrix,1,sum)
    }
    
    
    if(length(ip.bam.file) == 1){
        ip.rt.read.count <- ip.rt.read.count.matrix
    }else{
        ip.rt.read.count <- apply(ip.rt.read.count.matrix,1,sum)
    }
    
    n1 <-  sum(lib.size[ip.bam.file])
    n2 <-  sum(lib.size[input.bam.file])
    m  <-  (log2(ip.rt.read.count + 1) - log2(input.rt.read.count + 1)) - (log2(n1+1) - log2(n2+1))
    a  <-  (log2(ip.rt.read.count + 1) + log2(input.rt.read.count + 1))/2
    df         <- data.frame(M=m,A=a)
    loess.fit  <- loess(formula = M ~ A,data = df,span = 0.3)
    a          <- loess.fit$x %>% c
    m          <- loess.fit$residuals
    z.score    <-  get.zscore(m,a,n1,n2)
    p.value    <- pnorm(z.score,lower.tail = FALSE)
    fdr        <- p.adjust(p.value,method='fdr')
    data.frame(M=m,A=a,z.score=z.score,p.value=p.value,fdr=fdr,rt.id=names(m))
}
names(MA.list) <- as.character(exp.meta$exp.name) 


bam.file.list   <- system('ls alignment|grep all.rmdup.bam',wait=TRUE,intern=TRUE)
lib.read.length <- foreach(bam.file = bam.file.list,.combine='c') %do% {
    cmd <- sprintf("samtools view alignment/%s|head -n 1|cut -f 10",bam.file)
    seq <- system(cmd,wait=TRUE,intern=TRUE)
    l   <- ifelse(identical(seq, character(0)),0,nchar(seq))
    l
    
}
names(lib.read.length) <- bam.file.list


rs.file <- sprintf("RData/%s.chip.seq.results.RData",exp.file)
save(file=rs.file,list=c('MA.list','lib.size','rt.read.count.matrix','exp.meta','lib.read.length'))

