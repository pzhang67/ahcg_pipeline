#-------------
# open package
#-------------- 
library(ggplot2)
#------------
# read input
#-----------
args <- commandArgs(TRUE)
dt <- read.table(args[1],colClasses= c("character", rep("integer", 8)) ,header = TRUE,row.names=NULL)
m<-dim(dt)[2]
#----------
# function
#---------
plots.function = function (dt) 
{
##------------------------------------------
## adjust columns names and read genes names
##------------------------------------------
  genes.name<-(colnames(dt)[5:m])
  n<-length(genes.name)
  colnames(dt)<-c("CHR","POSITION","DEPTH",genes.name)
  dt<-dt[,1:9]
  dt.gene<-dt[,1:3]
  exons.list<-NA
  
  for (k in 1:n){ #for each gene
    dt.gene$EXON<-unlist(dt[,k+3])
    dt.gene$EXON<-as.numeric(dt.gene$EXON)
    chr.coo.null<-which(dt.gene$EXON!=0)
##------------
## gene cutoff 
##------------   
    quantile<-quantile(dt.gene$DEPTH[chr.coo.null], probs=0.1)
    median<-quantile(dt.gene$DEPTH[chr.coo.null], probs=0.5)
##--------------
## exons cutoff 
##--------------    
    exons.cutoff = function(dtt){
      n.exons<-unique(dtt$EXON)
      max.exons<-max(n.exons)
      min.exons<-min(n.exons)
      cutoff.exons.min<-list()
      cutoff.exons.max<-list()
      for (i in min.exons:max.exons){
        index<-which(dtt$EXON == i)
        data<-dtt$DEPTH[index]
        cutoff.exons.min<-append(cutoff.exons.min,boxplot.stats(data, coef = 1.5, do.conf = TRUE, do.out = TRUE)$stats[1])
        cutoff.exons.max<-append(cutoff.exons.max,boxplot.stats(data, coef = 1.5, do.conf = TRUE, do.out = TRUE)$stats[5])
      }
      out<-list(cutoff.max=cutoff.exons.max,max=max.exons,min=min.exons,cutoff=cutoff.exons.min)
      return(out)
    }
    cut.off.exons<-exons.cutoff(dt.gene[chr.coo.null,])
##----------------
## boxplots
##----------------
    maximum.y<-max(unlist(cut.off.exons$cutoff.max)) +20
    p<-ggplot(data = dt.gene[chr.coo.null,], aes(x= factor(EXON), y= DEPTH)) + geom_boxplot(outlier.colour = "green", outlier.size = 1 ,aes(group = cut_width(EXON, 0.25)))
    p+ ggtitle("Boxplots of Depth by Exons") + xlab("Exons") + ylab("Depth") + scale_y_continuous(breaks=seq(0,maximum.y,20))
    ggsave(filename=paste("SNP_report/boxplots/boxplot.exons",genes.name[k],"jpg",sep="."), width = 6.7, height = 6.7)
##-----------------------------------
## gene cutoff and exons cutoff table
##---------------------------------- 
    n.exons<-unique(as.factor(dt.gene$EXON[chr.coo.null]))
    n.exons<-as.list(n.exons)
    stat<-matrix(c(median,quantile, cut.off.exons$cutoff),ncol=2+(cut.off.exons$max-cut.off.exons$min +1),byrow=TRUE)
    rownames(stat)<-genes.name[k]
    colnames(stat)<-c("median","90th_quantile",c(cut.off.exons$min:cut.off.exons$max))
    write.csv(stat,file=paste("SNP_report/cutoffs/cutoff",genes.name[k],"csv",sep="."))
    
    list.cutoffs<-c(genes.name[k],rep("NULL",cut.off.exons$min-1))
    list.cutoffs<-append(list.cutoffs,cut.off.exons$cutoff)
    exons.list[k]<-list(list.cutoffs)
    if(k==1){
      max.n.exons<-cut.off.exons$max
    }else{
      max.n.exons<-max(max.n.exons,cut.off.exons$max)
    }
    rm(list.cutoffs)
  }
##--------------
## cutoffs table
##---------------
 max<-max(as.vector(max.n.exons))
 for (k in 1:n){
   length(exons.list[[k]])<-max
 }
 table.exons<-do.call("rbind", exons.list) 
 exons.table<-matrix(table.exons,nrow=n,ncol=max)
 write.table(table.exons,file="scripts/cutoff/P1cutoff.final.csv",row.names = FALSE,col.names = FALSE,sep=",")
}
answer<-plots.function(dt)




