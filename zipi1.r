#test
#read otu table
setwd("~/Desktop/R/410243571_robustness/test2/")
setwd(readClipboard())
getwd()
rm(list=ls())

install.packages("picante")
install.packages("dplyr")
install.packages("reshape2")
library(picante)
library(dplyr)
library(reshape2)

alltables<-c("tax.csv")
for(input1 in alltables){
 # input1<-"Asteraceae_otu_table.csv"
  print(input1)
  r_thresh<-0.6
  p_thresh<-0.05
  otutab<-read.delim(input1,header = T,row.names=1,sep=",")
  
  # mygroup<-read.delim(group1,header = T)
  # colnames(mygroup)[1]<-'sample'
  # otutab<-otutab[,mygroup$sample]
  
  otutab[is.na(otutab)]<-0
  ##keep 3 otus
  counts<-rowSums(otutab>0)
  
  otutab<-otutab[counts>=6,]
  
  sd1<-apply(otutab, MARGIN = 1, sd)
  otutab<-otutab[sd1>0,]
  comm<-t(otutab)
  colSums(otutab)
  
  sp.ra<-colMeans(comm)/mean(colSums(otutab))  #relative abundance of each species
  
  ###### there are two choices to get the correlation matrix #######
  ###### choice 1 (slow): calculate correlation matrix from OTU table
  
  result<-cor.table(comm,cor.method  = 'spearman')# spearman, pearson
  
  
  cormatrix=result$r
  
  #row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm) # if processed using MENAP, OTU order should match in the original OTU table and the correlation matrix downloaded from MENAP.
  
  cormatrix2<-cormatrix*(abs(cormatrix)>=r_thresh)  #only keep links above the cutoff point
  cormatrix2[is.na(cormatrix2)]<-0
  cormatrix2[result$P>0.05]<-0
  diag(cormatrix2)<-0    #no links for self-self    
  
  sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
  sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.
  
  network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
  #write.table(cbind(id=rownames(network.raw),network.raw),paste0(input1,'.network.raw.txt'),row.names = F,sep = '\t',quote = F)
  
  sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
  sum(row.names(network.raw)==names(sp.ra2))  #check if matched
  
  ## robustness simulation 
  #input network matrix, percentage of randomly removed species, and ra of all species
  #return the proportion of species remained
  
  rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
    id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
    net.Raw=netRaw #don't want change netRaw
    net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
    if (abundance.weighted){
      net.stength= net.Raw*sp.ra
    } else {
      net.stength= net.Raw
    }
    
    sp.meanInteration<-colMeans(net.stength)
    
    id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
    remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
    #for simplicity, I only consider the immediate effects of removing the
    #'id.rm' species; not consider the sequential effects of extinction of
    # the 'id.rm2' species.
    
    #you can write out the network pruned
    #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
    # write.csv( net.Raw,"network pruned.csv")
    
    remain.percent
  }
  
  rm.p.list=seq(0.05,0.2,by=0.05)
  rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=50){
    t(sapply(rm.p.list,function(x){
      remains=sapply(1:nperm,function(i){
        rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
      })
      remain.mean=mean(remains)
      remain.sd=sd(remains)
      remain.se=sd(remains)/(nperm^0.5)
      result<-c(remain.mean,remain.sd,remain.se)
      names(result)<-c("remain.mean","remain.sd","remain.se")
      result
    }))
  }
  
  
  #Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
  #Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)
  Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=rep(0.5,10), sp.ra=sp.ra2, abundance.weighted=T,
                        nperm=50)
  
  dat1<-data.frame(#Proportion.removed=rep(seq(0.05,1,by=0.05),1),
    Proportion.removed=0.5,
                   rbind(Weighted.simu),
                   weighted=rep(c("weighted"),each=10))
  
  currentdat<-dat1
  
  currentdat$network<-input1
  
  write.table(currentdat,paste0(input1,".robustness.txt"),
              row.names = F,sep = '\t',quote = F)
  library(reshape2)
  
  corr_r<-result$r
  corr_r[upper.tri(corr_r,diag = T)]<-NA
  corr_R<-melt(corr_r)
  colnames(corr_R)<-c('Source','Target','correlation')
  
  corr_p<-result$P
  corr_p[upper.tri(corr_p,diag = T)]<-NA
  corr_P<-melt(corr_p)
  corr_R_P<-corr_R
  corr_R_P$pvalue<-corr_P$value
  corr_R_P<-na.omit(corr_R_P)
  
  corr_R_P$fdr<-p.adjust(corr_R_P$pvalue,method = 'fdr')
  corr_R_P<-as.data.frame(corr_R_P %>% filter(abs(correlation)>=r_thresh,pvalue<0.05))
  corr_R_P$Weight<-abs(corr_R_P$correlation)
  corr_R_P$po_ne<-ifelse(corr_R_P$correlation>0,'positive','negative')
  corr_R_P$color<-ifelse(corr_R_P$correlation>0,'#E41A1C','#377EB8')
  write.table(corr_R_P,paste0(input1,'.network.txt'),row.names = F,sep = '\t',quote = F)
}
