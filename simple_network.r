library(picante)
library(reshape2)
library(dplyr)
rm(list=ls())


#导入表格
otu<-read.csv("tax.csv", row.names=1)
otu
otu1<-otu
otu1[otu1>0]<-1
otu<-otu[which(rowSums(otu1)>=6),]
otu<-t(otu)


#做相关
result<-cor.table(otu,cor.method  = 'spearman')# spearman, pearson

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
corr_R_P<-as.data.frame(corr_R_P %>% filter(abs(correlation)>0.6,pvalue<0.05))
corr_R_P$Weight<-abs(corr_R_P$correlation)
corr_R_P$po_ne<-ifelse(corr_R_P$correlation>0,'positive','negative')
corr_R_P$color<-ifelse(corr_R_P$correlation>0,'#E41A1C','#377EB8')


#导出gephi用的文件
write.table(corr_R_P,'network_gephi.csv',row.names = F,quote = F,sep = ',')


#导出gephi中node需要的数据，同时添加属性。边的数据已经在corr_R_P里了
node<-as.data.frame(corr_R_P)
taxonomy<-read.csv('species.csv',row.name=1)
taxonomy<-taxonomy[as.character(corr_R_P$Source),]
node$phylum<-taxonomy$phylum
node$class<-taxonomy$class
node$order<-taxonomy$order
node$family<-taxonomy$family
node$genus<-taxonomy$genus
node$species<-taxonomy$species

#重新做了一张表，删除原来表格中不需要的列
node1<-node
node1 <- subset(node1, select = -(3:8))

write.csv(node1, 'network.node.csv',row.names =F)








