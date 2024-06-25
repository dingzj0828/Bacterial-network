
rm(list=ls())

install.packages("reshape2")
install.packages("ggrepel")


library(reshape2)
library(igraph)
library(ggrepel)
library(dplyr)
alltables<-c('tax.csv.network.txt')

#netname='Optimal_group.txt.network.csv'
for(netname in  alltables){
  df1<-read.delim(netname,header = T)
  df2<-df1
  df2[,c(1,2)]<-df2[,c(2,1)]
  df3<-rbind(df1,df2)
  df2<-dcast(Source~ Target, value.var = 'correlation',data = df3)
  rownames(df2)<-df2$Source
  df2<-df2[,-1]
  df2<-df2[,rownames(df2)]
  adjacency_unweight<-df2
  adjacency_unweight<-ifelse(is.na(df2),0,1)
  write.table(cbind(id=rownames(adjacency_unweight),adjacency_unweight),paste0(netname,'.adjacency_unweight.txt'),row.names = F,sep = '\t',quote = F)
  #邻接矩阵 -> igraph 的邻接列表，获得非含权的无向网络
  igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
  igraph    #igraph 的邻接列表
  #计算节点度
  V(igraph)$degree <- igraph::degree(igraph)
  #模块划分，详情 ?cluster_fast_greedy，有多种模型
  set.seed(123)
  V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))
  #输出各节点（微生物 OTU）名称、节点度、及其所划分的模块的列表
  nodes_list <- data.frame(
    nodes_id = V(igraph)$name, 
    degree = V(igraph)$degree, 
    #   modularity = V(igraph)$modularity
    modularity = as.numeric(V(igraph)$modularity)
    
  )
  head(nodes_list)    #节点列表，包含节点名称、节点度、及其所划分的模块
  write.table(nodes_list, paste0(netname,'.nodes_list.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
  
  ##计算模块内连通度（Zi）和模块间连通度（Pi）
  source('3.zi_pi.r')
  
  #上述的邻接矩阵类型的网络文件
  adjacency_unweight <- read.delim(paste0(netname,'.adjacency_unweight.txt'), row.names = 1, sep = '\t', check.names = FALSE)
  
  #节点属性列表，包含节点所划分的模块
  nodes_list <- read.delim(paste0(netname,'.nodes_list.txt'), row.names = 1, sep = '\t', check.names = FALSE)
  
  #两个文件的节点顺序要一致
  nodes_list <- nodes_list[rownames(adjacency_unweight), ]
  
  #计算模块内连通度（Zi）和模块间连通度（Pi）
  #指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
  zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
  head(zi_pi)
  zi_pi$netname<-netname
  
  zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
  withinset<-2.5#2.5
  amongset<-0.62#0.62
  
  zi_pi[which(zi_pi$within_module_connectivities <= withinset & zi_pi$among_module_connectivities <= amongset),'type'] <- 'Peripherals'
  zi_pi[which(zi_pi$within_module_connectivities < withinset & zi_pi$among_module_connectivities > amongset),'type'] <- 'Connectors'
  zi_pi[which(zi_pi$within_module_connectivities > withinset & zi_pi$among_module_connectivities < amongset),'type'] <- 'Module hubs'
  zi_pi[which(zi_pi$within_module_connectivities > withinset & zi_pi$among_module_connectivities > amongset),'type'] <- 'Network hubs'
  
  
  write.table(zi_pi, paste0(netname,'.zi_pi_result.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
  
  
  df4<-as.data.frame(zi_pi %>%dplyr::filter(type!='Peripherals'))
  

  
  
}



###

options(stringsAsFactors = F)
rm(list=ls())
library(ggplot2)
library(ggrepel)
library(dplyr)

alltables<-c('tax.csv.network.txt')
alltables<-paste(alltables,'.zi_pi_result.txt',sep = '')
#list.files(pattern = 'zi_pi_result.txt$')
df3<-NULL
for(i in 1:length(alltables)){
  if(i==1){df1<-read.delim(alltables[1],header = T)}else{
    df2<-read.delim(alltables[i],header = T)
    df1<-rbind(df1,df2)
  }
  
}
df3<-df1


df3$netname<-gsub("\\..*",'',df3$netname)
#下边这行代码是我更改过的，能运行出来，原始给出的是下下边标注的，运行不出来
df3$netname<-factor(df3$netname,levels = unique(df3$netname))
#下边这行代码是原始给出的是下下边标注的，运行不出来
df3$netname<-factor(df3$netname,levels = unique())

withinset<-2.5#2.5
amongset<-0.62#0.62

df4<-as.data.frame(df3 %>%filter(type!='Peripherals'))


mycol<-c("#F8766D") 
ggplot(df3, aes(among_module_connectivities, within_module_connectivities,group=netname)) +
  geom_point( aes(color=netname,shape=netname),alpha = 0.5, size = 2) +
  scale_color_manual(values=mycol)+
  #labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = amongset,lty=2) +
  geom_hline(yintercept =withinset,lty=2)+
  #geom_text(mapping=aes(among_module_connectivities, within_module_connectivities,label=nodes_id,color=netname),data=df4,
                 #size=3,show.legend = F)+
  annotate("text", label = "bold(Peripherals)", x = 0.2, y = -3, size = 4, colour = "black", parse = TRUE)+
  annotate("text", label = "bold(module_hubs)", x = 0.2, y = 5.5, size = 4, colour = "black", parse = TRUE)+
  annotate("text", label = "bold(network_hubs)", x = 0.8, y = 5.5, size = 4, colour = "black", parse = TRUE)+
  annotate("text", label = "bold(Connector)", x = 0.8, y = -3, size = 4, colour = "black", parse = TRUE)+
 # labs(x='Pi',y='Zi')+
  labs(x = 'Among-module connectivities(Pi)', y = 'Within-module connectivities(Zi)') +
  scale_x_continuous(limits = c(0.0,1))+
  scale_y_continuous(limits=c(-3.5,6))+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(face='bold'),
        axis.title = element_text(face='bold'),
        text = element_text(face='bold'))

ggsave('step43.all_zp.pdf',width = 7,height = 6)

