setwd(readClipboard())
rm(list = ls())

library(psych)
library(reshape2)
library(ggplot2)

env <- read.table("nutrient.csv",header = T,row.names=1,sep=",")
spe <- read.table("module.csv",header = T,row.names=1,sep=",")

pearson <- corr.test(env, spe, method = 'pearson', adjust = 'none')
r <- data.frame(pearson$r)  
p <- data.frame(pearson$p)  

r$env <- rownames(r)
p$env <- rownames(p)
r <- melt(r, id = 'env')
p <- melt(p, id = 'env')
pearson <- cbind(r, p$value)
colnames(pearson) <- c('env', 'spe', 'pearson_correlation', 'p.value')
pearson$spe <- factor(pearson$spe, levels = colnames(spe))
head(pearson)  


library(dplyr)
pearson$env <- factor(pearson$env, levels = unique(pearson$env))
pearson$spe <- factor(pearson$spe, levels = unique(pearson$spe))
 


p1 <- ggplot() +
  geom_tile(data = pearson, aes(x = env, y = spe, fill = pearson_correlation)) +
  scale_fill_gradientn(colors = c('#4DBBD5FF', 'white', '#E64B35FF'), limit = c(-1, 1)) 

p1

pearson[pearson$p.value<0.05, 'sig'] <- '*'


head(pearson)  


p2 <- p1 +
  geom_text(data = pearson, aes(x = env, y = spe, label = sig), size = 4)

p2



p3 <- ggplot(pearson, aes(x = env, y = spe, fill = pearson_correlation)) +
  geom_tile() +
  scale_fill_gradientn(colors = c('#4DBBD5FF', 'white', '#E64B35FF'), limit = c(-1, 1))  +
  coord_polar() +
  theme_bw()+
  theme(panel.border = element_blank(),  
        panel.grid = element_blank(),   
        axis.title = element_blank()) + 
  geom_text(aes(label = sig, angle = 0, y = 1.1), size = 4,
            position = position_stack(vjust = 0.5)) +  
  scale_y_discrete(expand=expansion(mult=c(1,0)))+
  scale_x_discrete(expand=expansion(mult=c(0.4,0)))

p3

