
rm(list = ls())
setwd(readClipboard())
getwd()
#读表
tax2 <- read.table("TAX2.csv",header = T,row.names=1,sep=",")
#表转置
trans_tax2 <- t(tax2)
#安装，这个包只能通过github安装
devtools::install_github("Hy4m/linkET", force = TRUE)
library(linkET)
library(dplyr)
library(ggplot2)

varechem <- read.table("nutrient.csv",header = T,row.names=1,sep=",")
varespec <- trans_tax2
mantel <- mantel_test(varespec, varechem,
                      spec_select = list(Spec01 = 34990:34994,
                                         Spec02 = 34995:35337, 
                                         Spec03 = 35338:35416,
                                         Spec04 = 35417:37734,
                                         Spec05 = 1:34929,
                                         Spec06 = 34930:34989)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), # 对相关系数进行分割，便于映射大小
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), # 对P值进行分割，便于映射颜色
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

## `mantel_test()` using 'bray' dist method for 'spec'.
## `mantel_test()` using 'euclidean' dist method for 'env'.


#自定义颜色
color1 <- "#4DBBD5FF"
color2 <- "#FFFFFF"
color3 <- "#E64B35FF"



# 画图
qcorrplot(correlate(varechem), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), # 这行代码是关键
              data = mantel, 
              curvature = nice_curvature()) +
    scale_fill_gradientn( colours = c(color1, color2, color3),
  values = scales::rescale(c(-1, 0, 1)), 
  limits = c(-1, 1)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
write.csv(mantel, 'mantel.csv',row.names =F)
