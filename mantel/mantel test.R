rm(list=ls())

#
library(linkET)
library(dplyr)
library(ggplot2)

#
Otu <- as.data.frame(t(read.csv("otu_filter.csv",row.names = 1)))
Env <- read.csv("env.csv",row.names = 1)

Alp <- read.csv("alpha_diversity.csv",row.names = 1)
Otu[,ncol(Otu)+1] <- Alp$Shannon
Otu[,ncol(Otu)+1] <- Alp$Simpson
Otu[,ncol(Otu)+1] <- Alp$Chao1
Otu[,ncol(Otu)+1] <- Alp$ACE

a <- 1
#
# for (a in 1:2) {
  if(a==1){
    otu <- Otu[1:31,]
    env <- Env[1:31,]
  }else{
    otu <- Otu[32:68,]
    env <- Env[32:68,]
  }
# }

# Delete all 0 columns
otu <- otu[,colSums(otu != 0) > 0]

env <- as.data.frame(scale(env))
# otu <- as.data.frame(scale(otu))
if(a==1){
env <- env[,c("avde","pH","ORP","NO3.N","NO2.N","DON","Fe2","TP","T.Fe","T.Mn",
              "Cl.","SO42.","Na.","Ca2.","Mg2.","C.N","Soiltype")]
}else{
  env <- env[,c("avde","pH","ORP","moisture","NH4.N","NO3.N","NO2.N","DON","Fe2",
                "T.Mn","CEC", "Cl.","Na.","HCO32.","C.N","Soiltype")]
}

# calculation
mantel <- mantel_test(otu, env, method="pearson",
                      spec_select = list(OTU_abundance = 1:c(ncol(otu)-4),
                                         Alpha_diversity = c(ncol(otu)-3):ncol(otu)))
# p_values <- mantel$p
# 
# # Correction of p-values using the Benjamini-Hochberg method
# p_adjusted <- p.adjust(p_values, method = "BH")
# 
# # Add the corrected p-value back into the result
# mantel$p <- p_adjusted

mantel <- mutate(mantel,r_abs_d = cut(abs(r), breaks = c(-Inf, 0.2, 0.4, Inf), # 对相关系数进行分割，便于映射大小
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(abs(p), breaks = c(-Inf, 0.01, 0.05, Inf), # 对P值进行分割，便于映射颜色
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
         rd = ifelse(r > 0, ">0", ifelse(r == 0, "=0", "<0")))
## `mantel_test()` using 'bray' dist method for 'spec'.
## `mantel_test()` using 'euclidean' dist method for 'env'.

# 
qcorrplot(correlate(env, method="spearman"), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = r_abs_d,linetype= rd), # 这行代码是关键
              data = mantel, 
              curvature = nice_curvature()) +
  geom_mark(# 添加r值与显著性标记
    only_mark = TRUE,
    sep = '\n', 
    size = 3, 
    sig_level = c(0.05, 0.01, 0.001), # Significance level setting
    sig_thres = 0.05 # significance threshold, correlation coefficients with p-values greater than the threshold will not be plotted.
  ) +
  
  # 
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),limits=c(-1,1)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("< 0.01"="#B2192C","0.01 - 0.05"="#377E47", ">= 0.05"="grey")) +
  scale_linetype_manual(values = c("dashed","solid"))+
  guides(size = guide_legend(title = "abs of Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         linetype = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
