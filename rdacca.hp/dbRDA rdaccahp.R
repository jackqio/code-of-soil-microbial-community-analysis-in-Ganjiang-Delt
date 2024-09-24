###db-RDA
rm(list=ls())

#two group
area <- c("CY","HH")
i <- 1 #1 or 2 for two group results

library(rdacca.hp)
library(vegan)

mite <- data.frame(t(read.csv("otu.csv",row.names = 1)))
mite[1:6,1:6]

avd <- read.csv(file = "avd.csv",row.names = 1)

mite.env <- read.csv("env.csv",row.names = 1)
mite.env <- as.data.frame(scale(mite.env))
mite.env <- as.data.frame(c(avd,mite.env),row.names = rownames(avd))
# Screening of indicators. Indicator screening is referenced to the calculation of the VIF
if(i==1){
  mite.env <- mite.env[,c("avde","pH","ORP","NO3.N","NO2.N","DON","Fe2","TP","T.Fe","T.Mn",
                          "Cl.","SO42.","Na.","Ca2.","Mg2.","C.N","Soiltype")]
}else{
  mite.env <- mite.env[,c("avde","pH","ORP","moisture","NH4.N","NO3.N","NO2.N","DON","Fe2",
                          "T.Mn","CEC", "Cl.","Na.","HCO32.","C.N","Soiltype")]
}

mite <-  mite[rownames(mite) %in% rownames(mite.env), ]

mite_CY <- mite[1:31,]
mite_HH <- mite[32:68,]
mite.env_CY <- mite.env[1:31,-1]
mite.env_HH <- mite.env[32:68,-1]
avd_CY <- avd[1:31,]
avd_HH <- avd[32:68,]


if(i==1){
o <- mite_CY
e <- mite.env_CY
d <- avd_CY
}else if(i==2){
  o <- mite_HH
  e <- mite.env_HH
  d <- avd_HH
}

####db-RDA ####
### Bray-curtis dissamillty ###

mite.bray <- vegdist(o, method = 'bray')
mite.bray0 <- as.data.frame(as.matrix(mite.bray))

#db-RDA 
mite.cap <- dbrda(mite.bray~., e)
summary(mite.cap)  #db-RDA summary

# Significance of the effect of each indicator on subgroups. 
# ANOVA analysis of rda results no longer needs to consider normal distribution
value_anova <- anova(mite.cap,by="term")

# Significance-adjusted values for the effect of each indicator on subgroups. 
# When there are more than 10 environmental factors, try to adjust the p-value to avoid false positives
value_anova_adj <- value_anova
value_anova_adj$`Pr(>F)` <- p.adjust(value_anova$`Pr(>F)`, method = 'BH')


## Significance of the correlation of whole environmental factors on community change
dbRDA.perm=permutest(mite.cap,permu=999)

## Examine the correlation between each environmental factor and community change
dbRDA.env <- envfit(mite.cap,e,permu=999, choices = c(1,2), display = 'sites')
#adjusted-p
dbRDA.env_adj <- dbRDA.env
dbRDA.env_adj$vectors$pvals <- p.adjust(dbRDA.env$vectors$pvals, method = 'bonferroni')




exp_adj <- RsquareAdj(mite.cap)$adj.r.squared * mite.cap$CCA$eig/sum(mite.cap$CCA$eig)  #Get the adjusted R2
# exp_adj <- mite.cap$CCA$eig/sum(mite.cap$CCA$eig)  #Get the un-adjusted R2

cap1_exp <- paste('dbRDA1:', round(exp_adj[1]*100, 2), '%')
cap2_exp <- paste('dbRDA2:', round(exp_adj[2]*100, 2), '%')

##
library(ggplot2)
library(ggbiplot)
library(ggnewscale)
library(dplyr)
# 
rm(pcaxy,pcaaro)
#
pcaxy <- as.data.frame(mite.cap$CCA$wa[,1:2])
pcaxy$avde <- d
pcaxy$Group1 <- rep(area[i],time=length(d))
# colnames(pcaxy)[1] <- "PC1"
# colnames(pcaxy)[2] <- "PC2"
#
pcaaro <- data.frame(pc1=mite.cap$CCA$biplot[,1],pc2=mite.cap$CCA$biplot[,2])
pcaaro$lab <- row.names(mite.cap$CCA$biplot)

if(i==1){
p <- ggplot() +
  geom_point(data=subset(pcaxy,Group1=='CY'), 
             mapping=aes(dbRDA1, dbRDA2, color=subset(pcaxy,Group1=='CY')[['avde']],
                         size = subset(pcaxy, Group1 == 'CY')[['avde']])) + 
  scale_color_gradient('CY', low = "#ffbfb8", high = "#c4323f")+
  geom_segment(data = pcaaro,
               aes(x = 0, xend = pc1*0.6, y = 0,
                   yend = pc2*0.6),
               # size = 2,
               arrow = arrow(length = unit(0.05, "npc")),
               alpha = 1)+
  geom_text(data = pcaaro,
            aes(x = pc1*0.6, y = pc2*0.6, label = lab),
            hjust = 1, vjust = 1,  # 调整标签的位置
            nudge_x = 0.05, nudge_y = 0.05,  # 进一步微调标签的位置
            size = 4) + # 设置标签的大小
  scale_size_continuous(name="avde", range = c(2, 6))+
  labs(x=cap1_exp,y=cap2_exp)+
  # stat_ellipse(data=pcaxy,
  #              geom = "polygon",level=0.95,
  #              linetype = 1,size=0.8,
  #              aes(x = dbRDA1, y = dbRDA2,fill=Group1),
  #              alpha=0.1, show.legend = F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")  # 添加0轴虚线# 添加0轴虚线
} else if(i==2){
p <- ggplot()+
  geom_point(data=subset(pcaxy,Group1=='HH'),
             mapping=aes(dbRDA1, dbRDA2, color=subset(pcaxy,Group1=='HH')[['avde']],
                         size = subset(pcaxy, Group1 == 'HH')[['avde']])
  ) + 
  scale_color_gradient('HH', low = "#00FFC4", high = "#0e606b")+theme_bw()+
  geom_segment(data = pcaaro,
               aes(x = 0, xend = pc1*0.6, y = 0,
                   yend = pc2*0.6),
               # size = 2,
               arrow = arrow(length = unit(0.05, "npc")),
               alpha = 1)+
  # geom_text(data = pcaaro,
  #           aes(x = pc1*0.6, y = pc2*0.6, label = lab),
  #           hjust = 1, vjust = 1,  # 调整标签的位置
  #           nudge_x = 0.05, nudge_y = 0.05,  # 进一步微调标签的位置
  #           size = 4) + # 设置标签的大小
  scale_size_continuous(name="avde", range = c(2, 6))+
  labs(x=cap1_exp,y=cap2_exp)+
# stat_ellipse(data=pcaxy,
#              geom = "polygon",level=0.95,
#              linetype = 1,size=0.8,
#              aes(x = dbRDA1, y = dbRDA2,fill=Group1),
#              alpha=0.1, show.legend = F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # 移除主要网格线
           panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")  # 添加0轴虚线# 添加0轴虚线
}
p


value_anova_adj#difference-in-difference test

dbRDA.perm#

dbRDA.env#Fitting environmental variables to PCA unconstrained axes with multiple regression

dbRDA.env_adj#Adjusted p-value

