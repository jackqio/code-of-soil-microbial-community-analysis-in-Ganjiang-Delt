rm(list=ls())

#### A. Microeco ####

#Important package
library(microeco)
library(magrittr)
library(ggplot2)
library(ape)
theme_set(theme_bw())


##1.Import data##

otutable <- read.table(file = 'otu table.txt', header = TRUE)
taxtable <- read.table(file = 'tax table.txt', fill = TRUE, na.strings = "", header = TRUE, row.names = 1)#fill=true是填充空白值
taxtable %<>% tidy_taxonomy
sampletable <- read.table(file = 'sample table.txt', fill = TRUE, na.strings = "", header = TRUE, row.names = 1)
tree <- read.tree("otus.tree")

# Create a microtable object with more information
dataset <- microtable$new(sample_table = sampletable, otu_table = otutable, tax_table = taxtable, phylo_tree=tree)
dataset


##2.data processing##

#Modifying the order of grouping
dataset$sample_table$Group3 %<>% factor(., levels = c("X1",	"X2",	"Y2",	"Y4",	"Y3",	"Z7",	"Z11",	"Z8",	"Z2",	"Z9",	"Z4",	"Z10"))
#dataset$sample_table$Group4 %<>% factor(., levels = c("X1",	"X2",	"Y2",	"Y4",	"Y3",	"Z1",	"Z2",	"Z3",	"Z4",	"Z5",	"Z6",	"Z7"))
str(dataset$sample_table)

# use R subset function to filter taxa in tax_table
# Deleting of OTUs outside of archaea and bacteria
dataset$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
# another way with grepl function
dataset$tax_table %<>% .[grepl("Bacteria|Archaea", .$Kingdom), ]

# Deleting of OTUs for mitochondria and chloroplasts
dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))

# Data pruning ensures that the OTU table and the Samples table are consistent.
dataset$tidy_dataset()

# Maximum and minimum values of the number of sample sequences
dataset$sample_sums() %>% range

# Rarefying Samples Sequence 
dataset$rarefy_samples(sample.size = 63336)

# To ensure that unidentified taxa are involved in the analysis, adjustments are made to the taxa table
library(microeco)
library(magrittr)
# replace g__ with g__unclassified
dataset$tax_table$Genus[grepl("__$", dataset$tax_table$Genus)] %<>% paste0(., "unclassified")
# use default parameters, rel = T means relative abundance
dataset$cal_abund(rel=T)
# replace p__ with p__unclassified
dataset$tax_table$Phylum[grepl("__$", dataset$tax_table$Phylum)] %<>% paste0(., "unclassified")
# use default parameters
dataset$cal_abund(rel=T)
# save the abundance
dataset$save_abund(dirpath = "taxa abund")

# Calculating the alpha diversity
# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = F)
# save dataset$alpha_diversity to a directory
dataset$save_alphadiv(dirpath = "alpha_diversity")

# Calculating the beta diversity. If no method argument is provided, the function automatically computes the Bray-curtis, Jaccard, weighted Unifrac, and unweighted Unifrac matrices
# unifrac = FALSE means do not calculate unifrac metric
# require GUniFrac package installed
dataset$cal_betadiv(unifrac = FALSE)
# save dataset$beta_diversity to a directory
dataset$save_betadiv(dirpath = "beta_diversity")

# Extract a subset of samples CY
# remember first clone the whole dataset
# see https://chiliubio.github.io/microeco_tutorial/notes.html#clone-function
group_CY <- clone(dataset)
# select 'CY'
group_CY$sample_table <- subset(group_CY$sample_table, Group1 == "CY")
# or: group_CW$sample_table <- subset(group_CW$sample_table, grepl("CY", Group))
# use tidy_dataset to trim all the basic files
group_CY$tidy_dataset()
group_CY$cal_abund(rel=TRUE)
# #calculating the alpha diversity, but it not nessisary, because the results are same between Universal Set and subset
# group_CY$cal_alphadiv(PD = F)
# group_CY$save_alphadiv(dirpath = "alpha_diversity")

# Extract a subset of samples HH
group_HH <- clone(dataset)
# select 'HH'
group_HH$sample_table <- subset(group_HH$sample_table, Group1 == "HH")
# or: group_CW$sample_table <- subset(group_CW$sample_table, grepl("HH", Group))
# use tidy_dataset to trim all the basic files
group_HH$tidy_dataset()
group_HH$cal_abund(rel=TRUE)
# group_HH$cal_alphadiv(PD = F)
# group_HH$save_alphadiv(dirpath = "alpha_diversity")

# Filtering tree for beta NTI calcultaion 
dataset_filter <- clone(dataset)
dataset_filter
# In this example, mean relative abundance threshold 0.0001
dataset_filter$filter_taxa(rel_abund = 0.0001)
# emport the abundance and tree file for beta NTI calculation
write.tree(dataset_filter$phylo_tree,"Poyang202208_OTU_0.001per.tree")
write.csv(dataset_filter$otu_table,"Poyang202208_OTU_0.001per.csv")

#### B. Structural information####
library(microeco)
library(ggalluvial)
library(ggh4x)


## 2. Diversity ##
library(FSA)

## a.Alpha diversity ##
dtalpha <- dataset
t1 <- trans_alpha$new(dataset = dtalpha, group = "Group1")
# The transformed diversity data is stored in object$data_alpha ...
# The group statistics are stored in object$data_stat ...

library(vegan)
## Statistical Significance of differences
d_alpha <- data.frame(shannon = subset(t1$data_alpha, Measure == "Shannon")$Value, chao = subset(t1$data_alpha, Measure == "Chao1")$Value, group = subset(t1$data_alpha, Measure == "Chao1")$Group1, row.names = subset(t1$data_alpha,Measure == "Shannon")$Sample)
#shannon
# Calculating the Euclidean distance（baced Shannon index）
distance_matrix <- dist(d_alpha$shannon, method = "euclidean")
# PERMANOVA test of Shannon using adonis2 
adonis2(distance_matrix ~ d_alpha$group, data = d_alpha, permutations = 999)#用矩阵算，一般是多向量
wilcox.test(d_alpha$shannon ~ d_alpha$group, data = d_alpha, permutations = 999)#用向量算，看的是中位数，非正态分布

# Calculating the Euclidean distance（baced Chao1 index）
distance_matrix <- dist(d_alpha$chao, method = "euclidean")
# PERMANOVA test of Chao1 using adonis2 
adonis2(distance_matrix ~ d_alpha$group, data = d_alpha, permutations = 999)#用矩阵算，一般是多向量
wilcox.test(d_alpha$chao ~ d_alpha$group, data = d_alpha, permutations = 999)#用向量算，看的是中位数，非正态分布


## b.Beta diverstiy ##

dtbeta <- dataset

# create an trans_beta object
# measure parameter must be one of names(dataset$beta_diversity)
t1 <- trans_beta$new(dataset = dtbeta, measure = "bray")

# PCoA
t1$cal_ordination(ordination = "PCoA",scale_species = TRUE)
# t1$res_ordination is the ordination result list

##Gradient color fill according to average depth
library(ggplot2)
library(ggbiplot)
library(ggnewscale)
library(dplyr)
# Gradient Fill Scatter
rm(pcaxy,pcaaro)
pcaxy <- as.data.frame(t1$res_ordination$scores[,1:2])
pcaxy$avde <- t1$res_ordination$scores$avde
pcaxy$Group1 <- t1$res_ordination$scores$Group1
colnames(pcaxy)[1] <- "PC1"
colnames(pcaxy)[2] <- "PC2"
pcaaro <- data.frame(pc1=t1$res_ordination$loading[1:10,1])
pcaaro$pc2 <- t1$res_ordination$loading[1:10,2]
pcaaro$lab <- t1$res_ordination$loading$Genus[1:10]
# Calculate the mean and standard deviation for each subgroup
group_summary <- pcaxy %>%
  group_by(Group1) %>%
  summarize(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2),
            sd_PC1 = sd(PC1), sd_PC2 = sd(PC2))
# Calculation of PC1 and PC2 explanatory degrees
exp_adj <- t1$res_ordination$eig/sum(t1$res_ordination$eig)  #获取未校正后的 R2
cap1_exp <- paste('CAP1:', round(exp_adj[1]*100, 2), '%')
cap2_exp <- paste('CAP2:', round(exp_adj[2]*100, 2), '%')

pcaxy$depth_level <- sampletable$Depthlevel1[1:68]


# ##Different shap for different average depth
# p <- ggplot() +
#   geom_point(data=subset(pcaxy,Group1=='CY'), 
#              mapping=aes(PC1, PC2, color="#f66f69",shape = depth_level),
#              size = 3) + 
#   geom_point(data=subset(pcaxy,Group1=='HH'),
#              mapping=aes(PC1, PC2, color="#1597a5",shape = depth_level),
#              size = 3
#   ) + 
# 
#   stat_ellipse(data=pcaxy,
#                geom = "polygon",level=0.95,
#                linetype = 1,size=0.8,
#                aes(x = PC1, y = PC2,fill=Group1),
#                alpha=0.1, show.legend = F)+
#   scale_shape_manual('depth_level', values = c(4, 6, 8, 10, 12, 24, 25))+
#   theme_bw() 
# p

# Gradient size
p <- ggplot() +
  geom_point(data=subset(pcaxy,Group1=='CY'), 
             mapping=aes(PC1, PC2, color=subset(pcaxy,Group1=='CY')[['avde']], 
                         size = subset(pcaxy,Group1=='CY')[['avde']]), 
             shape = 16) +
  scale_color_gradient('CY', low = "#ffbfb8", high = "#f66f69") +
  new_scale("color") +
  geom_point(data=subset(pcaxy,Group1=='HH'),
             mapping=aes(PC1, PC2, color=subset(pcaxy,Group1=='HH')[['avde']], 
                         size = subset(pcaxy,Group1=='HH')[['avde']]), 
             shape = 16) +
  scale_color_gradient('HH', low = "#CCFAFF", high = "#1597a5") + 
  theme_bw() +
  stat_ellipse(data=pcaxy,
               geom = "polygon", level=0.95,
               linetype = 1, size=0.8,
               aes(x = PC1, y = PC2, fill=Group1),
               alpha=0.1, show.legend = F) +
  scale_size_continuous(name="avde", range = c(2, 6))  # 设置大小范围

p

# Permutational multivariate analysis of variance，PERMANOVA
t1 <- trans_beta$new(dataset = dtbeta, group = "Group1", measure = "bray")
# manova for all groups when manova_all = TRUE
t1$cal_manova(manova_all = TRUE)
t1$res_manova
# The result is stored in object$res_manova ...
# manova for each paired groups
t1$cal_manova(manova_all = FALSE)
t1$res_manova



####C. Lefse#####
library(microeco)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggsignif)

dtmode <- dataset

t1 <- trans_diff$new(dataset = dtmode, 
                     method = "lefse", 
                     alpha = 0.01, 
                     taxa_level="Genus",
                     group = "Group1", 
                     lefse_subgroup = NULL)
# see t1$res_diff for the result
p1 <- t1$plot_diff_bar(color_values = c("#f66f69","#1597a5"),use_number = 1:20)#+theme(text = element_text(size= 14, family = "A"))#group_order = c("X", "Y", "Z")可以用来调整顺序
p1

# #展示Lefse结果中物种的丰度
p1 <- t1$plot_diff_abund(color_values = c("#f66f69","#1597a5"),use_number = 1:20,add_sig = TRUE)#+theme(text = element_text(size= 14, family = "A"))
p1
