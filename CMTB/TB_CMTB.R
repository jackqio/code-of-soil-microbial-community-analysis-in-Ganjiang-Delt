library(vegan)
library(radiant.data)
library(FD)

rm(list=ls())

# setwd("~/Desktop/")
# input env data, standardize, calculate Eucldiean distance
## Rarefy
# Import
otu = read.csv('otu.csv', row.names=1)
# Sum
colSums(otu)
# Rarefy
otu_Flattening = as.data.frame(t(rrarefy(t(otu), min(colSums(otu)))))
# Sum rerafy
colSums(otu_Flattening)
# save
write.csv (otu_Flattening, file ="otu_Flattening.csv") #结果导出

env = read.csv("env.csv", header = T, row.names = 1)
env.std = decostand(env,method = "standardize")
env.dist = vegdist(env.std,"euclidean")

# input original otutable
otutable.16s.all = read.csv("otu_Flattening.csv", header = T, row.names = 1)
# otutable.16s.all = read.csv("otu.csv", header = T, row.names = 1)
otutable.16s.all = as.matrix(otutable.16s.all)

# use grouping file to align the sample, sample.id in grouping file should be the same as those in env file
grouping = read.csv("grouping.csv", header = T, row.names = 1)
otutable = t(t(otutable.16s.all)[match(rownames(grouping),rownames(t(otutable.16s.all))),]) 

# calculate total reads in otutable
sum.reads = sum(otutable[,1])

# use cutoff of abundance to select OTUs
otu.1 = ifelse(otutable < sum.reads*0.0001, 0, otutable)

# use occurrence (in >=3 samples) to further select OTUs
otu.1.sub = otu.1[rowSums(otu.1 != 0) >= 3, ]

sum(otu.1[,1])
sum(otu.1.sub[,1])
dim(otu.1)
dim(otu.1.sub)


# to calculate the TB of each OTUs, we apply Functional dispersion (FDis)
# in original FDis (Laliberté, E. and P. Legendre (2010) Ecology), we need two matrices:
# first is the traits (columns) of every species (rows), in order to calculate the functional distance (gowdis) between species
# second is the (relative) abundance of every species (columns) in every community (rows), as the weights
# the original FDis indicates the weighted mean distance of each species to the centroid for each community

# in our case of TB, for the 1st matrix, we can consider samples as species and ENV variables as different traits
# we standardize the ENVs (i.e., traits) and calculate the Euclidean distance
# for 2nd matrix, we consider OTUs as communities and samples as species
# thus, for each OTU (i.e., community), the #reads show the where this OTU detected in which samples (i.e., species) 
# the FDis here indicates the weighted mean distance of samples (where a OTU detected) to the centroid for each OTU
# larger FDis suggests that the ENVs (like traits) of samples are more distinct
# and thus the target OTU (like a community) can persis along a wider range of env and has large TB

FD = fdisp(env.dist, otu.1.sub)
TB = as.data.frame(FD$FDis)

write.csv(TB,"TB.abund.1.df.csv")
tb = read.csv("TB.abund.1.df.csv", row.names = 1)
names(tb)[1] = "TB"

###

# subset OTUs in tb
otu.tb = merge(tb, otu.1.sub, by = "row.names", all.x = T)
rownames(otu.tb) = otu.tb[,1]
otu.tb = otu.tb[,-1]
reads.present = colSums(otu.tb[,-1])/sum(otu.1[,1]) # we resampled 9736 reads
min(reads.present)
max(reads.present)
mean(reads.present)
sum(reads.present>0.9)

# OTU richness for tb.abund.1
richness = colSums(ifelse(otu.tb == 0, NA, 1)[ ,-1], na.rm = T)
min(richness)
max(richness)
mean(richness)


### CMTB calculation
## Method 1
## abundance-based version (weighted or unweighted)
# rank OTUs by abundance in each sample and calculate the community-level TB
sample.num = dim(otu.1.sub)[2]
tb.null = as.data.frame(matrix(nrow = 1, ncol = 13))
colnames(tb.null) = c("sample.id", "mean0.9", "mean0.85","mean0.8","mean0.75","cv0.9", "cv0.85","cv0.8","cv0.75",
                      "weighted.mean0.9","weighted.mean0.85","weighted.mean0.8","weighted.mean0.75")
rownames(tb.null) = "blank"

for (i in 1:sample.num) {
  df = otu.tb[, c(1,1+i)]
  sample.id = colnames(df[2])
  
  rank.tb = df[order(df[,1],decreasing = T),]
  rank.abund = rank.tb[order(rank.tb[,2],decreasing = T),]
  rank.abund$cumsum = cumsum(rank.abund[[2]])
  abund = sum(rank.abund[,2])/sum(otu.1[,1])
  
  # for different abundance cutoffs
  rank.abund$cutoff0.9 = ifelse(rank.abund[[3]] >= (sum(otu.1[,1])*0.9), NA, 1)
  rank.abund$cutoff0.85 = ifelse(rank.abund[[3]] >= (sum(otu.1[,1])*0.85), NA, 1)
  rank.abund$cutoff0.8 = ifelse(rank.abund[[3]] >= (sum(otu.1[,1])*0.8), NA, 1)
  rank.abund$cutoff0.75 = ifelse(rank.abund[[3]] >= (sum(otu.1[,1])*0.75), NA, 1)
  
  # mean TB and community weighted mean TB
  rank.abund$tb0.9 <- as.numeric(rank.abund$TB)*rank.abund$cutoff0.9
  rank.abund0.9 = rank.abund[ which(rank.abund$cutoff0.9 != "NA"),] 
  mean0.9 = mean(rank.abund0.9$tb0.9, na.rm = T)
  cv0.9 = sd(rank.abund0.9$tb0.9, na.rm = T)/mean0.9
  weighted.mean0.9 = weighted.mean(rank.abund0.9$tb0.9,rank.abund0.9[[2]])
  
  rank.abund$tb0.85 <- as.numeric(rank.abund$TB)*rank.abund$cutoff0.85
  rank.abund0.85 = rank.abund[ which(rank.abund$cutoff0.85 != "NA"),] 
  mean0.85 = mean(rank.abund0.85$tb0.85, na.rm = T)
  cv0.85 = sd(rank.abund0.85$tb0.85, na.rm = T)/mean0.85
  weighted.mean0.85 = weighted.mean(rank.abund0.85$tb0.85,rank.abund0.85[[2]])
  
  rank.abund$tb0.8 <- as.numeric(rank.abund$TB)*rank.abund$cutoff0.8
  rank.abund0.8 = rank.abund[ which(rank.abund$cutoff0.8 != "NA"),] 
  mean0.8 = mean(rank.abund0.8$tb0.8, na.rm = T)
  cv0.8 = sd(rank.abund0.8$tb0.8, na.rm = T)/mean0.8
  weighted.mean0.8 = weighted.mean(rank.abund0.8$tb0.8,rank.abund0.8[[2]])
  
  rank.abund$tb0.75 <- as.numeric(rank.abund$TB)*rank.abund$cutoff0.75
  rank.abund0.75 = rank.abund[ which(rank.abund$cutoff0.75 != "NA"),] 
  mean0.75 = mean(rank.abund0.75$tb0.75, na.rm = T)
  cv0.75 = sd(rank.abund0.75$tb0.75, na.rm = T)/mean0.75
  weighted.mean0.75 = weighted.mean(rank.abund0.75$tb0.75,rank.abund0.75[[2]])

  # output results
  tb =  as.data.frame(matrix(nrow = 1, ncol = 13))
  colnames(tb) = c("sample.id", "mean0.9", "mean0.85","mean0.8","mean0.75","cv0.9", "cv0.85","cv0.8","cv0.75",
                        "weighted.mean0.9","weighted.mean0.85","weighted.mean0.8","weighted.mean0.75")
  tb[] <- c(sample.id, mean0.9, mean0.85, mean0.8, mean0.75, cv0.9, cv0.85, cv0.8, cv0.75,
            weighted.mean0.9, weighted.mean0.85, weighted.mean0.8, weighted.mean0.75)
  
  tb.merge = rbind(tb.null, tb)
  tb.null = tb.merge
}

tb.merge = tb.merge[-1,]
tb.merge[,-1] = apply(tb.merge[,-1], 2, function(x) as.numeric(as.character(x)))
write.csv(tb.merge, "tb.merge rarefy.csv")

# correlation test
cor.test(as.numeric(tb.merge$mean0.9), as.numeric(tb.merge$weighted.mean0.9), method = "pearson")
cor.test(as.numeric(tb.merge$mean0.75), as.numeric(tb.merge$weighted.mean0.75), method = "pearson")
cor.test(as.numeric(tb.merge$cv0.9), as.numeric(tb.merge$cv0.75), method = "pearson")
cor.test(as.numeric(tb.merge$weighted.mean0.9), as.numeric(tb.merge$weighted.mean0.75), method = "pearson")



# ## OTU-based (i.e., randomly selected a certain number of OTU per sample)
# ## Method 2
# # here, we also use null model to test if the weighted mean tb is depended on the richness of each sample
# sample.num = dim(otu.1.sub)[2]
# 
# random.null = as.data.frame(matrix(nrow = 1, ncol = 4))
# colnames(random.null) = c("sample.id", "resample.richness", "mean.random", "weighted.mean.random")
# rownames(random.null) = "blank"
# 
# for (j in 1:dim(otu.1.sub)[2]) {
#   df.null = otu.tb[, c(1,1+j)]
#   df.null = df.null[ which(df.null[2] > 0), ]
#   sample.id = colnames(df.null[2])
#   
#   # random sample 30 OTUs (the min(richness) for all samples is 35, selecting this number according to the given case)
#   resample.richness = 30  # can change it to different number, such as 10, 20
# 
#   for (k in 1:100) { # repeat 100 times
#     
#     df.random = df.null[sample(nrow(df.null), resample.richness), ]
#     mean.random = mean(df.random$TB)
#     weighted.mean.random = weighted.mean(df.random$TB, df.random[[2]])
#     
#     # output results
#     random.sum = as.data.frame(matrix(nrow = 1, ncol = 4))
#     colnames(random.sum) = c("sample.id", "resample.richness", "mean.random", "weighted.mean.random")
#     random.sum[] <- c(sample.id, resample.richness, mean.random, weighted.mean.random)
#     random.merge = rbind(random.null, random.sum)
#     random.null = random.merge
#   }
# }
# 
# random.merge = random.merge[-1,]
# random.merge$mean.random = as.numeric(as.character(random.merge$mean.random))
# random.merge$weighted.mean.random = as.numeric(as.character(random.merge$weighted.mean.random))
# random.merge.sample.mean = aggregate(.~sample.id, random.merge[,-2], mean)
# 
# write.csv(random.merge.sample.mean,"tb.random.merge_900_rarefy.csv")








