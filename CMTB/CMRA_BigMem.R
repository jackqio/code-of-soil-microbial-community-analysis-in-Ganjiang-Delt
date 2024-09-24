library(vegan)
library(spaa)

rm(list=ls())
### for the first run, save the correlation matrix using bigmemory
########################
# setwd("~/Desktop/")
# input env data, standardize, calculate Eucldiean distance
env = read.csv("env.csv", header = T, row.names = 1)
env.std = decostand(env,method = "standardize")
env.dist = vegdist(env.std,"euclidean")

# input original otutable
otutable.16s.all = read.csv("otu_Flattening.csv", header = T, row.names = 1)
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

write.csv(otu.1.sub, "otu.1.sub.csv")

# Calculate pairwise spearman correlation as correlation-based synchrony at OTU-level 
cor = cor(t(otu.1.sub),method = "spearman")

# Transfer the spearman correlation as asynchrony [-1,1], as -1 is perfect synchrony and 1 is perfect asynchrony
cor.async = -(cor)

### save the correlation matrix
library(bigmemory)
# save to big matrix
#cormatrix=as.matrix(cor.async) # the big matrix you got, have to be a matrix file and not data.frame file
row.number = nrow(otu.1.sub) #set the row.number of cormatrix, i.e., 17811 in this case 
col.number = nrow(otu.1.sub) #set the col.number of cormatrix, i.e., 17811 in this case
backfile.name="cor.spearman17811.bin"
descripfile.name="cor.spearman17811.desc"
mbig = bigmemory::big.matrix(nrow = row.number, ncol = col.number,
                             type = "double", backingfile = backfile.name, descriptorfile = descripfile.name,
                             shared = TRUE) #can alse use "backingpath =" to set ouptut 
options(bigmemory.allow.dimnames=TRUE) # allow change col and row names of big matrix
mbig[]=cor.async # assign values from cor.async to big matrix
dimnames = rownames(cor.async)
write.csv(dimnames, "dimnames_cor.spearman17811.csv")
##################


### this step of CMRA calculation by resampling different number of OTUs is very slow, 
### can use different new RStudio windows at the same time for a parallel calculation
### When open a new RStudio window, run codes from here to calculate CMRA
library(vegan)
library(spaa)
library(bigmemory)
library(NST)

# setwd("~/Desktop/")
## call the big matrix in current or new RStudio windows
m.dat=bigmemory::attach.big.matrix(dget("cor.spearman17811.desc"))
# m.dat is big.matrix format but can be treated as original matrix

dimnames = read.csv("dimnames_cor.spearman17811.csv", row.names = 2)
options(bigmemory.allow.dimnames=TRUE) # allow change col and row names of big matrix
rownames(m.dat)=rownames(dimnames) # rename rows
colnames(m.dat)=rownames(dimnames) # rename columns

dat = read.csv("otu.1.sub.csv", row.names = 1)

# use null model to test if the weighted mean asychrony is depended on the richness of each sample
random.null = as.data.frame(matrix(nrow = 1, ncol = 5))
colnames(random.null) = c("sample.id", "resample.richness", "random.time", "mean.random", "weighted.mean.random")
rownames(random.null) = "blank"

for (j in 1:68) { ## j set the range of sample id (a totall of 472 in our case) to be calculated, e.g., 1:100,101:200,...401:472
  
  df.j = dat[which(dat[,j] != 0),j,drop=F] #remove OTU with relative abundance == 0
  sample.id = colnames(df.j[1])

  # random sample 900/700/500/300 OTUs ( the min(richness) for all samples is 964 )
  resample.richness = 30
  
  for (k in 1:100) { # k is the random.time, we need 100 times
    
    df.random = df.j[sample(nrow(df.j), resample.richness), , drop=F]
    random.time = k
    
    # generate OTU list in null community j
    otu.list.random = as.factor(row.names(df.random))
    
    # extract pairwise values of spearman correlation for otus at different abundance cutoffs
    value.random = m.dat[which(colnames(m.dat) %in% otu.list.random), which(rownames(m.dat) %in% otu.list.random)]
    value.random = dist.3col(value.random)
    
    # assign relative abundance for otus in name1 and name2
    value.random = merge(value.random, df.random[,1,drop=F], by.x="name1", by.y = "row.names", all.x = T)
    value.random = merge(value.random,df.random[,1,drop=F], by.x="name2", by.y = "row.names", all.x = T) 
    value.random$abund12 = value.random[,4] * value.random[,5]
    
    # mean asynchrony and community weighted mean asynchrony
    mean.random = mean(value.random[,"dis"])
    weighted.mean.random = weighted.mean(value.random[,"dis"], value.random[,"abund12"])
    
    # output results
    random.sum = c(sample.id, resample.richness, random.time, mean.random, weighted.mean.random)
    random.merge = rbind(random.null, random.sum)
    random.null = random.merge
  }
  write.csv(random.merge,"random.Asynchrony.csv")
}

random.merge = random.merge[-1,]
random.merge$mean.random = as.numeric(as.character(random.merge$mean.random))
random.merge$weighted.mean.random = as.numeric(as.character(random.merge$weighted.mean.random))
random.merge.sample.mean = aggregate(.~sample.id, random.merge[,-(2:3)], mean)

# output results, need to save as different file names if using parallel calculation
write.csv(random.merge.sample.mean,"Asynchrony.random.merge_1-100_filting.csv")




  
  
  