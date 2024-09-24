rm(list = ls())

library(TITAN2)
for (a in 1:2) {

Env <- read.csv("env.csv",row.names = 1)

if(a==1){
  tax <- read.csv("Specie_CY.csv",row.names = 1)
  tax <- tax[,-ncol(tax)]
  env <- Env[1:31,]
  res1 <-  titan(env[,5],
                tax,
                minSplt = 5,
                numPerm = 999,
                boot = TRUE,
                nBoot = 500,
                imax = FALSE,
                ivTot = FALSE,
                pur.cut = 0.95,
                rel.cut = 0.95,
                ncpus = 1,
                memory = FALSE)
  }else{
    tax <- read.csv("Specie_HH.csv",row.names = 1)
    tax <- tax[,-ncol(tax)]
    env <- Env[32:68,]
    res2 <-  titan(env[,5],
                   tax,
                   minSplt = 5,
                   numPerm = 999,
                   boot = TRUE,
                   nBoot = 500,
                   imax = FALSE,
                   ivTot = FALSE,
                   pur.cut = 0.95,
                   rel.cut = 0.95,
                   ncpus = 1,
                   memory = FALSE)
    }

}