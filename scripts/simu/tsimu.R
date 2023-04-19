##################################################################################################
###                                                                                            ###
### Replication code for "L0 trend filtering" by Canhong Wen, Xueqin Wang and Aijun Zhang      ###
### This file contains the codes for simulation studies with t-distributed error               ###
### Updated on 13 April 2023                                                                   ###
###                                                                                            ###
##################################################################################################


source("utils_simu.R")

# Run the Monte Carlo replication with a server ------
if(!require(foreach)) install.packages('foreach')
if(!require(doMC)) install.packages('doMC')
library(foreach)
library(doMC)
registerDoMC(cores=10)
iter = 100
iter2 <- 10

# 1 Blocks example ----- 

nn=c(seq(32,224,32), seq(256, 2048, 256)); sigma=0.1;
ResultL0 <- ResultL1 <- ResultSp <- ResultPt <- ResultWs <- ResultNt <- ResultLc <- NULL
for (n in nn){
  cat("\nRunning n =", n)
  TabL0 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SingleRunL0TF(dgm = "Blocks", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabL1 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SingleRunL1TF(dgm = "Blocks", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabSp = foreach(mcseed = 1:iter2, .combine = rbind) %dopar%
    SingleRunSpline(dgm = "Blocks", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabPt = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SingleRunPeltTF(dgm = "Blocks", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabWs = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SingleRunwbsTF(dgm = "Blocks", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabNt = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SingleRunNotTF(dgm = "Blocks", n=n, sigma=sigma, seed=mcseed, dist="t")
  if (n<=224){
    TabLc = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
      SingleRunL0tfcTF(dgm = "Blocks", n=n, sigma=sigma, seed=mcseed, dist="t")
  }
  
  ResultL0 = rbind(ResultL0, c(n=n, colMeans(TabL0)))
  ResultL1 = rbind(ResultL1, c(n=n, colMeans(TabL1)))
  ResultSp = rbind(ResultSp, c(n=n, colMeans(TabSp)))
  ResultPt = rbind(ResultPt, c(n=n, colMeans(TabPt)))
  ResultWs = rbind(ResultWs, c(n=n, colMeans(TabWs)))
  ResultNt = rbind(ResultNt, c(n=n, colMeans(TabNt)))
  if (n<=224){
    ResultLc = rbind(ResultLc, c(n=n, colMeans(TabLc)))
  }else{
    ResultLc = rbind(ResultLc, c(n=n,Inf,Inf,Inf,Inf,Inf,Inf,Inf))
  }
}
save(ResultL0, ResultL1, ResultSp, ResultPt, ResultWs, ResultNt, ResultLc, file="tBlocks.RData")

# 2 Wave example ----- 
ResultL0 <- ResultL1 <- ResultSp <- ResultCp <- ResultNt <- ResultLc <- NULL
for (n in nn){
  cat("\nRunning n =", n)
  TabL0 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SingleRunL0TF(dgm = "Wave", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabL1 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SingleRunL1TF(dgm = "Wave", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabSp = foreach(mcseed = 1:iter2, .combine = rbind) %dopar%
    SingleRunSpline(dgm = "Wave", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabCp = foreach(mcseed = 1:iter, .combine = rbind) %dopar% 
    SingleRunCpopTF(dgm = "Wave", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabNt = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SingleRunNotTF(dgm = "Wave", n=n, sigma=sigma, seed=mcseed, dist="t")
  if (n<128){
    TabLc = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
      SingleRunL0tfcTF(dgm = "Wave", n=n, sigma=sigma, seed=mcseed, dist="t")
  }
  
  ResultL0 = rbind(ResultL0, c(n=n, colMeans(TabL0)))
  ResultL1 = rbind(ResultL1, c(n=n, colMeans(TabL1)))
  ResultSp = rbind(ResultSp, c(n=n, colMeans(TabSp)))
  ResultCp = rbind(ResultCp, c(n=n, colMeans(TabCp)))
  ResultNt = rbind(ResultNt, c(n=n, colMeans(TabNt)))
  if (n<128){
    ResultLc = rbind(ResultLc, c(n=n, colMeans(TabLc)))
  }else{
    ResultLc = rbind(ResultLc, c(n=n,Inf,Inf,Inf,Inf,Inf,Inf,Inf))
  }
  
}
save(ResultL0, ResultL1, ResultSp, ResultCp, ResultNt, ResultLc, file="tWave.RData")

# 3 Doppler example -----
nn <- seq(1280, 2048, 256)
ResultL0 <- ResultL1 <- ResultSp <- ResultNt <- ResultLc <- NULL
for (n in nn){
  cat("\nRunning n =", n)
  TabL0 = foreach(mcseed = 1:iter, .combine = rbind) %dopar% 
    SingleRunL0TF(dgm = "Doppler", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabL1 = foreach(mcseed = 1:iter, .combine = rbind) %dopar% 
    SingleRunL1TF(dgm = "Doppler", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabSp = foreach(mcseed = 1:iter2, .combine = rbind) %dopar% 
    SingleRunSpline(dgm = "Doppler", n=n, sigma=sigma, seed=mcseed, dist="t")
  TabNt = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SingleRunNotTF(dgm = "Doppler", n=n, sigma=sigma, seed=mcseed, dist="t")
  if (n<96){
    TabLc = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
      SingleRunL0tfcTF(dgm = "Doppler", n=n, sigma=sigma, seed=mcseed, dist="t")
  }
  
  ResultL0 = rbind(ResultL0, c(n=n, colMeans(TabL0)))
  ResultL1 = rbind(ResultL1, c(n=n, colMeans(TabL1)))
  ResultSp = rbind(ResultSp, c(n=n, colMeans(TabSp)))
  ResultNt = rbind(ResultNt, c(n=n, colMeans(TabNt)))
  if (n<96){
    ResultLc = rbind(ResultLc, c(n=n, colMeans(TabLc)))
  }else{
    ResultLc = rbind(ResultLc, c(n=n,Inf,Inf,Inf))
  }
}

save(ResultL0, ResultL1, ResultSp, ResultNt, ResultLc, file="tDoppler.RData")
