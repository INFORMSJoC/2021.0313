##################################################################################################
###                                                                                            ###
### Replication code for "L0 trend filtering" by Canhong Wen, Xueqin Wang and Aijun Zhang      ###
### This file contains the codes for simulation studies with normal-distributed error for the  ###
### l0-MIP method with large sample size                                                       ###
### Updated on 13 April 2023                                                                   ###
###                                                                                            ###
##################################################################################################


# Run the Monte Carlo replication with a server ------
library(foreach)
library(doMC)
registerDoMC(cores=5)

# 2.1 Blocks example ----- 
nn=seq(32,256,32); sigma=0.1;
ResultLc <- NULL
for (n in nn){
  cat("\nRunning n =", n)
  TabLc = foreach(mcseed = 1:10, .combine = rbind) %dopar%
    SingleRunL0tfcTF(dgm = "Blocks", n=n, sigma=sigma, seed=mcseed, dist="normal")
  ResultLc = rbind(ResultLc, c(n=n, colMeans(TabLc)))
}
save(ResultLc, file="lcBlocks.RData")

# 2.2 Wave example ----- 
nn=seq(32,160,32); sigma=0.1;
ResultLc <- NULL
for (n in nn){
  cat("\nRunning n =", n)
  TabLc = foreach(mcseed = 1:10, .combine = rbind) %dopar%
    SingleRunL0tfcTF(dgm = "Wave", n=n, sigma=sigma, seed=mcseed, dist="normal")
  ResultLc = rbind(ResultLc, c(n=n, colMeans(TabLc)))
}
save(ResultLc, file="lcWave.RData")

# 2.3 Doppler example -----
nn=seq(32,96,32); sigma=0.1;
ResultLc <- NULL
for (n in nn){
  cat("\nRunning n =", n)
  TabLc = foreach(mcseed = 1:10, .combine = rbind) %dopar%
    SingleRunL0tfcTF(dgm = "Doppler", n=n, sigma=sigma, seed=mcseed, dist="normal")
  ResultLc = rbind(ResultLc, c(n=n, colMeans(TabLc)))
}

save(ResultLc, file="lcDoppler.RData")
