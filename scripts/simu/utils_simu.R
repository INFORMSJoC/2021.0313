##################################################################################################
###                                                                                            ###
### Replication code for "L0 trend filtering" by Canhong Wen, Xueqin Wang and Aijun Zhang      ###
### This file contains the utiles for simulation studies                                       ###
### Updated on 19 April 2023                                                                   ###
###                                                                                            ###
##################################################################################################

# Load Packages Needed -------

if(!require(gurobi)) install.packages('gurobi')
if(!require(pracma)) install.packages('pracma')
if(!require(AMIAS)) install.packages('AMIAS')
if(!require(wbs)) install.packages('wbs')
if(!require(not)) install.packages('not')
if(!require(changepoint)) install.packages('changepoint')
if(!require(freeknotsplines)) install.packages('freeknotsplines')
if(!require(cpop)) install.packages('cpop')
if(!require(genlasso)) install.packages('genlasso')

library(gurobi)
library(pracma)
library(AMIAS)
library(wbs)
library(not)
library(changepoint)
library(freeknotsplines)
library(cpop)
library(genlasso)

## The following three functions are used to generate simulated data in the Appendix of the papaer "L0 trend filtering".

# blocks: piecewise constant trend
# n is the sample size
# sigma is the constant in the formula $y = y0 + sigma \epsilon$
# seed is the seed for generating random sample
# dist is the distribution of the noise $\epsilon$, only standard normal distribution and t_(4) distribution is provided
SimuBlocks <- function (n, sigma = 0.1, seed=NA, dist=c("normal", "t"))
{
  if (!is.na(seed)) set.seed(seed)
  x = seq(1/n, 1,length.out = n)
  tau <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65,
           0.76, 0.78, 0.81)
  h <- c(0, 4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)/5
  A = sapply(tau, function(z) which(x>=z)[1])
  tau1=c(0,tau,1)
  y0 = 0*x
  for (j in 1:(length(A)+1)) y0[x>tau1[j] & x<=tau1[j+1]] = h[j]
  if(dist=="normal"){
    y <- y0 + sigma*rnorm(n)
  } else{
    y <- y0 + sigma*rt(n, 4)
  }
  
  return(list(y = y, x = x, y0 = y0, tau = tau, SetA = A))
}

# wave: piecewise linear trend
# n is the sample size
# sigma is the constant in the formula $y = y0 + sigma \epsilon$
# seed is the seed for generating random sample
# dist is the distribution of the noise $\epsilon$, only standard normal distribution and t_(4) distribution is provided
SimuWave <- function (n, sigma = 0.1, seed=NA, dist=c("normal", "t"))
{
  if (!is.na(seed)) set.seed(seed)
  x = seq(1/n, 1,length.out = n)
  tau = c(256, 512, 768, 1024, 1152, 1280, 1344, 1408)/1472
  nknot = length(tau)
  A = sapply(tau, function(z) which(x>=z)[1])
  tau=x[A]
  tau1 = c(0, tau, 1)
  h =  cumsum(c(1, (-1)^(1:nknot)*(1:nknot)))
  phi = rep(1, n)
  for (j in 1:(nknot+1)) phi = cbind(phi,pmin(pmax(x-tau1[j], 0), tau1[j+1]-tau1[j]))
  y0 = as.vector(phi%*%c(0, h))
  if(dist=="normal"){
    y <- y0 + sigma*rnorm(n)
  } else{
    y <- y0 + sigma*rt(n, 4)
  }
  return(list(y = y, x = x, y0 = y0, tau = tau, SetA = A))
}

# doppler: piecewise smooth trend with time-inhomogeneous smoothness
# n is the sample size
# sigma is the constant in the formula $y = y0 + sigma \epsilon$
# seed is the seed for generating random sample
# dist is the distribution of the noise $\epsilon$, only standard normal distribution and t_(4) distribution is provided
SimuDoppler <- function(n, sigma = 0.1, seed=NA, dist=c("normal", "t")){
  if (!is.na(seed)) set.seed(seed)
  x <- seq(1/n, 1,length.out = n)
  y0 <- sqrt(x*(1-x))*sin(2*pi*(1+0.05)/(x+0.05))
  if(dist=="normal"){
    y <- y0 + sigma*rnorm(n)
  } else{
    y <- y0 + sigma*rt(n, 4)
  }
  return(list(y = y, y0 = y0, x=x))
}

# Evalulate metrics for benchmark comparison
# beta is the estimate
# y0 is the true signal
# cpts is the estimated change points
# tcpt is the true change points
EvalMetrics <- function(beta, cpts=NULL, y0, tcpt = NULL){
  mse <- mean((beta-y0)^2)
  mad <- mean(abs(beta-y0))
  if(is.null(tcpt)){
    tab <- data.frame(MSE=mse, MAD=mad)
  }else{
    n <- length(beta)
    n.cpts <- length(cpts)
    segments.endpoints.true <- sort(unique(tcpt))
    segments.endpoints.est <- sort(unique(cpts))
    distm <- abs(matrix(rep(segments.endpoints.est, length(segments.endpoints.true)), nrow=length(segments.endpoints.est))
                 -matrix(rep(segments.endpoints.true, length(segments.endpoints.est)), nrow=length(segments.endpoints.est), byrow=TRUE))
    screening.dist <- max(apply(distm, 2, min)) * 100 / n # min distance for each true cpt
    precision.dist <- max(apply(distm, 1, min)) * 100 / n # min distance for each estimated cpt
    haus.dist <- max(screening.dist, precision.dist)
    tab <- data.frame(MSE=mse, MAD=mad, dS=screening.dist, dP=precision.dist, dH=haus.dist, nknot=n.cpts)
    #names(tab) <- c("MSE", "MAD", "d_S", "d_P", "d_H", "n.cpts")
  }
  return(tab)
}

# The difference matrix
DiffMat1 <- function(n){
  D = cbind(-diag(n-1),0) + cbind(0, diag(n-1))
  return(D)
}

# The k-th-oder difference matrix
DiffMat <- function(n, k=1){
  D = DiffMat1(n)
  if (k>1){for (l in 1:(k-1)) D=DiffMat1(n-l)%*%D}
  return(D)
}

## Run L0TF for just one run, which can be used in parallel computing
# dgm is the example we choose, only three examples are considered.
# other parameters are the same as SimuBlocks
SingleRunL0TF <- function(dgm = "Blocks", n=300, sigma=0.1, seed=NA, dist=c("normal", "t")){
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed, dist=dist)
    q=0; kmax=20;
  }
  if (dgm=="Wave"){
    data = SimuWave(n=n, sigma=sigma, seed=seed, dist=dist)
    q=1; kmax=20;
  }
  if (dgm=="Doppler"){
    data = SimuDoppler(n=n, sigma=sigma, seed=seed, dist=dist)
    q=2; kmax=20;
  }
  tic = Sys.time()
  if(q==0){
    resL0 = samias(as.numeric(data$y), D_type="tf0", kmax=kmax,  tmax=10)
  }else{
    resL0 = samias(as.numeric(data$y), D_type="tfq", q=q, kmax=kmax,  tmax=10)
  }
  toc = Sys.time()
  metric = EvalMetrics(resL0$alpha, resL0$A/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}

# Wrapped-up L1TF
L1TF <- function(data, n, q, maxdf) {
  resL1  <- trendfilter(pos=data$x, y=data$y, ord=q)
  idx <- which(resL1$df<=maxdf)
  
  bicL1 <- apply(resL1$beta[,idx], 2, function(beta) n*log(mean((data$y-beta)^2))) + 2*log(n)*resL1$df[idx]
  betaL1 <- resL1$beta[,which.min(bicL1)]
  knotL1 = data$x[which(abs(diff(betaL1, diff=q+1))>1e-5)+1]
  return(list(beta=betaL1, knot=knotL1))
}


# Run L1TF for just one run, which can be used in parallel computing
# dgm is the example we choose, only three examples are considered.
# other parameters are the same as SimuBlocks
SingleRunL1TF <- function(dgm = "Blocks", n=300, sigma=0.1, seed=NA, dist=c("normal", "t")){
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed, dist=dist)
    q=0; kmax=50;
  }
  if (dgm=="Wave"){
    data = SimuWave(n=n, sigma=sigma, seed=seed, dist=dist)
    q=1; kmax=50;
  }
  if (dgm=="Doppler"){
    data = SimuDoppler(n=n, sigma=sigma, seed=seed, dist=dist)
    q=2; kmax=100;
  }
  tic = Sys.time()
  resL1 = L1TF(data, n, q, maxdf=kmax)
  toc = Sys.time()
  metric = EvalMetrics(resL1$beta, resL1$knot, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}


# Run wbs for just one run, which can be used in parallel computing
# dgm is the example we choose, only three examples are considered.
# other parameters are the same as SimuBlocks
SingleRunwbsTF <- function(dgm = "Blocks", n=300, sigma=0.1, seed=NA, dist=c("normal", "t")){
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed, dist=dist)
    q=0; kmax=20;
  }
  if (dgm=="Wave"){
    stop("Wrong dgm in wbs!")
  }
  if (dgm=="Doppler"){
    stop("Wrong dgm in wbs!")
  }
  tic = Sys.time()
  resWbs = wbs::changepoints(wbs(data$y))
  toc = Sys.time()
  obj <- resWbs$cpt.th[[1]]
  metric = EvalMetrics(means.between.cpt(data$y, obj), obj/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}


# Run PELT for just one run, which can be used in parallel computing
# dgm is the example we choose, only three examples are considered.
# other parameters are the same as SimuBlocks
SingleRunPeltTF <- function(dgm = "Blocks", n=300, sigma=0.1, seed=NA, dist=c("normal", "t")){
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed, dist=dist)
    q=0; kmax=20;
  }
  if (dgm=="Wave"){
    stop("Wrong dgm in pelt!")
  }
  if (dgm=="Doppler"){
    stop("Wrong dgm in pelt!")
  }
  tic = Sys.time()
  resPelt = cpt.mean(data$y,method="PELT",class=TRUE)
  toc = Sys.time()
  obj = resPelt@cpts
  lis = c()
  value = resPelt@param.est$mean
  pivot = obj[1]
  for(ii in 1:length(obj)){
    lis <- c(lis,rep(value[ii], pivot))
    pivot <- obj[ii+1] - obj[ii]
  }
  metric = EvalMetrics(lis, obj[1:(length(obj)-1)]/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}


# Run CPOP for just one run, which can be used in parallel computing
# dgm is the example we choose, only three examples are considered.
# other parameters are the same as SimuBlocks
SingleRunCpopTF <- function(dgm = "Blocks", n=300, sigma=0.1, seed=NA, dist=c("normal", "t")){
  if (dgm=="Blocks"){
    stop("Wrong dgm in Cpop!")
  }
  if (dgm=="Wave"){
    data = SimuWave(n=n, sigma=sigma, seed=seed, dist=dist)
    q=0; kmax=20;
    tic = Sys.time()
    resCpop = cpop(data$y)
    toc = Sys.time()
  }
  if (dgm=="Doppler"){
    stop("Wrong dgm in Cpop!")
  }
  cpt <- resCpop@changepoints
  if(length(cpt)==2){
    cpt <- c()
  }else{
    cpt <- cpt[2:(length(cpt)-1)]
  }
  metric = EvalMetrics(estimate(resCpop)[,2], cpt/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}


# Run NOT for just one run, which can be used in parallel computing
# dgm is the example we choose, only three examples are considered.
# other parameters are the same as SimuBlocks
SingleRunNotTF <- function(dgm = "Blocks", n=300, sigma=0.1, seed=NA, dist=c("normal", "t")){
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed, dist=dist)
    q=0; kmax=20;
    tic = Sys.time()
    resNot = not(data$y, contrast="pcwsConstMean")
    toc = Sys.time()
  }
  if (dgm=="Wave"){
    data = SimuWave(n=n, sigma=sigma, seed=seed, dist=dist)
    q=1; kmax=20;
    tic = Sys.time()
    resNot = not(data$y, contrast="pcwsLinContMean")
    toc = Sys.time()
  }
  if (dgm=="Doppler"){
    data = SimuDoppler(n=n, sigma=sigma, seed=seed, dist=dist)
    q=2; kmax=20;
    tic = Sys.time()
    resNot = not(data$y, contrast="pcwsQuadMean")
    toc = Sys.time()
  }
  obj = features(resNot)$cpt
  lis = predict(resNot)
  metric = EvalMetrics(lis, obj/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}

# Functions adopted from https://github.com/yshin12/sparseHP. This implement the l0-MIP method.
l0tfc <- function(y=y, D_matrix, l0constraint=l0constraint,M=M,l2penalty=l2penalty,...){
  
  n <- length(y)
  rD <- nrow(D_matrix)
  A1 <- cbind(D_matrix,-M*diag(rD))
  A2 <- cbind(-D_matrix,-M*diag(rD))
  A3 <- matrix(c(rep(0,n),rep(1,(rD))),nrow=1)
  A_matrix <- rbind(A1,A2,A3)
  A_c <- ncol(A_matrix)
  A_r <- nrow(A_matrix)
  
  b_vector <- c(rep(0,(A_r-1)),l0constraint)
  
  c_vector <- c(-2*y,rep(0,rD))
  
  Q_matrix <- matrix(0,nrow=A_c,ncol=A_c)
  Q_matrix[1:n,1:n] <- diag(n) + l2penalty*(t(D_matrix)%*%D_matrix)
  
  x_ub <- c(rep(max(y),n),rep(1,rD))
  x_lb <- c(rep(min(y),n),rep(0,rD))
  
  x_type <- rep('C',n)
  z_type <- rep('B',rD)
  
  model <- list()
  
  model$A          <- A_matrix
  model$rhs        <-  b_vector
  model$sense      <- c(rep('<',A_r-1),"<=")
  model$obj        <- c_vector
  model$Q          <- Q_matrix
  model$modelsense <- 'min'
  model$lb         <- x_lb
  model$ub         <- x_ub
  model$vtype      <- c(x_type,z_type)
  
  params <- list(OutputFlag=0)
  
  result <- gurobi(model, params)
  
  return(result)
  
  rm(model, result, params)
}

# Choosing the optimal tuning parameter for l0-MIP using BIC.
# y is the observed data
# D is the difference matrix we consider
# kmax is the maximum number of knots
# q is the order of difference matrix
BIC_l0tfc <- function(y, D, kmax, q){
  bwl=3
  n = length(y)
  Y = movavg(y, bwl, type="s")
  M = max(abs(diff(x = Y, differences = 2)))
  bic <- c()
  fine <- matrix(0,nrow=n,ncol=kmax)
  knots <- list()
  for(k in 1:kmax){
    print(k)
    temp <- l0tfc(y,D,l0constraint=k,M,l2penalty=0)
    fine[,k] <- temp$x[1:n]
    knots[[k]] <- which(temp$x[(n+1):(2*n-1)]==1)+1
    sighat <- median(abs(diff(y, diff = q+1)))/(qnorm(3/4)*sqrt(choose(2*(q+1), q+1)))
    bbb <- sum((y-fine[,k])^2)/sighat^2 + 2*log(n)*(length(knots[[k]])+q+1)
    bic <- c(bic, bbb)
  }
  idx <- which.min(bic)
  list(knot=knots[[idx]], fit=fine[,idx], bic=bic[idx])
}


# Run l0-MIP for just one run, which can be used in parallel computing
# Main function to perform $\ell$_0 trend filtering using $\ell_0$-MIP algorithm
# dgm is the example we choose, only three examples are considered.
# other parameters are the same as SimuBlocks
SingleRunL0tfcTF <- function(dgm = "Blocks", n=300, sigma=0.1, seed=NA, dist=c("normal", "t")){
  
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed, dist=dist)
    q=0; kmax=20;
    tic = Sys.time()
    D_matrix = DiffMat1(n)
    resL0tfc = BIC_l0tfc(data$y,D_matrix,kmax,q)
    toc = Sys.time()
  }
  if (dgm=="Wave"){
    data = SimuWave(n=n, sigma=sigma, seed=seed, dist=dist)
    q=1; kmax=20;
    tic = Sys.time()
    D_matrix = DiffMat(n,q)
    resL0tfc = BIC_l0tfc(data$y,D_matrix,kmax,q)
    toc = Sys.time()
  }
  if (dgm=="Doppler"){
    data = SimuDoppler(n=n, sigma=sigma, seed=seed, dist=dist)
    q=2; kmax=20;
    tic = Sys.time()
    D_matrix = DiffMat(n,q)
    resL0tfc = BIC_l0tfc(data$y,D_matrix,kmax,q)
    toc = Sys.time()
  }
  
  metric = EvalMetrics(resL0tfc$fit, resL0tfc$knot/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}

# Run Spline for just one run, which can be used in parallel computing
# dgm is the example we choose, only three examples are considered.
# other parameters are the same as SimuBlocks
SingleRunSpline <- function(dgm = "Blocks", n=300, sigma=0.1, seed=NA, dist=c("normal", "t")){
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed, dist=dist)
    q=0; kmax=20;
  }
  if (dgm=="Wave"){
    data = SimuWave(n=n, sigma=sigma, seed=seed, dist=dist)
    q=1; kmax=20;
  }
  if (dgm=="Doppler"){
    data = SimuDoppler(n=n, sigma=sigma, seed=seed, dist=dist)
    q=2; kmax=20;
  }
  tic = Sys.time()
  resSpline <- fit.search.numknots(data$x, data$y, alg = "PS", degree = q, maxknot = kmax, knotnumcrit = "BIC")
  SplineBeta <- fitted.freekt(resSpline)
  SplineKnot <- attr(resSpline,"optknot")
  toc = Sys.time()
  metric = EvalMetrics(SplineBeta, SplineKnot, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}
