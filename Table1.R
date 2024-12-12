#Note that this code snippet calls functions from GeneralFunctions.R
capT <- 1
nSteps <- 50*capT
dt <- capT/nSteps
K <- 40
r <- 0.06
sigma <- 0.2
S0 <- 36
nPaths <- 2^16

Simulated_American <- matrix(0, 20, 8, dimnames = list(NULL, c("LSM_IS", "std deviation","Delta_LSM_IS", "std deviation","LSM_OFS","std deviation","Delta_LSM_OFS","std deviation")))

for(i in 1:20){
  sigma = rep(c(rep(0.2,2),rep(0.4,2)),5)[i]
  S0 = sort(rep(seq(36,44,2),4))[i]
  capT = rep(1:2, 10)[i]
  nSteps = 50*rep(1:2, 10)[i]
  results <- matrix(0, 100, 4, dimnames = list(NULL, c("LSM_IS","Delta_LSM_IS","LSM_OFS","Delta_LSM_OFS")))
  
  for (j in 1:100) {
    simulated_paths_IS <-  sim_stock_paths(S0, r, sigma, capT, nSteps, nPaths)
    simulated_paths_OFS <- sim_stock_paths(S0, r, sigma, capT, nSteps, nPaths)
    LSM_IS <- Price_Amr_put_LSM(S0, K, r, sigma, capT, nSteps, nPaths, stock_paths = simulated_paths_IS)
    Delta_LSM_IS <- Price_Amr_put_DeltaLSM(S0, K, r, sigma, capT, nSteps, nPaths, stock_paths = simulated_paths_IS)
    LSM_OFS <- Price_Amr_put_LSM_Forward(S0, K, r, sigma, capT, nSteps, nPaths, beta = LSM_IS[[3]], stock_paths = simulated_paths_OFS)
    Delta_LSM_OFS <- Price_Amr_put_LSM_Forward(S0, K, r, sigma, capT, nSteps, nPaths, beta = Delta_LSM_IS[[3]], stock_paths = simulated_paths_OFS)
    results[j,] <- c(LSM_IS[[1]],Delta_LSM_IS[[1]],LSM_OFS[[1]],Delta_LSM_OFS[[1]])
  }
  prices_and_se <- rep(0,8)
  prices_and_se[c(1,3,5,7)] <- colMeans(results)
  prices_and_se[c(2,4,6,8)] <- c(
    sqrt(var(results[,1])),
    sqrt(var(results[,2])),
    sqrt(var(results[,3])),
    sqrt(var(results[,4]))
  )
  Simulated_American[i,] <- c(prices_and_se)
}

sigma <- rep(c(rep(0.2,2),rep(0.4,2)),5)
S0 <- sort(rep(seq(36,44,2),4))
capT <- rep(1:2, 10)
nSteps <- 50*capT
m <- 40
nCRR <- nSteps*40
CRR <- rep(0,20)

for (i in 1:20) {
  start.time <- Sys.time()
  CRR[i] <- CRR_Put(S0[i],capT[i],strike = K,r,sigma[i],n = nCRR[i],m)
  end.time <- Sys.time()
  print(end.time-start.time)
}

table1 <- cbind.data.frame(S = S0, sigma = sigma, T = capT, CRR, Simulated_American)
