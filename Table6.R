#Note that this code snippet calls functions from GeneralFunctions.R. Also note that the CRR-boundaries are produced just like in Figure5.R, hence left out.
#CRR_Ex_T1 and CRR_Ex_T2 represent the CRR exercise bounds.
fit_T1 <- matrix(0, nrow=10, ncol=50)
for (i in 1:10) {
  data_points <- CRR_Ex_T1[i,]
  X <- c(1:length(data_points))
  exp_model <- nls(data_points ~ a*exp(r*X) + c*exp(d*X), start = list(a = 20, r = 0.01, c=0.002, d=0.2))
  fit_T1[i,] <- predict(exp_model)
}
fit_T2 <-  matrix(0, nrow=10, ncol=100)
for (i in 1:10) {
  data_points <- CRR_Ex_T2[i,]
  X <- c(1:length(data_points))
  exp_model <- nls(data_points ~ a*exp(r*X) + c*exp(d*X), start = list(a = 20, r = 0.001, c=0.0001, d=0.1))
  fit_T2[i,] <- predict(exp_model)
}
capT <- 1
nSteps <- 50*capT
dt <- capT/nSteps
K <- 40
r <- 0.06
sigma <- 0.2
S0 <- 36
nPaths <- 2^16

Simulated_American <- matrix(0, 20, 4, dimnames = list(NULL, c("CRR as exercise", "SD", "CRR fit as exercise", "SD")))
set.seed(1)

for(i in 1:20){
  sigma = rep(c(rep(0.2,2),rep(0.4,2)),5)[i]
  S0 = sort(rep(seq(36,44,2),4))[i]
  capT = rep(1:2, 10)[i]
  nSteps = 50*rep(1:2, 10)[i]
  results <- matrix(0, 100, 2, dimnames = list(NULL, c("CRR as exercise", "CRR fit as exercise")))

  for (j in 1:100) {
    stock_paths <- sim_stock_paths(S0, r, sigma, capT, nSteps, nPaths)
    if (i%%2==0) {
      CRR_strategy <- CRR_Ex_T2[(i/2),]
      fit_strategy <- fit_T2[(i/2),]
    } else {
      CRR_strategy <- CRR_Ex_T1[((i+1)/2),]
      fit_strategy <- fit_T1[((i+1)/2),]
    }
    payoffs_CRR <- rep(0,nPaths)
    discount_factor <- exp(-r*dt)
    
    for (t in 2:nSteps) {
      choose <- which(payoffs_CRR == 0 & stock_paths[,t]<CRR_strategy[t-1])
      payoffs_CRR[choose] <- payoff_put(stock_paths[choose,t],K)*discount_factor^(t-1)
    }
    choose <- which(payoffs_CRR==0 & stock_paths[,nSteps+1]<CRR_strategy[nSteps])
    payoffs_CRR[choose] <- payoff_put(stock_paths[choose,nSteps+1],K)*discount_factor^(nSteps)
    option_price_CRR <- mean(payoffs_CRR)
    payoffs_fit <- rep(0,nPaths)
    discount_factor <- exp(-r*dt)
    
    for (t in 2:nSteps) {
      choose <- which(payoffs_fit == 0 & stock_paths[,t]<fit_strategy[t-1])
      payoffs_fit[choose] <- payoff_put(stock_paths[choose,t],K)*discount_factor^(t-1)
    }
    choose <- which(payoffs_fit==0 & stock_paths[,nSteps+1]<fit_strategy[nSteps])
    payoffs_fit[choose] <- payoff_put(stock_paths[choose,nSteps+1],K)*discount_factor^(nSteps)
    option_price_fit <- mean(payoffs_fit)
    results[j,] <- c(option_price_CRR, option_price_fit)
  }
  prices_and_se <- rep(0,4)
  prices_and_se[c(1,3)] <- colMeans(results)
  prices_and_se[c(2,4)] <- c(
    sqrt(var(results[,1])),
    sqrt(var(results[,2]))
  )
  Simulated_American[i,] <- c(prices_and_se)
}
