#Note that this code snippet calls functions from GeneralFunctions.R
S0=36
sigma=0.4
capT=1
nSteps=50
stock_paths <- sim_stock_paths(S0,r,sigma,capT,nSteps,nPaths=2^16)
LSM <- Price_Amr_put_LSM(S0, K, r, sigma, capT, nSteps, nPaths = 2^16, stock_paths)
beta_LSM <- LSM[[3]]
Delta_LSM <- Price_Amr_put_DeltaLSM(S0, K, r, sigma, capT, nSteps, nPaths = 2^16, stock_paths)
beta_Delta_LSM <- Delta_LSM[[3]]
exercise_bound_LSM <- c(rep(0,50))
exercise_bound_LSM_delta <- c(rep(0,50))

for (j in 1:50) {
  if (length(exercise_intervals(beta_LSM,j,K)) > 0) {
    exercise_bound_LSM[j] <- max(exercise_intervals(beta_LSM,j,K)[,2])
  }
  if (length(exercise_intervals(beta = beta_Delta_LSM,j,K)) > 0) {
    exercise_bound_LSM_delta[j] <- max(exercise_intervals(beta = beta_Delta_LSM,j,K)[,2])
  }
}
