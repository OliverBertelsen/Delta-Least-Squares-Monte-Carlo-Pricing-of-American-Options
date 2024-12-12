#Note that this code snippet calls functions from GeneralFunctions.R
Upper_limit_new <- function(S0, K, r, sigma, capT, nSteps, nPaths_outer, nPaths_nested, beta, K_outer) {
  dt <- capT/nSteps
  first_part <- exp((r-0.5*sigma^2)*dt)
  second_part <- exp(sigma*sqrt(dt))
  discount_factor <- exp(-r*dt)
  delta0 <- c(rep(0,nPaths_outer))
  
  for (i in 1:nPaths_outer) {
    L_kB_k <- c(rep(0,nSteps+1))
    E_kL_kplusone <- c(rep(0,nSteps+1))
    l_k <- c(rep(0,nSteps+1))
    
    for (j in 2:nSteps) {
      intervals <- exercise_intervals(beta, j, K)
      state <- K_outer[i,j]
      if(exercise(intervals, state) & 0<state & state<K) {
        l_k[j] <- 1
        L_kB_k[j] <- payoff_put(state,K)*discount_factor^(j-1)
        if (j<nSteps) {
          payoffs <- c(rep(0, nPaths_nested))
          random_part <- rnorm((nSteps +1 -j) * nPaths_nested / 2)
          state_new_random_part <- t(matrix(c(random_part,-random_part), ncol = nPaths_nested))
          state_new_random_part <- first_part*second_part^state_new_random_part
          state_new_ <- state_new_random_part
          state_new_[,1] <- state*state_new_[,1]
          
          for (h in 1:(nSteps-j)) {
            state_new_[,h+1] <- state_new_[,h]*state_new_random_part[,h]
          }
          intervals <- exercise_intervals(beta, j=j+1, K)
          for (k in 1:nPaths_nested) {
            for (l in which(state_new_[k,(1:(nSteps -j))]<K)) {
              state_new <- state_new_[k,l]
              if(exercise(intervals, state=state_new)){
                payoffs[k] <- payoff_put(state_new,K)*discount_factor^l
                break
              }
            }
          }
          payoffs[which(payoffs==0)] <- payoff_put(state_new_[which(payoffs==0),(nSteps +1 -j)],K)*discount_factor^(nSteps +1 -j)
          E_kL_kplusone[j] <- mean(payoffs)*discount_factor^(j-1)
        } 
      } else {
        payoffs <- c(rep(0, nPaths_nested))
        random_part <- rnorm((nSteps +1 -j) * nPaths_nested / 2)
        state_new_random_part <- t(matrix(c(random_part,-random_part), ncol = nPaths_nested))
        state_new_random_part <- first_part*second_part^state_new_random_part
        state_new_ <- state_new_random_part
        state_new_[,1] <- state*state_new_[,1]
        if(j<nSteps){
        
          for (h in 1:(nSteps-j)) {
            state_new_[,h+1] <- state_new_[,h]*state_new_random_part[,h]
          }
        }
        
        for (k in 1:nPaths_nested) {
          for (l in which(state_new_[k,(1:(nSteps -j))]<K)) {
            state_new <- state_new_[k,l]
            if(exercise(intervals, state=state_new)){
              payoffs[k] <- payoff_put(state_new,K)*discount_factor^l
              break
            }
          }
        }
        payoffs[which(payoffs==0)] <- payoff_put(state_new_[which(payoffs==0),(nSteps +1 -j)],K)*discount_factor^(nSteps +1 -j)
        L_kB_k[j] <- mean(payoffs)*discount_factor^(j-1)
        E_kL_kplusone[j] <- L_kB_k[j]
      }
    }
    pi_k <- c(rep(0,nSteps+1))
    pi_k[2] <- L_kB_k[2]
    for (k in 3:nSteps) {
      pi_k[k] <- pi_k[k-1] + L_kB_k[k] - L_kB_k[k-1] - l_k[k-1]*(E_kL_kplusone[k-1]-L_kB_k[k-1])
    }
    k <- nSteps+1
    pi_k[nSteps+1] <- pi_k[k-1] + payoff_put(K_outer[i,nSteps+1],K)*discount_factor^nSteps - L_kB_k[k-1] - l_k[k-1]*(E_kL_kplusone[k-1]-L_kB_k[k-1])
    delta0[i] <- max(c(payoff_put(K_outer[i,],K)*discount_factor^(0:nSteps) - pi_k)[-1])
  }
  D0 <- mean(delta0)
  return(list(D0,delta0))
}

capT <- 1
nSteps <- 50*capT
dt <- capT/nSteps
K <- 40
r <- 0.06
sigma <- 0.2
S0 <- 36
nPaths <- 2^16
nPaths_outer <- 500
nPaths_nested <- 500

i_present <- #insert scenario number

Dual_Primal <- matrix(0, 1, 2, dimnames = list(NULL, c("LSM_D0", "Delta_LSM_D0")))
All_delta <- matrix(0,nPaths_outer,2,dimnames = list(NULL, c("LSM_all_D0", "Delta_LSM_all_D0")))

for(i in i_present){
  sigma = rep(c(rep(0.2,2),rep(0.4,2)),5)[i]
  S0 = sort(rep(seq(36,44,2),4))[i]
  capT = rep(1:2, 10)[i]
  nSteps = 50*rep(1:2, 10)[i]
  train <- sim_stock_paths(S0, r, sigma, capT, nSteps, nPaths)
  LSM <- Price_Amr_put_LSM(S0, K, r, sigma, capT, nSteps, nPaths, stock_paths=train)
  beta_LSM <- LSM[[3]]
  Delta_LSM <- Price_Amr_put_DeltaLSM(S0, K, r, sigma, capT, nSteps, nPaths, stock_paths=train)
  beta_Delta_LSM <- Delta_LSM[[3]]
  K_outer_for_both <- sim_stock_paths(S0, r, sigma, capT, nSteps, nPaths_outer)
  set.seed(2)
  PD_LSM <- Upper_limit_new(S0, K, r, sigma, capT, nSteps, nPaths_outer, nPaths_nested, beta = beta_LSM, K_outer = K_outer_for_both)
  Dual_Primal[,1] <- PD_LSM[[1]]
  set.seed(2)
  PD_Delta_LSM <- Upper_limit_new(S0, K, r, sigma, capT, nSteps, nPaths_outer, nPaths_nested, beta = beta_Delta_LSM, K_outer = K_outer_for_both)
  Dual_Primal[,2] <- PD_Delta_LSM[[1]]
  All_delta[,1] <- PD_LSM[[2]]
  All_delta[,2] <- PD_Delta_LSM[[2]]
}
print(Dual_Primal)
