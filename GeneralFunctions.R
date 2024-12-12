payoff_put <- function(S, K) {
  pmax(K-S, 0)
}

sim_stock_paths <- function (S0, r, sigma, capT, nSteps, nPaths) {
  dt <- capT/nSteps
  stock_paths <- matrix(0, nrow = nPaths, ncol = nSteps + 1)
  stock_paths[,1] <- S0
  Z <- t(matrix(rnorm((nSteps) * nPaths/2), ncol = nPaths/2))
  for (t in 1:nSteps) {
    stock_paths[1:(nPaths/2), t+1] <- stock_paths[1:(nPaths/2),t]*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*Z[,t])
    stock_paths[((nPaths/2)+1):nPaths, t+1] <- stock_paths[((nPaths/2)+1):nPaths,t]*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*(-Z[,t]))
  }
  return(stock_paths)
}

Price_Amr_put_DeltaLSM <- function(S0, K, r, sigma, capT, nSteps, nPaths, stock_paths) {
  payoffs <- payoff_put(stock_paths[, nSteps+1], K)
  dt <- capT/nSteps
  discount_factor <- exp(-r*dt)
  continuation_values <- matrix(0, nrow = nPaths, ncol = nSteps+1)
  continuation_values[,nSteps+1] <- ifelse(payoff_put(stock_paths[,nSteps+1],K)>continuation_values[,nSteps+1], 1, 0)
  STauStar <- ifelse(continuation_values[,nSteps+1]==1,stock_paths[,nSteps+1],0)
  betaMatrix <- matrix(0, nrow = 4, ncol = nSteps)
  
  for (t in nSteps:2) {
    choose <- which(payoff_put(stock_paths[,t],K)>0)
    if (length(choose)>0) {
      STauStar <- STauStar*discount_factor
      ZS <- (STauStar[choose]/stock_paths[choose,t])*ifelse(K>STauStar[choose],-1,0)
      lambda <- sum(((payoffs[choose]*discount_factor))^2)/sum((ZS)^2)
      X <- stock_paths[choose,t]
      Y <- payoffs[choose]*discount_factor
      L1 <- X
      L2 <- X^2
      L3 <- X^3
      phi <- as.matrix(data.frame(L0 = X^0, L1 = L1, L2 = L2, L3 = L3))
      phid <- as.matrix(data.frame(L0= X*0, L1 = (X^0), L2 = (2*X), L3 = (3*X^2)))
      beta <- try(chol2inv(chol(t(phi)%*%phi+lambda*t(phid)%*%phid)) %*% (t(phi)%*%(payoffs[choose]*discount_factor)+lambda*(t(phid)%*%ZS)), silent = TRUE)
      if (inherits(beta, "try-error")  || any(is.na(beta))) {
        betaMatrix[,t] <- betaMatrix[,(t+1)]
      } else {
        betaMatrix[,t] <- beta
      }
      fitted_vals <- phi%*%betaMatrix[,t]
      continuation_values[choose,t] <- fitted_vals
    }
    payoffs <- ifelse(payoff_put(stock_paths[,t],K) > continuation_values[,t], payoff_put(stock_paths[,t],K), payoffs*discount_factor)
    continuation_values[,t] <- ifelse(payoff_put(stock_paths[,t],K) > continuation_values[,t], 1, 0)
    continuation_values[which(continuation_values[,t]==1),(t+1):(nSteps+1)] <- 0
    STauStar <- ifelse(continuation_values[,t]==1,stock_paths[,t],STauStar)
  }
  option_price <- mean(payoffs)*discount_factor 
  std_deviation <- sqrt(var((payoffs*discount_factor)))
  return(list(option_price,std_deviation, betaMatrix))
}

Price_Amr_put_LSM <- function(S0, K, r, sigma, capT, nSteps, nPaths, stock_paths) {
  payoffs <- payoff_put(stock_paths[, nSteps+1], K)
  dt <- capT/nSteps
  discount_factor <- exp(-r*dt)
  continuation_values <- matrix(0, nrow = nPaths, ncol = nSteps+1)
  betaMatrix <- matrix(0, nrow = 4, ncol = nSteps)
  
  for (t in nSteps:2) {
    choose <- which(payoff_put(stock_paths[,t],K)>0)
    if (length(choose)>0) {
      X <- stock_paths[choose,t]
      Y <- payoffs[choose]*discount_factor
      L1 <- X
      L2 <- X^2
      L3 <- X^3
      phi <- as.matrix(data.frame(L0 = X^0, L1 = L1, L2 = L2, L3 = L3))
      beta <- try(chol2inv(chol(t(phi)%*%phi))%*%t(phi)%*%payoffs[choose]*discount_factor, silent = TRUE)
      if (inherits(beta, "try-error")) {
        betaMatrix[,t] <- betaMatrix[,(t+1)]
      } else {
        betaMatrix[,t] <- beta
      }
      fitted_vals <- phi%*%betaMatrix[,t]
      continuation_values[choose,t] <- fitted_vals
    }
    payoffs <- ifelse(payoff_put(stock_paths[,t],K) > continuation_values[,t], payoff_put(stock_paths[,t],K), payoffs*discount_factor)
  }
  option_price <- mean(payoffs)*discount_factor 
  std_deviation <- sqrt(var((payoffs*discount_factor)))
  return(list(option_price,std_deviation,betaMatrix))
}

Price_Amr_put_LSM_Forward <- function(S0, K, r, sigma, capT, nSteps, nPaths, beta, stock_paths) {
  continuation_values <- matrix(0, nrow = nPaths, ncol = nSteps)
  payoffs <- rep(0,nPaths)
  discount_factor <- exp(-r*dt)
  
  for (t in 2:nSteps) {
    choose <- which(payoffs==0 & payoff_put(stock_paths[,t],K)>0)
    X <- stock_paths[choose,t]
    L1 <- X
    L2 <- X^2
    L3 <- X^3
    phi <- as.matrix(data.frame(L0 = X^0, L1 = L1, L2 = L2, L3 = L3))
    continuation_values[choose,t] <- phi%*%beta[,t]
    payoffs[choose] <- ifelse(payoff_put(stock_paths[choose,t],K) > continuation_values[choose,t], payoff_put(stock_paths[choose,t],K)*discount_factor^(t-1), payoffs[choose])
  }
  choose <- which(payoffs==0 & payoff_put(stock_paths[,nSteps+1],K)>0)
  payoffs[choose] <- payoff_put(stock_paths[choose,nSteps+1],K)*discount_factor^(nSteps)
  option_price <- mean(payoffs)
  std_deviation <- sqrt(var(payoffs))
  return(list(option_price,std_deviation))
}

CRR_Put<- function(S0,capT,strike,r,sigma,n,m){
  dt<-capT/n
  u<-exp(sigma*sqrt(dt))
  d<-1/u
  R<-exp(r*dt)
  q<-(R-d)/(u-d)
  PutB<-pmax(strike-S0*u^(0:n)*d^(n:0),0)
  
  for (i in n:1) {
    St<-S0*u^(0:(i-1))*d^((i-1):0)
    dummy<-pmax(strike-St,0)*(i%%m == 0) 
    temp<-pmax((q*PutB[2:(i+1)]+(1-q)*PutB[1:i])/R,dummy)
    PutB<-temp
  }
  return(PutB[1])
}

library(SobolSequence)

sim_sobol_paths <- function(S0, r, sigma, capT, nSteps, nPaths) {
  dt <- capT/nSteps
  sobol_seq <- sobolSequence.points(dimR=nSteps, dimF2=31, count=nPaths+1)[-1,]
  norm_seq <- qnorm(sobol_seq)
  S <- matrix(0, nrow = nPaths, ncol = nSteps + 1)
  S[, 1] <- S0
  for (i in 1:nSteps) {
    S[, i + 1] <- S[, i] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * norm_seq[, i])
  }
  return(S)
}

exercise_intervals <- function(beta, j, K) {
  coefficients <- c((beta[1,j] - K), (beta[2,j] + 1), beta[3,j], beta[4,j])
  roots <- polyroot(coefficients)
  exercise_state <- sort(Re(roots)[abs(Im(roots)) < 1e-6])
  exercise_state <- exercise_state[0<exercise_state & exercise_state<K]
  deriv <- 3*beta[4,j]*exercise_state^2+2*beta[3,j]*exercise_state+beta[2,j]
  
  if(length(deriv)==0) {
    if(K-mean(c(0,K))>beta[1,j]+beta[2,j]*mean(c(0,K))+beta[3,j]*mean(c(0,K))^2+beta[4,j]*mean(c(0,K))^3) {
      exercise <- matrix(c(0,K),nrow=1)
    } else {
      exercise <- integer(0)
    }
  } else if (length(deriv)==1) {
    if(deriv>-1) {
      exercise <- matrix(c(0,exercise_state),nrow=1)
    } else {
      exercise <- matrix(c(exercise_state,K),nrow=1)
    }
  } else if (length(deriv)==2) {
    if(deriv[1]>-1) {
      exercise <- matrix(c(0,exercise_state[1],exercise_state[2],K), nrow=2, byrow=TRUE)
    } else {
      exercise <- matrix(c(exercise_state[1],exercise_state[2]),nrow=1)
    }
  } else if (length(deriv)==3) {
    if (deriv[1]>-1) {
      exercise <- matrix(c(0,exercise_state[1],exercise_state[2],exercise_state[3]), nrow=2, byrow=TRUE)
    } else {
      exercise <- matrix(c(exercise_state[1],exercise_state[2],exercise_state[3],K), nrow=2, byrow=TRUE)
    }
  }
  return(exercise)
}

exercise <- function(intervals, state) {
  within_interval <- FALSE
  if(length(intervals)>0){
    within_interval <- any(apply(intervals, 1, function(row) {
      row[1] < state && state < row[2]
    }))
  }
  return(within_interval)
}
