normal_compound_poisson_process <- function(lambda, T, m, delta) {
  t <- 0
  jumps <- 0
  event_values <- c(0)
  event_times <- c(0)
  
  while (t < T) {
    t <- t + rexp(1, rate = lambda)
    jumps <- jumps + rnorm(1, mean = m, sd = delta)
    
    if (t < T) {
      event_values <- c(event_values, jumps)
      event_times <- c(event_times, t)
    }
  }
  list(event_values = event_values, event_times = event_times)
}

sim_MJD_paths <- function (S0, r, sigma, capT, nSteps, nPaths, m_jumps, sd_jumps, lambda) {
  dt <- capT/nSteps
  stock_paths <- matrix(0, nrow = nPaths, ncol = nSteps + 1)
  stock_paths[,1] <- S0
  k <- exp(m_jumps + 0.5 * sd_jumps^2) - 1
  Z <- t(matrix(rnorm((nSteps) * nPaths/2), ncol = nPaths/2))
  jumps <- matrix(0, nrow = nPaths/2, ncol=nSteps)
  
  for (t in 1:length(jumps)) {
    jumps[t] <- sum(normal_compound_poisson_process(lambda, dt, m_jumps, sd_jumps)$event_values)
  }
  
  for (t in 1:nSteps) {
    stock_paths[1:(nPaths/2), t+1] <- stock_paths[1:(nPaths/2),t]*exp((r-0.5*sigma^2-lambda*k)*dt+sigma*sqrt(dt)*Z[,t] + jumps[,t])
    stock_paths[((nPaths/2)+1):nPaths, t+1] <- stock_paths[((nPaths/2)+1):nPaths,t]*exp((r-0.5*sigma^2-lambda*k)*dt+sigma*sqrt(dt)*(-Z[,t]) - jumps[,t])
  }
  return(stock_paths)
}
