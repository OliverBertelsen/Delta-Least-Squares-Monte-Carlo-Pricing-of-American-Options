#Note that this code snippet only creates the data used in creating Figure 5
CRR_Put<- function(S0,capT,strike,r,sigma,n,m){
  dt<-capT/n
  u<-exp(sigma*sqrt(dt))
  d<-1/u
  R<-exp(r*dt)
  q<-(R-d)/(u-d)
  PutB<-pmax(strike-S0*u^(0:n)*d^(n:0),0)
  exercise_area_below <- c(rep(0,n/40))
  for (i in n:1) {
    St<-S0*u^(0:(i-1))*d^((i-1):0)
    dummy<-pmax(strike-St,0)*(i%%m == 0) # in L&S-case: m=#steps per year/50
    if ((i%%m == 0)) {
      exercise_area_below[i/m] <- max(St[which((q*PutB[2:(i+1)]+(1-q)*PutB[1:i])/R<dummy)])
    }
    temp<-pmax((q*PutB[2:(i+1)]+(1-q)*PutB[1:i])/R,dummy)
    PutB<-temp
  }  
  return(list(PutB[1],exercise_area_below))
}
sigma <- rep(c(rep(0.2,2),rep(0.4,2)),5)
S0 <- sort(rep(seq(36,44,2),4))
capT <- rep(1:2, 10)
nSteps <- 50*capT
m <- 40
nCRR <- nSteps*40
CRR_Ex_T1 <- matrix(0, nrow = 10, ncol = 50)
for (i in c(1,3,5,7,9,11,13,15,17,19)) {
  CRR <- CRR_Put(S0[i],capT[i],strike = K,r,sigma[i],n = nCRR[i],m)
  CRR_Ex_T1[which(c(1,3,5,7,9,11,13,15,17,19)==i),] <- CRR[[2]]
}
