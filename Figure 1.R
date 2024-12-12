set.seed(100)
BlackScholesFormula  <- function (spot,timetomat,strike,r, q=0, sigma, opttype, greektype)
{
  d1<-(log(spot/strike)+ ((r-q)+0.5*sigma^2)*timetomat)/
    (sigma*sqrt(timetomat))
  d2<-d1-sigma*sqrt(timetomat)
  d3<-(-(r-sigma^2/2)*timetomat+log(strike/spot))/
    (sigma*sqrt(timetomat))
  d4<-(-(r-sigma^2/2)*timetomat-log(strike/spot))/
    (sigma*sqrt(timetomat))
  #european call (V)
  if (opttype==1 && greektype==1) result<-spot*exp(-q*timetomat)*pnorm(d1)-strike*exp(-r*timetomat)*pnorm(d2)
  #put (V)
  if (opttype==2 && greektype==1) result<-spot*exp(-q*timetomat)*pnorm(d1)-strike*exp(-r*timetomat)*pnorm(d2)-spot*exp(-q*timetomat)+strike*exp(-r*timetomat)
  #Greek = delta,  call (dV/dS)
  if (opttype==1 && greektype==2) result<-exp(-q*timetomat)*pnorm(d1)
  #Greek = delta,  put (dV/dS)
  if (opttype==2 && greektype==2) result<-exp(-q*timetomat)*(pnorm(d1)-1)
  #Greek = gamma  (d^2V/dS^2)
  if (greektype==3) result<-exp(-q*timetomat)*dnorm(d1)/
    (spot*sigma*sqrt(timetomat))
  #Greek = vega  (dV/dSigma)
  if (greektype==4) result<-exp(-q*timetomat)*spot*dnorm(d1)*sqrt(timetomat)
  #Digital (V)
  if (opttype==3 && greektype==1) result <- exp(-r*timetomat)*pnorm(-d2)
  #Digital delta
  if (opttype==3 && greektype==2) result <- exp(-r*timetomat)*(-dnorm(d2))/(spot*sigma*sqrt(timetomat))
  #No touch
  if (opttype==4 && greektype==1) result <- exp(-r*timetomat)*(pnorm(d3)-(spot/strike)^(1-2*r/(sigma^2))*pnorm(d4))
  #No touch delta
  if (opttype==4 && greektype==2) result <- exp(-r*timetomat)*((-1/(spot*sigma*sqrt(timetomat)))*dnorm(d3)-(1-(2*r/(sigma^2)))*spot^(-2*r/(sigma^2))*strike^(2*r/(sigma^2)-1)*pnorm(d4)-(spot/strike)^(1-2*r/(sigma^2))*dnorm(d4)*(1/(spot*sigma*sqrt(timetomat))))
  BlackScholesFormula<-result
}

S0<-100
r<-0.02
mu<-0.02
sigma<-0.2
sigma_hedge<-0.2
capT<-1
strike<-100
Nhedge<-252
Nrep<-10000

# HEDGE

St<-rep(S0, length=Nrep)
dt<-capT/Nhedge
initialoutlay<-BlackScholesFormula(S0,capT,strike, r,0,sigma,2,1)
Vpf<-rep(initialoutlay,length=Nrep)
a<-BlackScholesFormula(St,capT,strike, r,0,sigma_hedge,2,2)
b<-Vpf-a*St
for(i in 2:Nhedge){
  St<-St*exp((mu-0.5*sigma^2)*dt +sigma*sqrt(dt)*rnorm(Nrep))
  Vpf<-a*St+b*exp(dt*r)
  a<-BlackScholesFormula(St,(capT-(i-1)*dt),strike,r,0,sigma_hedge,2,2)
  b<-Vpf-a*St
}
St<-St*exp((mu-0.5*sigma^2)*dt +sigma*sqrt(dt)*rnorm(Nrep))
Vpf<-a*St+b*exp(dt*r)
optionpayoff<-pmax(strike-St,0)
hedgeerror<-Vpf-optionpayoff

# SUMMARY STATS & GRAPHS

print(paste("Initial investment =",round(initialoutlay,4)))
print(paste("Average discounted option payoff =",round(exp(-r*capT)*mean(optionpayoff),4)))
print(paste("Average discounted portfolio value =",round(exp(-r*capT)*mean(Vpf),4)))
plot(St,Vpf,col="blue",xlab="S(T)",ylab="Value of portfolio",main="Discrete hedging of a put-option,
ATM K=100",xlim=c(50,180),ylim=c(-5,60))
text(115,55,paste("# hegde points =",Nhedge),adj=0)
text(115,48, paste("Initial Investment=", round(initialoutlay,4)),adj=0)
text(115, 41, paste("Mean time0 option payoff=", round(exp(-r*capT)*mean(optionpayoff),4)),adj=0)
text(115, 34, paste("Mean time0 portfolio value=",round(exp(-r*capT)*mean(Vpf),4)),adj=0)
points(50:200,pmax(strike - 50:200,0),type='l',lwd=3)
print(paste("Mean and average hedgeerror=",mean(hedgeerror),sd(hedgeerror)))
