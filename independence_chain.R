set.seed(451)
#
# Independence Chain Example
#

#
# Variables:
#
# y     = observed data
# n     = number of iterations 
# x.val = contains mixture parameter values
# f     = posterior density
# R     = computes Metropolis-Hastings ratio
# g     = proposal density
#

## Mixture data
# mixture.dat = read.table(file.choose(),header=TRUE)
mixture.dat = read.table(file='mixture.dat', header=TRUE);
y = mixture.dat$y

#
# Construct histogram of the mixture distribution
#

delta = .7;
par(mfrow=c(1,1))
x=seq(5,14,by=.01)
d = delta*dnorm(x,7,.5) + delta*dnorm(x,10,.5)
hist(y,breaks=20,freq=FALSE,main="Histogram of mixture data \n See Fig 7.1 in Givens and Hoeting",ylab="Density")
points(x,d,type="l")

#
# Set initial values
#

n           = 10000
x.val1      = NULL
x.val2      = NULL


#
# Define density functions: target and Metropolis-Hastings Ratio
#

f = function(x)   {prod(x*dnorm(y,7,0.5) + (1-x)*dnorm(y,10,0.5))}
R = function(xt,x){f(x)*g(xt)/(f(xt)*g(x))}


# Proposal density for Beta(1,1)
g = function(x){dbeta(x,1,1)}

x.val1[1] = rbeta(1,1,1)
for(i in 1:n){
      xt = x.val1[i]
      x = rbeta(1,1,1)
      p = min(R(xt,x),1)
      d = rbinom(1,1,p)
      x.val1[i+1] = x*d + xt*(1-d)
}

par(mfrow=c(2,2))
plot(x.val1[201:(n+1)],ylim=c(0,1),type="l",ylab="delta",xlab="t")
title("Sample path for Beta(1,1) Proposal Dist.")
hist(x.val1[201:(n+1)],breaks=20,xlab="delta", main="Hist. for Beta(1,1) Proposal Dist.")

#
# Proposal density for Beta(2,10)
#

g = function(x){dbeta(x,2,10)}
x.val2[1] = rbeta(1,2,10)
for(i in 1:n){
      xt = x.val2[i]
      x = rbeta(1,2,10)
      p = min(R(xt,x),1)
      d = rbinom(1,1,p)
      x.val2[i+1] = x*d + xt*(1-d)
}

plot(x.val2[201:(n+1)],ylim=c(0,1),type="l",ylab="delta",xlab="t")
title("Sample path for Beta(2,10) Proposal Dist.")
hist(x.val2[201:(n+1)],breaks=20,xlab="delta", main="Hist. for Beta(2,10) Proposal Dist.")
