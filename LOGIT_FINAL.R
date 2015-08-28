#####################
## Bayesian Logit  ##
#####################

library(MCMCpack)
library(mvtnorm)
library(boot)


###############################
## log of posterior function ##
###############################

# Note: 
# log posterior = log binom prob  
# when inverse of sigma.0 prior is zero matrix
#
# else 
# log posterior = log mvnormal prob + log binom prob 

log.posterior <- function(beta, y, X, beta.0, sigma.0.inv)
{
 
 if (all(sigma.0.inv==0)){
    log.post<-sum(dbinom(x=y,size=1,prob=inv.logit(X%*%beta),log=T))
    }
  else{  
  log.post<-sum(dbinom(x=y,size=1,prob=inv.logit(X%*%beta),log=T))+
    dmvnorm(x=beta, mean=beta.0, sigma=solve(sigma.0.inv), log=T)
    }
  return(log.post) 
}

###########################################
## Proposal Variance Adaptation Function ##
###########################################

proposal.V <- function(y, X, beta.0, sigma.0.inv)
{
  #number of iterations
  nadapt<-5000
  accepted<-0
  #initial matrix
  V <- solve(t(X)%*%X)
  #constant to scale V matrix; will be adjusted in this function
  vsd=.1
  #initial beta
  beta.curr<-coef(glm(y~X+0, family="binomial"))
  lp.curr<-log.posterior(beta=as.numeric(beta.curr), y=y, X=X, beta.0=beta.0, sigma.0.inv=sigma.0.inv)
  #change vsd
  for(i in 1:nadapt){
    #propose beta star matrix - choose poropsal variance carefully 
    beta.star <- rmvnorm(n=1, mean = beta.curr, sigma=vsd*V)
    lp.curr<-log.posterior(beta=as.numeric(beta.curr), y=y, X=X, beta.0=beta.0, sigma.0.inv=sigma.0.inv)
    #log uniform
    log.uniform <- log(runif(n=1, min=0, max=1))
    #log posterior probabilities under beta t and beta star
    lp.star <- log.posterior(beta=as.numeric(beta.star), y=y, X=X, beta.0=beta.0, sigma.0.inv=sigma.0.inv)
    #test this value against log-uniform 
    prob <- lp.star-lp.curr
    #test criteria
    if(log.uniform<=prob)
      {
      lp.curr<-lp.star
      accepted <- accepted+1
      beta.curr<-beta.star
     }
  #adjustment to vsb
   if (i%%100 == 0){
      cat(paste0("Finished adapt ",i,"...\n"))
        if(accepted > 80){vsd<-vsd*4} 
        if(accepted < 80 && accepted > 50){vsd<-vsd*2}
        if(accepted < 30 && accepted > 10){vsd<-vsd*.5}
        if(accepted < 10 ){vsd<-vsd * .25}
      accepted<-0
    }
  }
  V <- vsd*V
  print(vsd)
  return(V)
}

##################################################
## Logit Function: Metropolis Sampler For Betas ##
##################################################

mh.logit<-function(y, X, beta.0, sigma.0.inv, samples, burnin, p, V)
{
 
  #intitialize beta samples matrix
  beta.t<-matrix(NA, ncol=p, nrow=samples+1)
  
  #starting state - normal logit coeffecients
  b.logit<-coef(glm(y~X+0, family="binomial"))
  beta.t[1,] <- as.vector(b.logit)
  
  #proposal beta 
  beta.star <- matrix(NA, ncol=p, nrow=1)
  
  #acceptance count
  accept <- 0
  
 
  
  #now sample beta
  for(i in  1:samples)
    {
    #propose beta.star
    beta.star[1,] <- rmvnorm(n=1, mean = beta.t[i,], sigma=V)
    #log uniform
    log.uniform <- log(runif(n=1, min=0, max=1))
    #log posterior probabilities under beta t and beta star
    p.star <- log.posterior(beta=beta.star[1,], y=y, X=X, beta.0=beta.0, sigma.0.inv=sigma.0.inv)
    p.t <- log.posterior(beta=beta.t[i,], y=y, X=X, beta.0=beta.0, sigma.0.inv=sigma.0.inv)
    #test this value against log-uniform 
    prob <- p.star-p.t
    #test criteria
    if(p.t ==-Inf && p.star!=-Inf)
      {
      beta.t[i+1,] <- beta.star[1,] 
      }
    else if(log.uniform<=prob && p.star !=-Inf && p.t !=-Inf)
    {
      beta.t[i+1,] <- beta.star[1,]
      if(i>burnin){
        accept <- accept+1
        }
    }
    else
    { 
      beta.t[i+1,] <- beta.t[i,]
    }
    #note this sampleing algorithm runs somewhat slow, even on a good computer
    if(i%%1000==0)
      {
      print(i)
      }
    }
  rate<- accept/(samples-burnin)
  cat(paste0("\n ACCEPTANCE RATE: ", rate))
  return(beta.t)
}

##########
## DATA ##
##########

  #load data 
  diabetes_subset <- read.table("E:/R/diabetes_subset.txt", header=T, quote="\"")
  y <- as.matrix(diabetes_subset$diabetes)
  #DROP THE BAD VARIABLE
  X <- diabetes_subset[,c(2:5, 8, 10, 11, 13:17)]
  A <- model.matrix( ~ location + gender + frame, data = diabetes_subset)
  intercept<-matrix(1, ncol=1, nrow =length(y[,1]))
  X <- as.matrix(cbind(intercept, X, A[,2:5]))
  #parameters in design matrix
  p<-as.numeric(ncol(X))
  #scale X - delete this if you don't want a scaled run
  for(i in 2:(p-4))
  {
    avg <- mean(X[,i])
    std <- sqrt(var(X[,i]))
    X[,i] <- (X[,i] - avg)/std 
  }

############
## PRIORS ##
############

  beta.0<-matrix(0, ncol=1, nrow=p)
  sigma.0.inv<-diag(0,p)

#######################
## PROPOSAL VARIANCE ##
#######################

  V <- proposal.V(y=y, X=X, beta.0=beta.0, sigma.0.inv=sigma.0.inv)

#################################
## SAMPLING ROUTINE PARAMETERS ##
#################################

  samples <- 100000
  burnin<- 1000

##########################
## RUN SAMPLING ROUTINE ##
##########################

  sample.beta<-mh.logit(y=y, X=X, beta.0=beta.0, sigma.0.inv=sigma.0.inv, samples=samples, burnin=burnin, p=p, V=V)

#########################################
## SUMMARY, PLOTS AND FURTHER ANALYSIS ##
#########################################

  summary(sample.beta)
  effectiveSize(sample.beta[burnin:samples,])
  plot(mcmc(sample.beta[,1:17]))


