library(ape) # read.tree branching.times vcv.phylo
library(R2jags) # jags

rhino.dat<-read.csv("rhino.csv")
rhino.tree<-read.tree("rhino.tree")

rhino.tree$edge.length<-rhino.tree$edge.length/max(branching.times(rhino.tree))

rhino.vcv <- vcv.phylo(rhino.tree)
ID<-diag(100)
VCV<-rhino.vcv

nc<-3  #number of chains
ni<-12000 #number of iterations
nb<-2000 #number of burnin
nt<-10 #thinning number

set.seed(12345) #we set a seed so results from the mcmc are perfectly replicable

sem.data<-list(BM=rhino.dat$BM,NL=rhino.dat$NL,DD=rhino.dat$DD,LS=rhino.dat$LS, RS=rhino.dat$RS,VCV=VCV,ID=ID,Nspec=100)

sem2.jg<-function(){
  #Structural equations 
  for (i in 1:Nspec) {
    muLS[i] <- alphaLS+betaBM*BM[i]
    muNL[i] <- alphaNL+betaBM2*BM[i]
    muDD[i] <- alphaDD+betaNL*NL[i]
    muRS[i] <- alphaRS+betaDD*DD[i]+betaLS*LS[i]
  }
  #Multivariate normal likelihoods
  LS[1:Nspec] ~ dmnorm(muLS[],TAUls)
  NL[1:Nspec] ~ dmnorm(muNL[],TAUnl)
  DD[1:Nspec] ~ dmnorm(muDD[],TAUdd)
  RS[1:Nspec] ~ dmnorm(muRS[],TAUrs)
  #Priors
  alphaLS ~ dnorm(0,1.0E-06)
  alphaNL ~ dnorm(0,1.0E-06)
  alphaDD ~ dnorm(0,1.0E-06)
  alphaRS ~ dnorm(0,1.0E-06)
  betaBM ~ dnorm(0,1.0E-06)
  betaBM2 ~ dnorm(0,1.0E-06)
  betaNL ~ dnorm(0,1.0E-06)
  betaDD ~ dnorm(0,1.0E-06)
  betaLS ~ dnorm(0,1.0E-06)
  lambdaLS ~ dunif(0,1)
  lambdaNL ~ dunif(0,1)
  lambdaDD ~ dunif(0,1)
  lambdaRS ~ dunif(0,1)
  tauLS ~ dgamma(1,1)
  tauNL ~ dgamma(1,1)
  tauDD ~ dgamma(1,1)
  tauRS ~ dgamma(1,1)
  sigmaLS <- 1/sqrt(tauLS)
  sigmaNL <- 1/sqrt(tauNL)
  sigmaDD <- 1/sqrt(tauDD)
  sigmaRS <- 1/sqrt(tauRS)
  #lambda computation
  MlamLS <- lambdaLS*VCV+(1-lambdaLS)*ID
  TAUls <- tauLS*inverse(MlamLS)
  MlamNL <- lambdaNL*VCV+(1-lambdaNL)*ID
  TAUnl <- tauNL*inverse(MlamNL)
  MlamDD <- lambdaDD*VCV+(1-lambdaDD)*ID
  TAUdd <- tauDD*inverse(MlamDD)
  MlamRS <- lambdaRS*VCV+(1-lambdaRS)*ID
  TAUrs <- tauRS*inverse(MlamRS)
}

params2 <- c("betaBM","betaBM2","betaNL","betaDD","betaLS")

sem2.mcmc<-jags(data=sem.data,model.file=sem2.jg,n.chains=3,n.iter=12000,n.burnin=2000,n.thin=10,parameters.to.save=params2)
save(file="sem2.Rdata", list="sem2.mcmc")
