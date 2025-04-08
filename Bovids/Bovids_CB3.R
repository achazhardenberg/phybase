library(ape)
library(rjags)
library(R2jags)
bovids.dat<-read.csv("Bovidsdata.csv",header=T)
keeps<-c("Species","Brain","Gestation","Adultwt","Maxlongev","Groupsz")
bovidsCut.dat<-na.omit(bovids.dat[keeps]) #datafile with complete cases only
bovidsCut.dat$Adultwt<-log(bovidsCut.dat$Adultwt)
bovidsCut.dat$Brain<-log(bovidsCut.dat$Brain)
bovidsCut.dat$Gestation<-log(bovidsCut.dat$Gestation)
bovidsCut.dat$Maxlongev<-log(bovidsCut.dat$Maxlongev)

bovidsCut.dat$Greg<-ifelse(bovidsCut.dat$Groupsz>=6,1,0)
bovidsCut.dat$Groupsz<-log(bovidsCut.dat$Groupsz)

bovidsCutSc.dat<-data.frame(scale(bovidsCut.dat[2:6]))
bovidsCutSc.dat$Greg<-bovidsCut.dat$Greg

bovids.tree<-read.nexus("BovidsMultiphyCutnew.nex")
bovidsCut.tree<-bovids.tree[1:100] #100 trees only
for(i in 1:100) {
  bovidsCut.tree[[i]]$edge.length<-bovidsCut.tree[[i]]$edge.length/max(branching.times(bovidsCut.tree[[i]])) #all trees rescaled to total length of 1 to get correct lambda
}
for(i in 1:length(bovidsCut.tree)){
  multiVCV <- sapply(bovidsCut.tree, vcv.phylo, simplify="array")
}  #this creates the multi VCV array for all trees
ID <- diag(67)

#Set seed so results from the mcmc are perfectly replicable
set.seed(12345)

#Alternative model CB3

bovidsCB3.jg<-function(){
  #Linear regression and multivariate normal likelihood
  for (i in 1:Nspec) {
    muBR[i] <- alpha1+beta1*BM[i]+beta4*L[i]
    muS[i] <- alpha2+beta2*BM[i]
    muG[i] <- alpha3+beta3*BR[i]
  }
  BR[1:Nspec]~dmnorm(muBR[],TAUbr[,])
  S[1:Nspec]~dmnorm(muS[],TAUs[,])
  G[1:Nspec]~dmnorm(muG[],TAUg[,])
  
  #Priors
  alpha1 ~ dnorm(0,1.0E-06)
  alpha2 ~ dnorm(0,1.0E-06)
  alpha3 ~ dnorm(0,1.0E-06)
  beta1 ~ dnorm(0,1.0E-06)
  beta2 ~ dnorm(0,1.0E-06)
  beta3 ~ dnorm(0,1.0E-06)
  beta4 ~ dnorm(0,1.0E-06)
  lambdaBR ~ dunif(0,1)
  lambdaS ~ dunif(0,1)
  lambdaG ~ dunif(0,1)
  tauBR ~ dgamma(1,1)
  tauS ~ dgamma(1,1)
  tauG ~ dgamma(1,1)
  #Tree sampling and lambda computation
  for (k in 1:Ntree) {
    p[k] <- 1/Ntree
  }
  K~dcat(p[])
  #LAMBDA is a matrix with off-diagonal lambda value and 1 in the diagonal
  MlamBR <- lambdaBR*multiVCV[,,K]+(1-lambdaBR)*ID
  MlamS <- lambdaS*multiVCV[,,K]+(1-lambdaS)*ID
  MlamG <- lambdaG*multiVCV[,,K]+(1-lambdaG)*ID
  TAUbr <- tauBR*inverse(MlamBR)
  TAUs <- tauS*inverse(MlamS)
  TAUg <- tauG*inverse(MlamG)
}


params <- c("beta1","beta2","beta3","beta4","lambdaBR","lambdaS","lambdaG")

bovids.data<-list(BM=bovidsCutSc.dat$Adultwt,S=bovidsCutSc.dat$Groupsz,G=bovidsCutSc.dat$Gestation,L=bovidsCutSc.dat$Maxlongev,BR=bovidsCutSc.dat$Brain,multiVCV=multiVCV,ID=ID,Nspec=67,Ntree=100)

bovidsCB3.mcmc<-jags(data=bovids.data, model.file=bovidsCB3.jg,n.chains=3,n.iter=24000,n.burnin=4000,n.thin=10,parameters.to.save=params)

samples.CB3 <- jags.samples(bovidsCB3.mcmc$model, 
                            c("WAIC","deviance"), 
                            type = "mean", 
                            n.iter = 24000,
                            n.burnin = 4000,
                            n.thin = 10)

samples.CB3$p_waic <- samples.CB3$WAIC
samples.CB3$waic <- samples.CB3$deviance + samples.CB3$p_waic
tmp <- sapply(samples.CB3, sum)
waic.CB3 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)

save(file="bovidsCB3.Rdata", list=c("bovidsCB3.mcmc", "samples.CB3","waic.CB3"))
