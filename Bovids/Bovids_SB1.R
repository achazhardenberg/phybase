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

set.seed(12345) #we set a seed so results from the mcmc are perfectly replicable

#model SB1 

bovidsSB1.jg<-function(){
  #Linear regression and multivariate normal likelihood
  for (i in 1:Nspec) {
    muBR[i] <- alpha1+beta1*BM[i]+beta2*S[i]
    muG[i] <- alpha2+beta3*BR[i]
    muL[i] <- alpha3+beta4*S[i]
  }
  BR[1:Nspec]~dmnorm(muBR[],TAUbr[,])
  G[1:Nspec]~dmnorm(muG[],TAUg[,])
  L[1:Nspec]~dmnorm(muL[],TAUl[,])
  #Priors
  alpha1 ~ dnorm(0,1.0E-06)
  alpha2 ~ dnorm(0,1.0E-06)
  alpha3 ~ dnorm(0,1.0E-06)
  beta1 ~ dnorm(0,1.0E-06)
  beta2 ~ dnorm(0,1.0E-06)
  beta3 ~ dnorm(0,1.0E-06)
  beta4 ~ dnorm(0,1.0E-06)
  lambdaBR ~ dunif(0,1)
  lambdaL ~ dunif(0,1)
  lambdaG ~ dunif(0,1)
  tauBR ~ dgamma(1,1)
  tauL ~ dgamma(1,1)
  tauG ~ dgamma(1,1)
  #Tree sampling and lambda computation
  for (k in 1:Ntree) {
    p[k] <- 1/Ntree
  }
  K~dcat(p[])
  #LAMBDA is a matrix with off-diagonal lambda value and 1 in the diagonal
  MlamBR <- lambdaBR*multiVCV[,,K]+(1-lambdaBR)*ID
  MlamL <- lambdaL*multiVCV[,,K]+(1-lambdaL)*ID
  MlamG <- lambdaG*multiVCV[,,K]+(1-lambdaG)*ID
  TAUbr <- tauBR*inverse(MlamBR)
  TAUl <- tauL*inverse(MlamL)
  TAUg <- tauG*inverse(MlamG)
}

params <- c("beta1","beta2","beta3","beta4","lambdaBR","lambdaL","lambdaG")

bovids.data<-list(BM=bovidsCutSc.dat$Adultwt,S=bovidsCutSc.dat$Groupsz,G=bovidsCutSc.dat$Gestation,L=bovidsCutSc.dat$Maxlongev,BR=bovidsCutSc.dat$Brain,multiVCV=multiVCV,ID=ID,Nspec=67,Ntree=100)

bovidsSB1.mcmc<-jags(data=bovids.data, model.file=bovidsSB1.jg,n.chains=3,n.iter=24000,n.burnin=4000,n.thin=10,parameters.to.save=params)

samples.SB1<- jags.samples(bovidsSB1.mcmc$model, 
                            c("WAIC","deviance"), 
                            type = "mean", 
                            n.iter = 24000,
                            n.burnin = 4000,
                            n.thin = 10)

samples.SB1$p_waic <- samples.SB1$WAIC
samples.SB1$waic <- samples.SB1$deviance + samples.SB1$p_waic
tmp <- sapply(samples.SB1, sum)
waic.SB1 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)

save(file="bovidsSB1.Rdata", list=c("bovidsSB1.mcmc", "samples.SB1","waic.SB1"))
