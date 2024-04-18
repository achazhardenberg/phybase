### Variability in trait model 

# As an example of a model taking account of data for which repeated measures 
# are available for each species we add white noise to the rhino data. 
# To do this each trait will be transformed to a distribution with mean equal 
# to the original species trait value and SD equal to a proportion (25%) of the among-species variance for that trait. 

##Simulation of repeated measures for each trait
set.seed(12345)

BM.new<-data.frame()
for(i in 1:length(rhino.dat$BM)){
  BM.mult<-sample(rnorm(500,mean=rhino.dat$BM[i],sd=sd(rhino.dat$BM)/4),10)
  BM.new<-rbind(BM.new,BM.mult)
}
names(BM.new)<-c("BM1","BM2","BM3","BM4","BM5","BM6","BM7","BM8","BM9","BM10")

NL.new<-data.frame()
for(i in 1:length(rhino.dat$NL)){
  NL.mult<-sample(rnorm(500,mean=rhino.dat$NL[i],sd=sd(rhino.dat$NL)/4),10)
  NL.new<-rbind(NL.new,NL.mult)
}
names(NL.new)<-c("NL1","NL2","NL3","NL4","NL5","NL6","NL7","NL8","NL9","NL10")

LS.new<-data.frame()
for(i in 1:length(rhino.dat$LS)){
  LS.mult<-sample(rnorm(500,mean=rhino.dat$LS[i],sd=sd(rhino.dat$LS)/4),10)
  LS.new<-rbind(LS.new,LS.mult)
}
names(LS.new)<-c("LS1","LS2","LS3","LS4","LS5","LS6","LS7","LS8","LS9","LS10")

DD.new<-data.frame()
for(i in 1:length(rhino.dat$DD)){
  DD.mult<-sample(rnorm(500,mean=rhino.dat$DD[i],sd=sd(rhino.dat$DD)/4),10)
  DD.new<-rbind(DD.new,DD.mult)
}
names(DD.new)<-c("DD1","DD2","DD3","DD4","DD5","DD6","DD7","DD8","DD9","DD10")

RS.new<-data.frame()
for(i in 1:length(rhino.dat$RS)){
  RS.mult<-sample(rnorm(500,mean=rhino.dat$RS[i],sd=sd(rhino.dat$RS)/4),10)
  RS.new<-rbind(RS.new,RS.mult)
}
names(RS.new)<-c("RS1","RS2","RS3","RS4","RS5","RS6","RS7","RS8","RS9","RS10")

data<-cbind(BM.new,NL.new,LS.new,DD.new,RS.new)

BM<-subset(data, select=BM1:BM10)
BM<-data.matrix(BM, rownames.force=NA)
BM2<-t(BM)
BMMulti<-as.vector(c(BM2[,1:100]))

NL<-subset(data, select=NL1:NL10)
NL<-data.matrix(NL, rownames.force=NA)
NL2<-t(NL)
NLMulti<-as.vector(c(NL2[,1:100]))

LS<-subset(data, select=LS1:LS10)
LS<-data.matrix(LS, rownames.force=NA)
LS2<-t(LS)
LSMulti<-as.vector(c(LS2[,1:100]))

DD<-subset(data, select=DD1:DD10)
DD<-data.matrix(DD, rownames.force=NA)
DD2<-t(DD)
DDMulti<-as.vector(c(DD2[,1:100]))

RS<-subset(data, select=RS1:RS10)
RS<-data.matrix(RS, rownames.force=NA)
RS2<-t(RS)
RSMulti<-as.vector(c(RS2[,1:100]))

dataMulti<-as.data.frame(cbind(BMMulti,NLMulti,LSMulti,DDMulti,RSMulti))

SP<-paste("s",sort(rep(1:100,10)),sep="")

dataMulti<-cbind(SP,dataMulti)

write.csv(dataMulti,file="RhinoMulti.csv")


#We read in the data and rescale the tree

library(ape)
library(caper)
library(nlme)
library(rjags)
library(R2jags)

rhinomulti.dat<-read.csv("RhinoMulti.csv")
names(rhinomulti.dat)<-c("X","SP","repBM","repNL","repLS","repDD","repRS")

rhino.tree<-read.tree("http://mpcm-evolution.com/OPM/Chapter8_OPM/download/rhino.tree")
rhino.tree$edge.length<-rhino.tree$edge.length/max(branching.times(rhino.tree)) 
#we need to rescale the total tree length to 1 to get a correct estimate of lambda


#And here the measurement error model:

sem8ME.jg<-function(){
  #Multivariate normal likelihoods
  LS[1:Nspec] ~ dmnorm(muLS[],TAUls[,])
  NL[1:Nspec] ~ dmnorm(muNL[],TAUnl[,])
  DD[1:Nspec] ~ dmnorm(muDD[],TAUdd[,])
  BM[1:Nspec] ~ dmnorm(muBM[],TAUbm[,])
  RS[1:Nspec] ~ dmnorm(muRS[],TAUrs[,])
  for (i in 1:Nspec) {
    #Prior mean of the specific levels muBM and muRS (not defined by the structural equations). Villemereur 2012, in his example on variable X fixes it to a single value (-10); we prefer sampling it from a uninformative prior distribution. As a working alternative (may be more stable): 
    #muBM[i]<-0
    #muRS[i]<-0
    muBM[i] ~ dnorm(0,1.0E-06)
    muRS[i] ~ dnorm(0,1.0E-06)
    #Structural equations 
    muLS[i] <- alphaLS+betaBM*BM[i]
    muNL[i] <- alphaNL+betaBM2*BM[i]+betaRS*RS[i]
    muDD[i] <- alphaDD+betaNL*NL[i]
    #Replicates normal independent distributions
    for (n in 1:NrepBM[i]) {
      repBM[i,n] ~ dnorm(BM[i], tauBM)
    }
    for (t in 1:NrepRS[i]) {
      repRS[i,t] ~ dnorm(RS[i], tauRS)
    }
    for (j in 1:NrepLS[i]) {
      repLS[i,j] ~ dnorm(LS[i], tauLS)
    }
    for (l in 1:NrepNL[i]) {
      repNL[i,l] ~ dnorm(NL[i], tauNL)
    }
    for (m in 1:NrepDD[i]) {
      repDD[i,m] ~ dnorm(DD[i], tauDD)
    }
  }
  #Priors
  alphaLS ~ dnorm(0,1.0E-06)
  alphaNL ~ dnorm(0,1.0E-06)
  alphaDD ~ dnorm(0,1.0E-06)
  betaBM ~ dnorm(0,1.0E-06)
  betaBM2 ~ dnorm(0,1.0E-06)
  betaRS ~ dnorm(0,1.0E-06)
  betaNL ~ dnorm(0,1.0E-06)
  lambdaLS ~ dunif(0,1)
  lambdaNL ~ dunif(0,1)
  lambdaDD ~ dunif(0,1)
  lambdaBM ~ dunif(0,1)
  lambdaRS ~ dunif(0,1)
  tauLS ~ dgamma(1,1)
  tauNL ~ dgamma(1,1)
  tauDD ~ dgamma(1,1)
  tauBM ~ dgamma(1,1)
  tauRS ~ dgamma(1,1)
  sigmaLS <- 1/sqrt(tauLS)
  sigmaNL <- 1/sqrt(tauNL)
  sigmaDD <- 1/sqrt(tauDD)
  sigmaBM <- 1/sqrt(tauBM)
  sigmaRS <- 1/sqrt(tauRS)
  #lambda computation
  MlamLS <- lambdaLS*VCV+(1-lambdaLS)*ID
  MlamNL <- lambdaNL*VCV+(1-lambdaNL)*ID
  MlamDD <- lambdaDD*VCV+(1-lambdaDD)*ID
  MlamBM <- lambdaBM*VCV+(1-lambdaBM)*ID
  MlamRS <- lambdaRS*VCV+(1-lambdaRS)*ID
  TAUls <- tauLS*inverse(MlamLS)
  TAUnl <- tauNL*inverse(MlamNL)
  TAUdd <- tauDD*inverse(MlamDD)
  TAUbm <- tauBM*inverse(MlamBM)
  TAUrs <- tauRS*inverse(MlamRS)
}

#here are the parameters to estimate:
params4 <- c("alphaLS", "alphaNL","alphaDD","betaBM","betaBM2","betaRS","betaNL","lambdaLS","lambdaNL","lambdaDD","lambdaBM","lambdaRS")

rhino.vcv <- vcv.phylo(rhino.tree)
ID<-diag(100)
VCV<-rhino.vcv

#For thew data list, we need to arrange the data in arrays with species in rows and replicates in columns:
sem8ME.data<-list(repBM=t(matrix(rhinomulti.dat$repBM,10,100)),repNL=t(matrix(rhinomulti.dat$repNL,10,100)),repDD=t(matrix(rhinomulti.dat$repDD,10,100)),repLS=t(matrix(rhinomulti.dat$repLS,10,100)),repRS=t(matrix(rhinomulti.dat$repRS,10,100)),VCV=VCV,ID=ID,NrepNL=rep(10,100),NrepLS=rep(10,100),NrepDD=rep(10,100),NrepBM=rep(10,100),NrepRS=rep(10,100),Nspec=100)

#Run the measurement error model: 
sem8ME.mcmc<-jags(model.file=sem8ME.jg,data=sem8ME.data,n.chains=3,n.iter=13000,n.burnin=3000,n.thin=10,parameters.to.save=params4)
sem8ME.mcmc


## Here below we show instead an example for the situation in which we do not have repeated measures for each species, 
## but we have average and standard error. This example code is based on the Rhinograd data for which mean and se 
## are known for each species. Data simulated using the following code: 
  
library(ape)
library(caper)
library(nlme)
library(rjags)
library(R2jags)

rhino.tree<-read.tree("http://mpcm-evolution.com/OPM/Chapter8_OPM/download/rhino.tree")

rhino.tree$edge.length<-rhino.tree$edge.length/max(branching.times(rhino.tree))
ID<-diag(100)
rhino.vcv <- vcv.phylo(rhino.tree)
VCV<-rhino.vcv

# Using the data file RhinoMulti.csv to estimate the mean and standard deviation of each trait value per species.
#This is for the alternative means of running a model with within-species variation, including mean and sd instead of multiple samples per species.

#Add new names for the file with just numbers without the "s" for species.
#this is necessary for the first loop over the various samples per species

data<-read.csv("RhinoMulti.csv",header=T)
NUM<-sort(rep(1:100,10))
data2<-data.frame(NUM,data[,3:7])

write.csv(data2,file="RhinoMulti2.csv")

data<-read.csv("RhinoMulti2.csv", header=TRUE)

BMmu<-tapply(data$BMMulti, data$NUM, mean)
BMsd<-tapply(data$BMMulti, data$NUM, sd)

NLmu<-tapply(data$NLMulti, data$NUM, mean)
NLsd<-tapply(data$NLMulti, data$NUM, sd)

LSmu<-tapply(data$LSMulti, data$NUM, mean)
LSsd<-tapply(data$LSMulti, data$NUM, sd)

DDmu<-tapply(data$DDMulti, data$NUM, mean)
DDsd<-tapply(data$DDMulti, data$NUM, sd)

RSmu<-tapply(data$RSMulti, data$NUM, mean)
RSsd<-tapply(data$RSMulti, data$NUM, sd)

dataSD<-data.frame(BMmu,BMsd,NLmu,NLsd,LSmu,LSsd,DDmu,DDsd,RSmu,RSsd)

# Here the actual model: 

sem8ME2.jg<-function(){
  # inclusion of trait variance
  for (i in 1:Nspec) {
    BMmu[i]~ dnorm(BM[i],tBM[i])
    LSmu[i]~ dnorm(LS[i],tLS[i])
    NLmu[i]~ dnorm(NL[i],tNL[i])
    DDmu[i]~ dnorm(DD[i],tDD[i])
    RSmu[i]~ dnorm(RS[i],tRS[i])
    muLS[i] <- alphaLS+betaBM*BM[i]
    muNL[i] <- alphaNL+betaBM2*BM[i]+betaRS*RS[i]
    muDD[i] <- alphaDD+betaNL*NL[i]
    muBM[i] ~ dnorm(0,1.0E-06)
    muRS[i] ~ dnorm(0,1.0E-06)
    tLS[i]<-(1/LSse[i])*(1/LSse[i])
    tNL[i]<-(1/NLse[i])*(1/NLse[i])
    tDD[i]<-(1/DDse[i])*(1/DDse[i])
    tRS[i]<-(1/RSse[i])*(1/RSse[i])
    tBM[i]<-(1/BMse[i])*(1/BMse[i])
  }
  #Multivariate normal likelihoods
  LS[1:Nspec] ~ dmnorm(muLS[],TAUls)
  BM[1:Nspec] ~ dmnorm(muBM[],TAUbm)
  NL[1:Nspec] ~ dmnorm(muNL[],TAUnl)
  DD[1:Nspec] ~ dmnorm(muDD[],TAUdd)
  RS[1:Nspec] ~ dmnorm(muRS[],TAUrs)
  #Priors
  alphaLS ~ dnorm(0,1.0E-06)
  alphaNL ~ dnorm(0,1.0E-06)
  alphaDD ~ dnorm(0,1.0E-06)
  betaBM ~ dnorm(0,1.0E-06)
  betaBM2 ~ dnorm(0,1.0E-06)
  betaNL ~ dnorm(0,1.0E-06)
  betaRS ~ dnorm(0,1.0E-06)
  lambdaLS ~ dunif(0,1)
  tauLS ~ dgamma(1,1)
  lambdaBM ~ dunif(0,1)
  tauBM ~ dgamma(1,1)
  lambdaNL ~ dunif(0,1)
  tauNL ~ dgamma(1,1)
  lambdaDD ~ dunif(0,1)
  tauDD ~ dgamma(1,1)
  lambdaRS ~ dunif(0,1)
  tauRS ~ dgamma(1,1)
  sigmaBM <- 1/sqrt(tauBM)
  sigmaLS <- 1/sqrt(tauLS)
  sigmaDD <- 1/sqrt(tauDD)
  sigmaNL <- 1/sqrt(tauNL)
  sigmaRS <- 1/sqrt(tauRS)
  #lambda computation
  MlamLS <- lambdaLS*VCV+(1-lambdaLS)*ID
  TAUls <- tauLS*inverse(MlamLS)
  MlamBM <- lambdaBM*VCV+(1-lambdaBM)*ID
  TAUbm <- tauBM*inverse(MlamBM)
  MlamDD <- lambdaDD*VCV+(1-lambdaDD)*ID
  TAUdd <- tauDD*inverse(MlamDD)
  MlamNL <- lambdaNL*VCV+(1-lambdaNL)*ID
  TAUnl <- tauNL*inverse(MlamNL)
  MlamRS <- lambdaRS*VCV+(1-lambdaRS)*ID
  TAUrs <- tauRS*inverse(MlamRS)
}

params <- c("betaBM","betaBM2","betaRS", "betaNL", "lambdaBM","lambdaLS","lambdaNL","lambdaRS","lambdaDD","sigmaBM","sigmaLS","sigmaNL","sigmaRS","sigmaDD")

sem8ME2.data<-list(BMmu=dataSD$BMmu,LSmu=dataSD$LSmu,NLmu=dataSD$NLmu,RSmu=dataSD$RSmu,DDmu=dataSD$DDmu,BMse=dataSD$BMsd/sqrt(9),LSse=dataSD$LSsd/sqrt(9),NLse=dataSD$NLsd/sqrt(9),RSse=dataSD$RSsd/sqrt(9),DDse=dataSD$DDsd/sqrt(9), VCV=VCV,ID=ID,Nspec=100)

sem8ME2.mcmc<-jags(data=sem8ME2.data,model.file=trait3.jg,n.chains=3,n.iter=13000,n.burnin=3000,n.thin=10,parameters.to.save=params)
sem8ME2.mcmc


