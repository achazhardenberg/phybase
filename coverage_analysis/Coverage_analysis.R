library(ape)
library(caper)
library(nlme)
library(rjags)
library(R2jags)
library(runjags)

sem_sim.datadevtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R")
x1 <- array(rep(1, 10*11*4), dim=c(10, 11, 4))
rhino.dat<-read.csv("http://mpcm-evolution.com/OPM/Chapter8_OPM/download/rhino.csv")
rhino.tree<-read.tree("http://mpcm-evolution.com/OPM/Chapter8_OPM/download/rhino.tree")

rhino.tree$edge.length<-rhino.tree$edge.length/max(branching.times(rhino.tree))

rhino.vcv <- vcv.phylo(rhino.tree)
ID<-diag(100)
VCV<-rhino.vcv

nc<-3  #number of chains
ni<-12000 #number of iterations
nb<-2000 #number of burnin
nt<-10 #thinning number

set.seed(12345) #we set a seed so results from the mcmc are perfectly replicable

rhino.tree[[i]]$edge.length<-rhino.tree[[i]]$edge.length/max(branching.times(rhino.tree[[i]])) #all trees rescaled to total length of 1 to get correct lambda


for(i in 1:10) {
datamodel8<- '
data{
  #Structural equations 
  for (i in 1:Nspec) {
    muLS[i] <- alphaLS+betaBM*BM[i]
    muNL[i] <- alphaNL+betaBM2*BM[i]+betaRS*RS[i]
    muDD[i] <- alphaDD+betaNL*NL[i]
  }
  #Multivariate normal likelihoods
  LS[1:Nspec] ~ dmnorm(muLS[],TAUls)
  NL[1:Nspec] ~ dmnorm(muNL[],TAUnl)
  DD[1:Nspec] ~ dmnorm(muDD[],TAUdd)
  #lambda computation
  MlamLS <- lambdaLS*VCV+(1-lambdaLS)*ID
  TAUls <- tauLS*inverse(MlamLS)
  MlamNL <- lambdaNL*VCV+(1-lambdaNL)*ID
  TAUnl <- tauNL*inverse(MlamNL)
  MlamDD <- lambdaDD*VCV+(1-lambdaDD)*ID
  TAUdd <- tauDD*inverse(MlamDD)
}
model{
fake<-0
}
'
# parameters for simulations 

sigmaLS <- 0.2 # residual sd
sigmaNL <- 0.2 # residual sd
sigmaDD <- 0.2 # residual sd
tauLS <- 1 / (sigmaLS * sigmaLS)
tauNL <- 1 / (sigmaNL * sigmaNL)
tauDD <- 1 / (sigmaDD * sigmaDD)

# parameters are treated as data for the simulation step
data<-list(Nspec=100,RS=rhino.dat$RS,BM=rhino.dat$BM,alphaLS=0,alphaNL=0, alphaDD=0,
           betaBM=0.5,betaBM2=0.5,betaRS=0.5,betaNL=0.5,tauLS=tauLS,tauNL=tauNL,tauDD=tauDD,
           lambdaLS=0.8,lambdaNL=0.8,lambdaDD=0.8,VCV=VCV,ID=ID)
# run jags
out <- run.jags(datamodel8, data = data,monitor=c("LS","NL","DD","BM","RS"),sample=1, n.chains=1, summarise=FALSE)

Simulated <- coda::as.mcmc(out)
Simulated
dim(Simulated)
dat = as.vector(Simulated)
dat.mat<-matrix(dat, ncol=5,nrow=100, byrow = FALSE)
dat.dat <- as.data.frame(dat.mat)
colnames(dat.dat)<-c("LS","NL","DD","BM","RS")


params8 <- c("alphaLS", "alphaNL","alphaDD","betaBM","betaBM2","betaRS","betaNL","lambdaNL","lambdaLS","lambdaDD")

sem_sim.data<-list(BM=dat.dat$BM,NL=dat.dat$NL,DD=dat.dat$DD,LS=dat.dat$LS, RS=dat.dat$RS,VCV=VCV,ID=ID,Nspec=100)


sem8_sim.mcmc<-jags(data=sem_sim.data,model.file=sem8.jg,n.chains=3,n.iter=12000,n.burnin=2000, n.thin=10,parameters.to.save=params8)


sem8_val<-as.matrix(mcmctab(sem8_sim.mcmc)[2:5])

x[i,,]<-sem8_val
}

saveRDS(x, file = "x1.RDS")

#Rerun above for 10 simulation blocks (x1-x10)
#Then run the following to get the mean of the 10 simulation blocks

x1<-readRDS("x1.RDS")
x2<-readRDS("x2.RDS")
x3<-readRDS("x3.RDS")
x4<-readRDS("x4.RDS")
x5<-readRDS("x5.RDS")
x6<-readRDS("x6.RDS")
x7<-readRDS("x7.RDS")
x8<-readRDS("x8.RDS")
x9<-readRDS("x9.RDS")
x10<-readRDS("x10.RDS")

library(abind)

array<-abind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10, along=1)

val_means<-t(apply(array,c(3,2),mean)) #get means for all parameters over all 100 sims
rownames(val_means)<-c("alphaDD","alphaLS","alphaNL","betaBM","betaBM2","betaNL","betaRS","deviance","lambdaDD","lambdaLS","lambdaNL")
colnames(val_means)<-c("mean","se","lower","upper")

#               mean      se      lower      upper
# alphaDD    -0.00209 0.08053   -0.15968    0.15619
# alphaLS     0.00709 0.07975   -0.14906    0.16362
# alphaNL    -0.00374 0.08131   -0.16347    0.15620
# betaBM      0.49856 0.01295    0.47318    0.52400
# betaBM2     0.49900 0.01313    0.47323    0.52474
# betaNL      0.49895 0.01804    0.46355    0.53415
# betaRS      0.50192 0.01217    0.47810    0.52573
# deviance -269.77824 9.27219 -285.31105 -249.23306
# lambdaDD    0.91276 0.04110    0.80573    0.96352
# lambdaLS    0.91109 0.04067    0.80491    0.96193
# lambdaNL    0.90829 0.04268    0.79684    0.96099


#function to calculate bias as percentage of simulations where the target value is within the 95% credible interval
bias<-function(array,parameter,target){
  var.dat<-as.data.frame(array[,parameter,])
  var.dat$bias<-ifelse(target>=var.dat[,3],ifelse(target<=var.dat[,4],1,0),0)
  out<-sum(var.dat$bias)/100
  out
}


bias_alphaDD<-bias(array,1,0) #0.99
bias_alphaLS<-bias(array,2,0) #1
bias_alphaNL<-bias(array,3,0) #1
bias_betaBM<-bias(array,4,0.5) #0.97
bias_betaBM2<-bias(array,5,0.5) #0.99
bias_betaNL<-bias(array,6,0.5) #0.98
bias_betaRS<-bias(array,7,0.5) #0.98
bias_lambdaDD<-bias(array,9,0.8) #0.43
bias_lambdaLS<-bias(array,10,0.8) #0.43
bias_lambdaNL<-bias(array,11,0.8) #0.48
