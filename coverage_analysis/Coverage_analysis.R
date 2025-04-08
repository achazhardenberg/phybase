library(ape) # read.tree branching.times vcv.phylo 
library(R2jags) # jags 
library(runjags) # run.jags 
library(truncnorm) # rtruncnorm
library(cowplot) # plot_grid  

devtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R") 
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

sem8_cov.jg<-function(){
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
  #Priors
  alphaLS ~ dnorm(0,1)
  alphaNL ~ dnorm(0,1)
  alphaDD ~ dnorm(0,1)
  betaBM ~ dnorm(0.5,1)
  betaBM2 ~ dnorm(0.5,1)
  betaRS ~ dnorm(0.5,1)
  betaNL ~ dnorm(0.5,1)
  lambdaLS ~ dnorm(0.8,25);T(0,1)
  lambdaNL ~ dnorm(0.8,25);T(0,1)
  lambdaDD ~ dnorm(0.8,25);T(0,1)
  tauLS ~ dgamma(1,1)
  tauNL ~ dgamma(1,1)
  tauDD ~ dgamma(1,1)
  sigmaLS <- 1/sqrt(tauLS)
  sigmaNL <- 1/sqrt(tauNL)
  sigmaDD <- 1/sqrt(tauDD)
  #lambda computation
  MlamLS <- lambdaLS*VCV+(1-lambdaLS)*ID
  TAUls <- tauLS*inverse(MlamLS)
  MlamNL <- lambdaNL*VCV+(1-lambdaNL)*ID
  TAUnl <- tauNL*inverse(MlamNL)
  MlamDD <- lambdaDD*VCV+(1-lambdaDD)*ID
  TAUdd <- tauDD*inverse(MlamDD)
}
x <- array(NA, dim = c(100, 11, 5))
for(i in 1:100) {
# parameters are treated as data for the simulation step
data<-list(Nspec=100,RS=rhino.dat$RS,BM=rhino.dat$BM,alphaLS=rnorm(1,0,1),alphaNL=rnorm(1,0,1), alphaDD=rnorm(1,0,1),
           betaBM=rnorm(1,0.5,1),betaBM2=rnorm(1,0.5,1),betaRS=rnorm(1,0.5,1),betaNL=rnorm(1,0.5,1),tauLS=rgamma(1,1),tauNL=rgamma(1,1),tauDD=rgamma(1,1),
           lambdaLS=rtruncnorm(1, a = 0, b = 1, mean = 0.8, sd = sqrt(1/25)),lambdaNL=rtruncnorm(1, a = 0, b = 1, mean = 0.8, sd = sqrt(1/25)),lambdaDD=rtruncnorm(1, a = 0, b = 1, mean = 0.8, sd = sqrt(1/25)),VCV=VCV,ID=ID)
# run jags
out <- run.jags(datamodel8, data = data,monitor=c("LS","NL","DD","BM","RS"),sample=1, n.chains=1, summarise=FALSE)
Simulated <- coda::as.mcmc(out)
dat = as.vector(Simulated)
dat.mat<-matrix(dat, ncol=5,nrow=100, byrow = FALSE)
dat.dat <- as.data.frame(dat.mat)
colnames(dat.dat)<-c("LS","NL","DD","BM","RS")
params8 <- c("alphaLS", "alphaNL","alphaDD","betaBM","betaBM2","betaRS","betaNL","lambdaNL","lambdaLS","lambdaDD")
sem_sim.data<-list(BM=dat.dat$BM,NL=dat.dat$NL,DD=dat.dat$DD,LS=dat.dat$LS, RS=dat.dat$RS,VCV=VCV,ID=ID,Nspec=100)
sem8_sim.mcmc<-jags(data=sem_sim.data,model.file=sem8_cov.jg,n.chains=3,n.iter=12000,n.burnin=2000, n.thin=10,parameters.to.save=params8)
sem8.dat<-mcmctab(sem8_sim.mcmc)
sem8.dat$target<-c(data$alphaDD,data$alphaLS,data$alphaNL,data$betaBM,data$betaBM2,data$betaNL,data$betaRS,NA,data$lambdaDD,data$lambdaLS,data$lambdaNL)
sem8_val<-as.matrix(sem8.dat[2:6])
x[i,,]<-sem8_val
}
saveRDS(x, file = "x1c.RDS")

array<-readRDS("x1c.RDS")

#function to calculate bias as percentage of simulations where the target value is within the 95% credible interval
bias<-function(array,parameter){
  var.dat<-as.data.frame(array[,parameter,])
  var.dat$bias<-ifelse(var.dat[,5]>=var.dat[,3],ifelse(var.dat[,5]<=var.dat[,4],1,0),0)
  out<-sum(var.dat$bias)/100
  out
}


bias_alphaDD<-bias(array,1) #0.95
bias_alphaLS<-bias(array,2) #0.93
bias_alphaNL<-bias(array,3) #0.98
bias_betaBM<-bias(array,4) #0.96
bias_betaBM2<-bias(array,5)  #0.97
bias_betaNL<-bias(array,6) #0.97
bias_betaRS<-bias(array,7) #0.96
bias_lambdaDD<-bias(array,9) #0.95
bias_lambdaLS<-bias(array,10) #0.93
bias_lambdaNL<-bias(array,11) #0.96

plot_bias<-function(array,parameter){
  var.dat<-as.data.frame(array[,parameter,])
  var.dat$bias<-ifelse(var.dat[,5]>=var.dat[,3],ifelse(var.dat[,5]<=var.dat[,4],1,0),0)
  plot(var.dat[,1],var.dat[,5])
}

library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics

plot_bias <- function(array, parameter, Target) {
  # Convert array slice to dataframe
  var.dat <- as.data.frame(array[, parameter, ])
  var.dat$bias <- ifelse(var.dat[, 5] >= var.dat[, 3] & var.dat[, 5] <= var.dat[, 4], 1, 0)
  
  # Add a color column based on bias
  var.dat$color <- ifelse(var.dat$bias == 1, "lightblue", "salmon")
  
  # Determine axis limits
  axis_limits <- range(c(var.dat$V1, var.dat$V5, var.dat$V3, var.dat$V4))
  
  # Plot with ggplot2
  ggplot(var.dat, aes(x = V1, y = V5)) +
    # Confidence intervals as transparent solid bars
    geom_segment(aes(x = V1, xend = V1, y = V3, yend = V4, color = color), size = 2, alpha = 0.6) +
    # Black points for posterior means
    geom_point(color = "black", size = 0.5) +
    # 45-degree dotted line
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", linewidth = 0.5) +
    # Use colors directly from the dataframe
    scale_color_identity() +
    # Labels
    labs(
      x = Target,
      y = "Posterior mean",
    ) +
    # Set equal scaling and synchronized axes
    coord_cartesian(xlim = axis_limits, ylim = axis_limits) +
    theme_classic(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")
    )
}

p1<-plot_bias(x,1,"alphaDD")
p2<-plot_bias(x,2,"alphaLS")
p3<-plot_bias(x,3,"alphaNL")  
p4<-plot_bias(x,4,"betaBM")
p5<-plot_bias(x,5,"betaBM2")
p6<-plot_bias(x,6,"betaNL")
p7<-plot_bias(x,7,"betaRS")
p8<-plot_bias(x,9,"lambdaDD")
p9<-plot_bias(x,10,"lambdaLS")
p10<-plot_bias(x,11,"lambdaNL")

plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,ncol=3)+theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 

