#' Code to run the Markov Switching regime shift models in JAGS
#'
#' @author Gavin Fay


#memory.limit(size=4084)
#setwd("gfay/sandbox/regime/")
if (Sys.info()[['sysname']] == "Windows") setwd("C:/Research/regime/")
if (Sys.info()[['sysname']] == "Darwin") setwd("research/regimeswitch/")

library(rjags)

for (isp in 1:5)
{
#isp=5  
  if (isp==1) infile = "gb_cod.dat"
  if (isp==2) infile = "gom_had.dat"
  if (isp==3) infile = "gb_ytf.dat"
  if (isp==4) infile = "sne_winf.dat" 
  if (isp==5) infile = "gom_cod2.dat"
  infile <- paste("data/",infile,sep="")

  Nobs = scan(infile,skip=1,n=1)
  Nreg = scan(infile,skip=0,n=1)
  Yobs = scan(infile,skip=2,n=Nobs)
  obsErr = scan(infile,skip=3,n=Nobs)
  Catch = scan(infile,skip=4,n=Nobs)
  model1.data = list(N=Nobs,I=Yobs,CV=obsErr,C=Catch)
 
  model2.jags2save = c("K","r","Sigma2","Tau2","TransProbs","Regime","B",
                     "RegimeNext","FMSY","NextCatch","q")

  inits <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),
                r2frac=runif(1,0.2,0.8), k=0.00004, isigma2=100, itau2=100,
                TransProbs=c(0.97,0.97),Regime=sample(1:2,Nobs,replace=TRUE),
                q=-0.01)
  inits2 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),
                 r2frac=runif(1,0.2,0.8), k=0.00003, isigma2=80, itau2=150,
                 TransProbs=c(0.7,0.8),Regime=sample(1:2,Nobs,replace=TRUE),
                 q=-0.04)
  inits3 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),
                 r2frac=runif(1,0.2,0.8), k=0.00005, isigma2=120, itau2=50,
                 TransProbs=c(0.9,0.7),Regime=sample(1:2,Nobs,replace=TRUE),
                 q=-0.006)

  inits <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),
                r2frac=runif(1,0.2,0.8), k=0.00004, isigma2=100, itau2=100,
                TransProbs=c(0.97,0.97),Regime=sample(1:2,Nobs,replace=TRUE),
                iq=800,idep=0)
  inits2 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),
                 r2frac=runif(1,0.2,0.8), k=0.00003, isigma2=80, itau2=150,
                 TransProbs=c(0.7,0.8),Regime=sample(1:2,Nobs,replace=TRUE),
                 iq=600,idep=-0.5)
  inits3 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),
                 r2frac=runif(1,0.2,0.8), k=0.00005, isigma2=120, itau2=50,
                 TransProbs=c(0.9,0.7),Regime=sample(1:2,Nobs,replace=TRUE),
                 iq=450,idep=-0.2) 
  # Some alternative initial q values for some of the stocks
  if (isp ==1 | isp ==4)
   {
    inits <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),
                  r2frac=runif(1,0.2,0.8), k=0.00004, isigma2=100, itau2=100,
                  TransProbs=c(0.97,0.97),Regime=sample(1:2,Nobs,replace=TRUE),
                  iq=8000,idep=0)
    inits2 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),
                   r2frac=runif(1,0.2,0.8), k=0.00003, isigma2=80, itau2=150,
                   TransProbs=c(0.7,0.8),Regime=sample(1:2,Nobs,replace=TRUE),
                   iq=6000,idep=-0.5)
    inits3 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),
                   r2frac=runif(1,0.2,0.8), k=0.00005, isigma2=120, itau2=50,
                   TransProbs=c(0.9,0.7),Regime=sample(1:2,Nobs,replace=TRUE),
                   iq=4500,idep=-0.2)
   }

  model2.init <- list(inits,inits2,inits3)

  #jags <- jags.model("C:/Research/regime/schaef_jags_notatK.txt",
  #jags <- jags.model("regime_schaef_jags_notatK.txt",
  jags <- jags.model("regime_schaef_jags_notatK_2.txt",
                   data = model1.data,
                   inits = model2.init,
                   n.chains = 3,
                   n.adapt = 100)

  update(jags,500000)  ####burn-in
  mod2 <- jags.samples(jags,model2.jags2save,n.iter=1000000,thin=1000)
  #mod2 <- jags.samples(jags,model2.jags2save,n.iter=10000,thin=10)
  #dic.mod2 <- dic.samples(jags,n.iter=10000,thin=10,type="pD")
  dic.mod2 <- dic.samples(jags,n.iter=1000000,thin=1000,type="pD")

  #xx is a list containing sampled vectors of desired nodes
  mod2.results <- NULL
  mod2.results$mod <- mod2
  mod2.results$dic <- dic.mod2

  ##save(mod2.results,file=paste(isp,"_mod2_posterior.Rdata",sep=""))
  ##save(mod2.results,file=paste(isp,"_mod2_big_posterior.Rdata",sep=""))
  ##save(mod2.results,file=paste(isp,"_mod2_test_posterior.Rdata",sep=""))
  ##save(mod2.results,file=paste(isp,"_mod2_20140604_2_posterior.Rdata",sep=""))
  save(mod2.results,file=paste(isp,"_mod2_big_posterior_20140618.Rdata",sep=""))

}


######## NO REGIMES #############

#setwd("C:/Research/regime/")
#setwd("gfay/sandbox/regime/")

#for (isp in c(2,3:5))
#{
  
# if (isp==1) infile = "gb_cod.dat"
# if (isp==2) infile = "gom_had.dat"
# if (isp==3) infile = "gb_ytf.dat"
# if (isp==4) infile = "sne_winf.dat" 
#  if (isp==5) infile = "gom_cod2.dat"

#Nobs = scan(infile,skip=1,n=1)
#Nreg = scan(infile,skip=0,n=1)
#Yobs = scan(infile,skip=2,n=Nobs)
#obsErr = scan(infile,skip=3,n=Nobs)
#Catch = scan(infile,skip=4,n=Nobs)
#model1.data = list(N=Nobs,I=Yobs,CV=obsErr,C=Catch)


inits <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r=exp(rnorm(1,-0.5,0.2)), k=0.00004, isigma2=100, itau2=100,q=-0.01)
inits2 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r=exp(rnorm(1,-0.5,0.2)), k=0.00003, isigma2=80, itau2=150,q=-0.04)
inits3 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r=exp(rnorm(1,-0.5,0.2)), k=0.00005, isigma2=120, itau2=50,q=-0.006)
inits <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r=exp(rnorm(1,-0.5,0.2)), k=0.00004, isigma2=100, itau2=100,iq=700,idep=0)
inits2 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r=exp(rnorm(1,-0.5,0.2)), k=0.00003, isigma2=80, itau2=150,iq=500,idep=-0.3)
inits3 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r=exp(rnorm(1,-0.5,0.2)), k=0.00005, isigma2=120, itau2=50,iq=400,idep=-0.1)
if (isp ==1 | isp ==4)
 {
inits <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r=exp(rnorm(1,-0.5,0.2)), k=0.00004, isigma2=100, itau2=100,iq=7000,idep=0)
inits2 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r=exp(rnorm(1,-0.5,0.2)), k=0.00003, isigma2=80, itau2=150,iq=5000,idep=-0.3)
inits3 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r=exp(rnorm(1,-0.5,0.2)), k=0.00005, isigma2=120, itau2=50,iq=4000,idep=-0.1)
 }
model1.init <- list(inits,inits2,inits3)


#CHANGES FOR JAGS include
#1. inits as a vector not a list
#2. Truncation syntax in model file different, not I(,) but  T(,)

require(rjags)

model1.jags2save = c("K","r","Sigma2","Tau2","B","FMSY","NextCatch")

#jags <- jags.model("C:/Research/regime/schaef_jags_notatK.txt",
#jags <- jags.model("schaef_jags_notatK.txt",
jags <- jags.model("schaef_jags_notatK_2.txt",
       data = model1.data,
		   inits = model1.init,
                   n.chains = 3,
                   n.adapt = 100)

update(jags,500000)  ####burn-in
mod1 <- jags.samples(jags,model1.jags2save,n.iter=1000000,thin=1000)
dic.mod1 <- dic.samples(jags,n.iter=1000000,thin=1000,type="pD")

#xx is a list containing sampled vectors of desired nodes
mod1.results <- NULL
mod1.results$mod <- mod1
mod1.results$dic <- dic.mod1

#save(mod1.results,file=paste(isp,"_mod1_posterior.Rdata",sep=""))
#save(mod1.results,file=paste(isp,"_mod1_big_posterior.Rdata",sep=""))
#save(mod1.results,file=paste(isp,"_mod1_test_posterior.Rdata",sep=""))
#save(mod1.results,file=paste(isp,"_mod1_20140604_posterior_2.Rdata",sep=""))
#save(mod1.results,file=paste(isp,"_mod1_big_posterior_20140610.Rdata",sep=""))
save(mod1.results,file=paste(isp,"_mod1_big_posterior_20140618.Rdata",sep=""))

#}



#################### 3 REGIMES ##################################


#xx <- jags(model="regime_schaef_jags_notatK.txt",
#           data = model1.data,
#           inits = model2.init,
#           n.chains = 3,
#           DIC=TRUE,n.iter=10000,n.burnin=500,n.thin=10,parameters.to.save=model2.jags2save)

#setwd("C:/Research/regime/")
#require(rjags)
#for (isp in c(2,3:5))
#{
  
  if (isp==1) infile = "gb_cod.dat"
  if (isp==2) infile = "gom_had.dat"
  if (isp==3) infile = "gb_ytf.dat"
  if (isp==4) infile = "sne_winf.dat" 
  if (isp==5) infile = "gom_cod2.dat"

  Nobs = scan(infile,skip=1,n=1)
  Nreg = scan(infile,skip=0,n=1)
  Yobs = scan(infile,skip=2,n=Nobs)
  obsErr = scan(infile,skip=3,n=Nobs)
  Catch = scan(infile,skip=4,n=Nobs)
  model1.data = list(N=Nobs,I=Yobs,CV=obsErr,C=Catch)
  
model3.jags2save = c("K","r","Sigma2","Tau2","TransProbs","Regime","B","RegimeNext","FMSY","NextCatch")

inits <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),r2frac=runif(1,0.2,0.8),r3frac=runif(1,0.2,0.8), k=0.00004, isigma2=100, itau2=100,TransProbs=c(0.97,0.01,0.97,0.01,0.97,0.01),Regime=sample(1:3,Nobs,replace=TRUE),q=-0.01)
inits2 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),r2frac=runif(1,0.2,0.8),r3frac=runif(1,0.2,0.8), k=0.00003, isigma2=80, itau2=150,TransProbs=c(0.7,0.2,0.8,0.1,0.7,0.2),Regime=sample(1:3,Nobs,replace=TRUE),q=-0.04)
inits3 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),r2frac=runif(1,0.2,0.8),r3frac=runif(1,0.2,0.8), k=0.00005, isigma2=120, itau2=50,TransProbs=c(0.9,0.05,0.7,0.2,0.6,0.2),Regime=sample(1:3,Nobs,replace=TRUE),q=-0.006)
inits <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),r2frac=runif(1,0.2,0.8),r3frac=runif(1,0.2,0.8), k=0.00004, isigma2=100, itau2=100,TransProbs=c(0.97,0.01,0.97,0.01,0.97,0.01),Regime=sample(1:3,Nobs,replace=TRUE),iq=500,idep=0)
inits2 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),r2frac=runif(1,0.2,0.8),r3frac=runif(1,0.2,0.8), k=0.00003, isigma2=80, itau2=150,TransProbs=c(0.7,0.2,0.8,0.1,0.7,0.2),Regime=sample(1:3,Nobs,replace=TRUE),iq=300,idep=-0.2)
inits3 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),r2frac=runif(1,0.2,0.8),r3frac=runif(1,0.2,0.8), k=0.00005, isigma2=120, itau2=50,TransProbs=c(0.9,0.05,0.7,0.2,0.6,0.2),Regime=sample(1:3,Nobs,replace=TRUE),iq=600,idep=-0.5)
if (isp ==1 | isp ==4)
 {
inits <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),r2frac=runif(1,0.2,0.8),r3frac=runif(1,0.2,0.8), k=0.00004, isigma2=100, itau2=100,TransProbs=c(0.97,0.01,0.97,0.01,0.97,0.01),Regime=sample(1:3,Nobs,replace=TRUE),iq=5000,idep=0)
inits2 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),r2frac=runif(1,0.2,0.8),r3frac=runif(1,0.2,0.8), k=0.00003, isigma2=80, itau2=150,TransProbs=c(0.7,0.2,0.8,0.1,0.7,0.2),Regime=sample(1:3,Nobs,replace=TRUE),iq=3000,idep=-0.2)
inits3 <- list(P=exp(rnorm(Nobs+1,-0.6,0.2)),r1=exp(rnorm(1,-0.5,0.2)),r2frac=runif(1,0.2,0.8),r3frac=runif(1,0.2,0.8), k=0.00005, isigma2=120, itau2=50,TransProbs=c(0.9,0.05,0.7,0.2,0.6,0.2),Regime=sample(1:3,Nobs,replace=TRUE),iq=6000,idep=-0.5)
 }
model3.init <- list(inits,inits2,inits3)

#jags <- jags.model("C:/Research/regime/schaef_jags_notatK.txt",
#jags <- jags.model("3regime_schaef_jags_notatK.txt",
jags <- jags.model("3regime_schaef_jags_notatK_2.txt",
                   data = model1.data,
                   inits = model3.init,
                   n.chains = 3,
                   n.adapt = 100)

#mod3 <- jags(model="3regime_schaef_jags_notatK.txt",
#           data = model1.data,
#           inits = model3.init,
#           n.chains = 3,
#           DIC=TRUE,n.iter=150000,n.burnin=50000,n.thin=100,parameters.to.save=model2.jags2save)



update(jags,500000)  ####burn-in
mod3 <- jags.samples(jags,model3.jags2save,n.iter=1000000,thin=1000)
dic.mod3 <- dic.samples(jags,n.iter=1000000,thin=1000,type="pD")

#xx is a list containing sampled vectors of desired nodes
mod3.results <- NULL
mod3.results$mod <- mod3
mod3.results$dic <- dic.mod3

##save(mod3.results,file=paste(isp,"_mod3_posterior.Rdata",sep=""))
##save(mod3.results,file=paste(isp,"_mod3_test_posterior.Rdata",sep=""))
#save(mod3.results,file=paste(isp,"_mod3_big_posterior.Rdata",sep=""))
#save(mod3.results,file=paste(isp,"_mod3_20140604_2_posterior.Rdata",sep=""))
#save(mod3.results,file=paste(isp,"_mod3_big_posterior_20140610.Rdata",sep=""))
save(mod3.results,file=paste(isp,"_mod3_big_posterior_20140618.Rdata",sep=""))

}



###################################
