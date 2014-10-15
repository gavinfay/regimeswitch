###########################
#' Code to plot the results of the Markov Switching regime shift models in JAGS
#'
#' @author Gavin Fay


setwd("/media/My\ Passport/NEFSC/research/regime/")
if (Sys.info()[['sysname']] == "Windows") setwd("J:/NEFSC/research/regime/")
if (Sys.info()[['sysname']] == "Darwin")
  setwd("/Volumes/MyPassport/NEFSC/research/regime/")

stock.names <- c("GB cod","GOM haddock","GB Yellowtail","SNE Winter flounder",
                 "GOM cod")

par(mfrow=c(3,2))
for (isp in 1:5)
{
  load(paste(isp,"_mod2_big_posterior.Rdata",sep=""))
  hist(mod2.results$mod$r[2,,1]/mod2.results$mod$r[1,,1],
       main=stock.names[isp],ylab="",xlab="")
}

NextCatch <- matrix(0,nrow=5,ncol=3)
DIC <- matrix(0,nrow=5,ncol=3)
for (isp in 1:5)
{
  if (isp==1) load(paste(isp,"_mod1_big_posterior_20140604.Rdata",sep=""))
  if (isp>1) load(paste(isp,"_mod1_big_posterior_20140610.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod1_big_posterior_20140618.Rdata",sep=""))
  load(paste(isp,"_mod2_big_posterior_20140604.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod2_big_posterior_20140618.Rdata",sep=""))
  if (isp==1) load(paste(isp,"_mod3_big_posterior_20140604.Rdata",sep=""))
  if (isp>1) load(paste(isp,"_mod3_big_posterior_20140610.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod3_big_posterior_20140618.Rdata",sep=""))  
  #if (isp==5) load(paste(isp,"_mod3_20140604_posterior.Rdata",sep=""))
  #load(paste(isp,"_mod3_big_posterior.Rdata",sep=""))
  #if (isp == 1 | isp == 4)
  # {
  #load(paste(isp,"_mod1_20140604_posterior_2.Rdata",sep=""))
  #load(paste(isp,"_mod2_20140604_2_posterior.Rdata",sep=""))
  #load(paste(isp,"_mod3_20140604_2_posterior.Rdata",sep=""))
  # }
  # if (isp==5)  load(paste(isp,"_mod3_20140604_posterior.Rdata",sep=""))
  print(isp)
  DIC[isp,1] <- sum(mod1.results$dic[[1]])+sum(mod1.results$dic[[2]])
  DIC[isp,2] <- sum(mod2.results$dic[[1]])+sum(mod2.results$dic[[2]])
  DIC[isp,3] <- sum(mod3.results$dic[[1]])+sum(mod3.results$dic[[2]])
  NextCatch[isp,1] <- median(mod1.results$mod$NextCatch[1,,1])
  NextCatch[isp,2] <- median(mod2.results$mod$NextCatch[1,,1])
  NextCatch[isp,3] <- median(mod3.results$mod$NextCatch[1,,1])  
}
write.table(cbind(stock.names,DIC),file="DIC.csv",col.names=FALSE,row.names=FALSE,quote=F,sep=",")
write.table(cbind(stock.names,NextCatch),file="NextCatch.csv",col.names=FALSE,row.names=FALSE,quote=F,sep=",")


NextCatch <- matrix(0,nrow=5,ncol=3)
for (isp in 1:5)
{
  load(paste(isp,"_mod1_big_posterior.Rdata",sep=""))
  load(paste(isp,"_mod2_big_posterior.Rdata",sep=""))
  load(paste(isp,"_mod3_big_posterior.Rdata",sep=""))
  if (isp == 1 | isp == 4)
  {
    load(paste(isp,"_mod1_20140604_posterior_2.Rdata",sep=""))
    load(paste(isp,"_mod2_20140604_2_posterior.Rdata",sep=""))
    load(paste(isp,"_mod3_20140604_2_posterior.Rdata",sep=""))
  }
  if (isp==5)  load(paste(isp,"_mod3_20140604_posterior.Rdata",sep=""))
  print(isp)
  NextCatch[isp,1] <- median(mod1.results$mod$NextCatch[1,,1])
  NextCatch[isp,2] <- median(mod2.results$mod$NextCatch[1,,1])
  NextCatch[isp,3] <- median(mod3.results$mod$NextCatch[1,,1])
  if (isp==5) NextCatch[isp,3] <- median(mod3.results$mod$NextCatch[1,,3])
}
write.table(cbind(stock.names,NextCatch),file="NextCatch.csv",col.names=FALSE,row.names=FALSE,quote=F,sep=",")


col.use <- c(
  "#66c2a5",
  "#fc8d62",
  "#8da0cb")


#windows(record=TRUE)
pdf(file="NextCatch.pdf")
par(mfrow=c(3,2),mar=c(2,2,3,1),oma=c(0.5,0.5,0.5,0))
for (isp in 1:5)
{
  #load(paste(isp,"_mod1_big_posterior.Rdata",sep=""))
  #load(paste(isp,"_mod2_big_posterior.Rdata",sep=""))
  #load(paste(isp,"_mod3_big_posterior.Rdata",sep=""))
  #if (isp == 1 | isp == 4)
  # {
  #load(paste(isp,"_mod1_20140604_posterior_2.Rdata",sep=""))
  #load(paste(isp,"_mod2_20140604_2_posterior.Rdata",sep=""))
  #load(paste(isp,"_mod3_20140604_2_posterior.Rdata",sep=""))
  # }
  # if (isp==5)  load(paste(isp,"_mod3_20140604_posterior.Rdata",sep=""))
  
  if (isp==1) load(paste(isp,"_mod1_big_posterior_20140604.Rdata",sep=""))
  if (isp>1) load(paste(isp,"_mod1_big_posterior_20140610.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod1_big_posterior_20140618.Rdata",sep=""))
  load(paste(isp,"_mod2_big_posterior_20140604.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod2_big_posterior_20140618.Rdata",sep=""))
  if (isp==1) load(paste(isp,"_mod3_big_posterior_20140604.Rdata",sep=""))
  if (isp>1) load(paste(isp,"_mod3_big_posterior_20140610.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod3_big_posterior_20140618.Rdata",sep=""))    
  
  ymax <- 1.05*max(c(max(density(mod1.results$mod$NextCatch[1,,1])$y),
  max(density(mod2.results$mod$NextCatch[1,,1])$y),
  max(density(mod3.results$mod$NextCatch[1,,1])$y)))
  
  plot(density(mod1.results$mod$NextCatch[1,,1]),main=stock.names[isp],ylab="",
       xlab="",lwd=3,col=col.use[1],ylim=c(0,ymax),xlim=c(0,
       1.2*quantile(mod1.results$mod$NextCatch[1,,1],probs=0.99)))
  lines(density(mod2.results$mod$NextCatch[1,,1]),lty=1,col=col.use[2],lwd=3)
  #if (isp!=5) 
  lines(density(mod3.results$mod$NextCatch[1,,1]),lty=1,col=col.use[3],lwd=3)
  #if (isp==5) lines(density(mod3.results$mod$NextCatch[1,,3]),lty=3,col="red",lwd=3)
}
plot(0,0,axes=F,col="white",xlab="",ylab="")
legend("center",col=col.use,lwd=5,cex=2,legend=c("No regime change",
      "Two regimes","Three regimes"),bty='n')
dev.off()

pdf(file="NextUMSY.pdf")
par(mfrow=c(3,2),mar=c(2,2,3,1),oma=c(0.5,0.5,0.5,0))
for (isp in 1:5)
{
  #  load(paste(isp,"_mod1_big_posterior.Rdata",sep=""))
  #  load(paste(isp,"_mod2_big_posterior.Rdata",sep=""))
  # load(paste(isp,"_mod3_big_posterior.Rdata",sep=""))
  
  if (isp==1) load(paste(isp,"_mod1_big_posterior_20140604.Rdata",sep=""))
  if (isp>1) load(paste(isp,"_mod1_big_posterior_20140610.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod1_big_posterior_20140618.Rdata",sep=""))
  load(paste(isp,"_mod2_big_posterior_20140604.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod2_big_posterior_20140618.Rdata",sep=""))
  if (isp==1) load(paste(isp,"_mod3_big_posterior_20140604.Rdata",sep=""))
  if (isp>1) load(paste(isp,"_mod3_big_posterior_20140610.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod3_big_posterior_20140618.Rdata",sep=""))
  
  ymax <- 1.05*max(c(max(density(mod1.results$mod$FMSY[1,,1])$y),
                     max(density(mod2.results$mod$FMSY[1,,1])$y),
                     max(density(mod3.results$mod$FMSY[1,,1])$y)))
  
  
  plot(density(mod1.results$mod$FMSY[1,,1]),main=stock.names[isp],ylab="",
       xlab="",lwd=3,col=col.use[1],ylim=c(0,ymax),xlim=c(0,1))
  lines(density(mod2.results$mod$FMSY[1,,1]),lty=1,col=col.use[2],lwd=3,
        xlim=c(0,1))
  lines(density(mod3.results$mod$FMSY[1,,1]),lty=1,col=col.use[3],lwd=3,
        xlim=c(0,1))
}
plot(0,0,axes=F,col="white",xlab="",ylab="")
legend("center",col=col.use,lwd=5,cex=2,legend=c("No regime change",
                                  "Two regimes","Three regimes"),bty='n')
dev.off()


#model fits, not suitable for presentation
#windows(record=TRUE)
for (isp in 1:5)
{
  #load(paste(isp,"_mod2_big_posterior.Rdata",sep=""))
  #load(paste(isp,"_mod2_big_posterior_20140618.Rdata",sep=""))
  load(paste(isp,"_mod2_big_posterior_20140604.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod2_big_posterior_20140618.Rdata",sep=""))
  par(mfrow=c(3,2))
  hist(mod2.results$mod$r[1,,1],xlab="r",ylab="",main="Regime 1 r",xlim=c(0,1.6)) #,add=TRUE,col="white")
  hist(mod2.results$mod$r[2,,1],xlab="r",ylab="",main="Regime 2 r",xlim=c(0,1.6)) #,col=gray(0.8))
  
  #load(paste(isp,"_mod2_20140604_posterior.Rdata",sep=""))
  Bio.summary <- apply(mod2.results$mod$B[,,1],1,quantile,prob=c(0.05,0.25,0.5,0.75,0.95))
  Nyrs <- ncol(Bio.summary)
  plot(Bio.summary[3,],ylim=c(0,max(Bio.summary)),type='l',lwd=3,main="Biomass",axes=F,ylab="",xlab="")
  box()
  axis(1,xlab="Time")
  axis(2,ylab="Biomass (t)")
  polygon(c(1:Nyrs,Nyrs:1),c(Bio.summary[1,1:Nyrs],Bio.summary[5,Nyrs:1]),col=gray(0.95),border=NA)
  polygon(c(1:Nyrs,Nyrs:1),c(Bio.summary[4,1:Nyrs],Bio.summary[2,Nyrs:1]),col=gray(0.8),border=NA)
  lines(Bio.summary[3,1:Nyrs],lwd=3)
  
  if (isp==1) infile = "gb_cod.dat"
  if (isp==2) infile = "gom_had.dat"
  if (isp==3) infile = "gb_ytf.dat"
  if (isp==4) infile = "sne_winf.dat" 
  if (isp==5) infile = "gom_cod2.dat"
  Nobs = scan(infile,skip=1,n=1)
  Nreg = scan(infile,skip=0,n=1)
  Yobs = scan(infile,skip=2,n=Nobs)
  
  points(Yobs/median(mod2.results$mod$q[1,,1]),pch=16,col="red",cex=0.6)
  
  
  #load(paste(isp,"_mod2_big_posterior.Rdata",sep=""))
  Reg.summary <- apply(mod2.results$mod$Regime[,,1],1,quantile,prob=c(0.05,0.25,0.5,0.75,0.95))
  Nyrs <- ncol(Reg.summary)
  plot(Reg.summary[3,],ylim=c(1,2),type='l',lwd=3,main="Regime membership",axes=F)
  box()
  axis(1,xlab="Time")
  axis(2,at=c(1,2),labels=c(1,2))
  polygon(c(1:Nyrs,Nyrs:1),c(Reg.summary[1,1:Nyrs],Reg.summary[5,Nyrs:1]),col=gray(0.95),border=NA)
  polygon(c(1:Nyrs,Nyrs:1),c(Reg.summary[4,1:Nyrs],Reg.summary[2,Nyrs:1]),col=gray(0.8),border=NA)
  #lines(Reg.summary[5,1:Nyrs],lwd=3)
  #lines(Reg.summary[1,1:Nyrs],lwd=3)
  lines(Reg.summary[3,1:Nyrs],lwd=3)
  
  hist(mod2.results$mod$RegimeNext[1,,1],col=gray(0.8),breaks=c(0.5,1.5,2.5),xlab="",ylab="",main="Regime next year",axes=F)
  axis(1,at=c(1,2),labels=c(1,2),cex.lab=1.5)
  axis(2)
  hist(mod2.results$mod$NextCatch[1,,1],xlab="",ylab="",main="Catch next year",breaks=20)
  abline(v=median(mod2.results$mod$NextCatch[1,,1]),lty=2)
}


### Transition probabilities
isp=5
load(paste(isp,"_mod2_big_posterior.Rdata",sep=""))
#windows(record=TRUE)
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(0.5,0.5,0.5,0.5))
hist(1.-mod2.results$mod$TransProbs[1,,1],xlab="",ylab="",main="P(1 -> 2)",
     xlim=c(0,1.0),col=col.use[2]) #,add=TRUE,col="white")
hist(1.-mod2.results$mod$TransProbs[2,,1],xlab="",ylab="",main="P(2 -> 1)",
     xlim=c(0,1.0),col=col.use[1]) #,add=TRUE,col="white")

#### regime plots ####
pdf(file="RegimeMembership.pdf")
for (isp in 1:5)
{
  #isp=5
  #load(paste(isp,"_mod2_big_posterior.Rdata",sep=""))
  #load(paste(isp,"_mod2_big_posterior_20140618.Rdata",sep=""))
  load(paste(isp,"_mod2_big_posterior_20140604.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod2_big_posterior_20140618.Rdata",sep=""))
  par(mfrow=c(2,2),mar=c(3,2,2,2),oma=c(0.5,0.5,0.5,0.5))
  hist(mod2.results$mod$r[1,,1],xlab="r",ylab="",main="High Regime r",
       col=col.use[1],xlim=c(0,1.6)) #,add=TRUE,col="white")
  hist(mod2.results$mod$r[2,,1],xlab="r",ylab="",main="Low Regime r",
       col=col.use[2],xlim=c(0,1.6)) #,col=gray(0.8))
  
  
  #load(paste(isp,"_mod2_big_posterior.Rdata",sep=""))
  Reg.summary <- apply(3+-1*mod2.results$mod$Regime[,,1],1,quantile,
                       prob=c(0.05,0.25,0.5,0.75,0.95))
  Nyrs <- ncol(Reg.summary)
  plot(Reg.summary[3,],ylim=c(1,2),type='l',lwd=3,main="Regime membership",axes=F)
  box()
  axis(1,xlab="Year")
  par(las=2)
  axis(2,at=c(1,2),labels=c("lo","hi"))
  par(las=0)  
  polygon(c(1:Nyrs,Nyrs:1),c(Reg.summary[1,1:Nyrs],Reg.summary[5,Nyrs:1]),
          col=gray(0.95),border=NA)
  polygon(c(1:Nyrs,Nyrs:1),c(Reg.summary[4,1:Nyrs],Reg.summary[2,Nyrs:1]),
          col=gray(0.8),border=NA)
  #lines(Reg.summary[5,1:Nyrs],lwd=3)
  #lines(Reg.summary[1,1:Nyrs],lwd=3)
  lines(Reg.summary[3,1:Nyrs],lwd=3)
  
  hist(mod2.results$mod$RegimeNext[1,,1],col=col.use[1:2],breaks=c(0.5,1.5,2.5),
       xlab="",ylab="",main="Regime next year",axes=F)
  axis(1,at=c(1,2),labels=c("high r","low r"),cex.lab=1.5)
  axis(2)
  #hist(mod2.results$mod$NextCatch[1,,1],xlab="",ylab="",main="Catch next year",breaks=20)
  #abline(v=median(mod2.results$mod$NextCatch[1,,1]),lty=2)
}
dev.off()

pdf(file="BiomassFits.pdf")
for (isp in 1:5)
{
  
  if (isp==1) load(paste(isp,"_mod1_big_posterior_20140604.Rdata",sep=""))
  if (isp>1) load(paste(isp,"_mod1_big_posterior_20140610.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod1_big_posterior_20140618.Rdata",sep=""))
  load(paste(isp,"_mod2_big_posterior_20140604.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod2_big_posterior_20140618.Rdata",sep=""))
  if (isp==1) load(paste(isp,"_mod3_big_posterior_20140604.Rdata",sep=""))
  if (isp>1) load(paste(isp,"_mod3_big_posterior_20140610.Rdata",sep=""))
  if (isp==5) load(paste(isp,"_mod3_big_posterior_20140618.Rdata",sep=""))  
  
  par(mfrow=c(1,1),mar=c(3,3,2,2),oma=c(2,2,0.5,0.5))
#load(paste(isp,"_mod2_20140604_posterior.Rdata",sep=""))
Bio.summary <- apply(mod2.results$mod$B[,,1],1,quantile,
                     prob=c(0.05,0.25,0.5,0.75,0.95))
Nyrs <- ncol(Bio.summary)
plot(Bio.summary[3,],ylim=c(0,max(Bio.summary)),type='l',lwd=3,
     axes=F,ylab="",xlab="")
box()
axis(1)
mtext("Year",side=1,line=3,cex=1.2)
axis(2)
mtext("Biomass (t)",side=2,line=3,cex=1.2)
polygon(c(1:Nyrs,Nyrs:1),c(Bio.summary[1,1:Nyrs],Bio.summary[5,Nyrs:1]),
        col=gray(0.95),border=NA)
polygon(c(1:Nyrs,Nyrs:1),c(Bio.summary[4,1:Nyrs],Bio.summary[2,Nyrs:1]),
        col=gray(0.8),border=NA)
lines(Bio.summary[3,1:Nyrs],lwd=3,col=col.use[2])

if (isp==1) infile = "gb_cod.dat"
if (isp==2) infile = "gom_had.dat"
if (isp==3) infile = "gb_ytf.dat"
if (isp==4) infile = "sne_winf.dat" 
if (isp==5) infile = "gom_cod2.dat"
Nobs = scan(infile,skip=1,n=1)
Nreg = scan(infile,skip=0,n=1)
Yobs = scan(infile,skip=2,n=Nobs)



Bio.summary <- apply(mod1.results$mod$B[,,1],1,quantile,
                     prob=c(0.05,0.25,0.5,0.75,0.95))
lines(Bio.summary[3,1:Nyrs],lwd=3,col=col.use[1])
Bio.summary <- apply(mod3.results$mod$B[,,1],1,quantile,
                     prob=c(0.05,0.25,0.5,0.75,0.95))
lines(Bio.summary[3,1:Nyrs],lwd=3,col=col.use[3])

points(Yobs/median(mod2.results$mod$q[1,,1]),pch=16,col="black",cex=0.9)

Bio.summary <- apply(mod2.results$mod$B[,,1],1,quantile,
                     prob=c(0.05,0.25,0.5,0.75,0.95))
lines(Bio.summary[3,1:Nyrs],lwd=3,col=col.use[2])

legend("topright",col=col.use,lwd=5,cex=1.2,
       legend=c("No regime change","Two regimes","Three regimes"),bty='n')

}
dev.off()
