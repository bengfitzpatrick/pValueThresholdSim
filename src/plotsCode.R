load("D:\\tempest\\projects\\integrity\\code\\rerunData\\ic_r_a005rm40.RData")

pPowerList <- powerList
pAlfaList  <- alfaList
pValueList <- valueList
pEffortList<- effortList
pRetryList <- rrList
pFPList    <- alfaList[,2]/(powerList[,2]+alfaList[,2])

yrs        <- seq(1,N)

print(30*mean(labs$value/labs$age))
print(hOrR)
print(num_labs)

finalVal005  <- matrix(finalVal,ncol=1)
finalEff005  <- matrix(finalEff,ncol=1)
finalRR005   <- matrix(finalRR,ncol=1)
finalDelt005 <- matrix(finalDelt,ncol=1)
finalAge005  <- matrix(finalDelt,ncol=1)
finalFP005   <-colSums(finalcFP)
finalTP005   <-colSums(finalcTP)
finalFPP005  <- finalFP005/(finalFP005+finalTP005)

load("D:\\tempest\\projects\\integrity\\code\\rerunData\\ic_r_a05rm1.RData")

sPowerList <- powerList
sAlfaList  <- alfaList
sValueList <- valueList
sEffortList<- effortList
sRetryList <- rrList
sFPList    <- alfaList[,2]/(powerList[,2]+alfaList[,2])

print(30*mean(labs$value/labs$age))
print(hOrR)
print(num_labs)

finalVal1  <- matrix(finalVal,ncol=1)
finalEff1  <- matrix(finalEff,ncol=1)
finalRR1   <- matrix(finalRR,ncol=1)
finalDelt1 <- matrix(finalDelt,ncol=1)
finalAge1  <- matrix(finalDelt,ncol=1)
finalFP1   <-colSums(finalcFP)
finalTP1   <-colSums(finalcTP)
finalFPP1  <- finalFP1/(finalFP1+finalTP1)

load("D:\\tempest\\projects\\integrity\\code\\rerunData\\ic_r_a05rm40.RData")

rPowerList <- powerList
rAlfaList  <- alfaList
rValueList <- valueList
rEffortList<- effortList
rRetryList <- rrList
rFPList    <- alfaList[,2]/(powerList[,2]+alfaList[,2])

print(30*mean(labs$value/labs$age))
print(hOrR)
print(num_labs)

finalVal05  <- matrix(finalVal,ncol=1)
finalEff05  <- matrix(finalEff,ncol=1)
finalRR05   <- matrix(finalRR,ncol=1)
finalDelt05 <- matrix(finalDelt,ncol=1)
finalAge05  <- matrix(finalDelt,ncol=1)
finalFP05   <-colSums(finalcFP)
finalTP05   <-colSums(finalcTP)
finalFPP05  <- finalFP05/(finalFP05+finalTP05)

dev.new()
plot(yrs,sFPList,ylim=c(0,1),xlab="steps",ylab="Fraction of *False Pos* Pubs")

points(yrs,rFPList,col="blue")
points(yrs,pFPList,col="red")
legend(x="topleft",legend=c("1 hypothesis","many hypotheses","p<0.005, many hypotheses"),lty=c(1,1,1),col=c("black","blue","red"),lwd=3)
dev.new()
plot(yrs,sAlfaList[,2],ylim=c(0,1),xlab="steps",ylab="Fraction of *False Pos* Results")

points(yrs,rAlfaList[,2],col="blue")
points(yrs,pAlfaList[,2],col="red")
legend(x="topleft",legend=c("1 hypothesis","many hypotheses","p<0.005, many hypotheses"),lty=c(1,1,1),col=c("black","blue","red"),lwd=3)
dev.new()
plot(yrs,sEffortList[,2],xlab="steps",ylab="Effort",ylim=c(0,100))

points(yrs,rEffortList[,2],col="blue")
points(yrs,pEffortList[,2],col="red")
legend(x="topleft",legend=c("1 hypothesis","many hypotheses","p<0.005, many hypotheses"),lty=c(1,1,1),col=c("black","blue","red"),lwd=3)
dev.new()
plot(yrs,sValueList[,2],xlab="steps",ylab="Value",ylim=c(0,400))

points(yrs,rValueList[,2],col="blue")
points(yrs,pValueList[,2],col="red")
legend(x="topleft",legend=c("1 hypothesis","many hypotheses","p<0.005, many hypotheses"),lty=c(1,1,1),col=c("black","blue","red"),lwd=3)
dev.new()
plot(yrs,sRetryList[,2],ylim=c(0,40),xlab="steps",ylab="re-try rate")

points(yrs,rRetryList[,2],col="blue")
points(yrs,pRetryList[,2],col="red")
legend(x="topleft",legend=c("1 hypothesis","many hypotheses","p<0.005, many hypotheses"),lty=c(1,1,1),col=c("black","blue","red"),lwd=3)

dev.new()
par(mfrow=c(2,2))
plot(yrs,sEffortList[,2],xlab="steps",ylab="Effort",ylim=c(0,100))
points(yrs,rEffortList[,2],col="blue")
points(yrs,pEffortList[,2],col="red")


plot(yrs,sValueList[,2],xlab="steps",ylab="Value",ylim=c(0,400))
points(yrs,rValueList[,2],col="blue")
points(yrs,pValueList[,2],col="red")


plot(yrs,sRetryList[,2],ylim=c(0,40),xlab="steps",ylab="Number of hypotheses")
points(yrs,rRetryList[,2],col="blue")
points(yrs,pRetryList[,2],col="red")
legend(x="topleft",legend=c("1 hypothesis","many hypotheses","p<0.005, many hypotheses"),lty=c(1,1,1),col=c("black","blue","red"),lwd=3)

plot(yrs,sFPList,ylim=c(0,1),xlab="steps",ylab="Fraction of *False Pos* Pubs")
points(yrs,rFPList,col="blue")
points(yrs,pFPList,col="red")

dev.new()
par(mfrow=c(3,1))
hist(finalEff1,breaks=seq(5,500,by=5),main="Effort histogram",xlab="Effort",ylab="Frequency",xlim=c(0,100))
hist(finalEff05,breaks=seq(5,500,by=5),col="blue",main="Effort histogram",xlab="Effort",ylab="Frequency",xlim=c(0,100))
hist(finalEff005,breaks=seq(5,500,by=5),col="red",main="Effort histogram",xlab="Effort",ylab="Frequency",xlim=c(0,100))

dev.new()
par(mfrow=c(3,1))
hist(finalVal1,breaks=seq(0,4400,by=20),main="Value histogram",xlab="Value",ylab="Frequency",xlim=c(0,800))
hist(finalVal05,col="blue",breaks=seq(0,4400,by=20),main="Value histogram",xlab="Value",ylab="Frequency",xlim=c(0,800))
hist(finalVal005,col="red",breaks=seq(0,4400,by=20),main="Value histogram",xlab="Value",ylab="Frequency",xlim=c(0,800))

dev.new()
par(mfrow=c(3,1))
hist(finalRR1,breaks=seq(0,80,by=1),main="Number of hypotheses histogram",xlab="Number of hypotheses",ylab="Frequency",xlim=c(0,20))
hist(finalRR05,col="blue",breaks=seq(0,80,by=1),main="Number of hypotheses histogram",xlab="Number of hypotheses",ylab="Frequency",xlim=c(0,20))
hist(finalRR005,col="red",breaks=seq(0,80,by=1),main="Number of hypotheses histogram",xlab="Number of hypotheses",ylab="Frequency",xlim=c(0,20))

dev.new()
par(mfrow=c(3,1))
hist(finalFPP1,breaks=seq(0,1,by=0.02),main="False positive rate histogram",xlab="False positive rate",ylab="Frequency",xlim=c(0,1))
hist(finalFPP05,breaks=seq(0,1,by=0.02),col="blue",main="False positive rate histogram",xlab="False positive rate",ylab="Frequency",xlim=c(0,1))
hist(finalFPP005,breaks=seq(0,1,by=0.02),col="red",main="False positive rate histogram",xlab="False positive rate",ylab="Frequency",xlim=c(0,1))