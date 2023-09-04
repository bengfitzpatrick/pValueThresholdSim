library(dplyr)
library("tictoc")
tic()
myFileName  <- "testFile.Rdata"
N           <- 5000
nMC         <- 100
e0          <- 20
eta         <- 0.01; #0.02; #0.3
LM          <- 1.0
sigma_d_t   <- 0.0;
sigma_d     <- 0.05
sigma_r     <- 2
sigma_e     <- 10
num_labs    <- 100
lab_id      <- 1:num_labs
d_sample    <- 10
tau   		<- 0.8
E     		<- 0.2
EP          <- 0.2

alfa 		<- 0.05
V_novel 	<- 1
V_neg       <- 0.5

pEff        <- 0.0
maxEff      <- 500
minEff      <- 5
rrRate      <- 0.05
effRate     <- 0.05
deathRate   <- 0.005
delayTime   <- 100

hOrR        <- 1

replication_rate  <- 1
maxRR             <- 1

testMat     <- matrix(0,ncol=num_labs,nrow=maxRR)
deltMat     <- testMat;

labs           <- data.frame(matrix(ncol=10, nrow= num_labs))
colnames(labs) <- c('age','id','effort','power','replication_rate','value','alpha','runExp','cTP','cFP')

powerList  <- matrix(0,ncol=3,nrow=N)
effortList <- matrix(0,ncol=3,nrow=N)
alfaList   <- matrix(0,ncol=3,nrow=N)
valueList  <- matrix(0,ncol=3,nrow=N)
rrList     <- matrix(0,ncol=3,nrow=N)
mutList    <- matrix(0,ncol=3,nrow=N)

indxMat    <- rep(seq(1,maxRR),times=num_labs);
indxMat    <- t(matrix(indxMat,ncol=num_labs))

finalRR   <- matrix(0,nrow=num_labs,ncol=nMC)
finalEff  <- matrix(0,nrow=num_labs,ncol=nMC)
finalVal  <- matrix(0,nrow=num_labs,ncol=nMC)
finalFP   <- matrix(0,nrow=num_labs,ncol=nMC)
finalTP   <- matrix(0,nrow=num_labs,ncol=nMC)
finalAge  <- matrix(0,nrow=num_labs,ncol=nMC)
finalDelt <- matrix(0,nrow=num_labs,ncol=nMC)
finalcFP  <- matrix(0,nrow=num_labs,ncol=nMC)
finalcTP  <- matrix(0,nrow=num_labs,ncol=nMC)

for (iMC in 1:nMC){

    print(paste("iMC = ",iMC,sep=""))
	labs$age               <- 0
	labs$id                <- seq(1,num_labs)
	labs$effort            <- e0
	labs$power             <- 0
	labs$replication_rate  <- replication_rate
	labs$value             <- 0
	labs$alpha             <- 0
	labs$runExp            <- 0
	labs$cTP               <- 0
	labs$cFP               <- 0
	
	#myDelta                <- matrix(rgamma(num_labs*maxRR,1,rate = -log(1-tau)/E),ncol=maxRR)
	delta1                 <- rgamma(num_labs,1,rate = -log(1-tau)/EP)
	#delta1                 <- rgamma(num_labs,2,scale=EP)
	myDelta                <- matrix(rep(delta1,times=maxRR),ncol=maxRR)

	for (i in 1:N) {

	  h       <- pEff+(1-pEff)*exp(-(labs$effort-minEff)*eta)
	  u       <- runif(num_labs,0,1) 
	  w       <- runif(num_labs,0,1)
	  x       <- runif(num_labs,0,1)
	  deaths  <- runif(num_labs,0,1)
	  rrs     <- runif(num_labs,0,1)
	  effs    <- runif(num_labs,0,1);
	  
	  runExp       <- ifelse(h>=w,1,0);
	  labs$runExp  <- labs$runExp+runExp;
	  
	  alfNow   <- 1-(1-alfa)^(1/labs$replication_rate)
	  t_alfNow <- qt(1-alfNow/2, 2*labs$effort-2)
	  t_alfa   <- qt(1-alfa/2, 2*labs$effort-2)

	  
	  detct   <- rep(0,num_labs);
	  truth   <- rep(0,num_labs);
	  detCorr <- rep(0,num_labs);
	  
	  rrMat    <- matrix(rep(labs$replication_rate,times=maxRR),ncol=maxRR);
	  testMat  <- ifelse(indxMat<=rrMat,1,0)
	  
	  myDelt_t <- matrix(rnorm(maxRR*num_labs,0,sigma_d_t),ncol=maxRR);
	  myDelt_t <- ifelse(myDelt_t>=0,myDelt_t,0)
	  
	  if (hOrR == 1) {
	  
		myDelta2 <- (myDelta+myDelt_t)*testMat
		
	  }else{
	  	#delta1                 <- rgamma(num_labs,1,rate = -log(1-tau)/EP)
	    #myDelta3               <- matrix(rep(delta1,times=maxRR),ncol=maxRR)
		myDelta3               <- matrix(rgamma(num_labs*maxRR,1,rate = -log(1-tau)/E),ncol=maxRR)
		myDelta2               <- (myDelta3)*testMat
		
	 }
	  
	  myTruth  <- ifelse(myDelta2<=E,0,1)
	  
	  effMat   <- matrix(rep(labs$effort,times=maxRR),ncol=maxRR);
	  myNcp    <- sqrt(effMat)*myDelta2/sqrt(2)
	  tVals    <- abs(rt(rrMat,2*effMat-2,myNcp)*testMat)
	  tMax     <- apply(tVals,1,max)
	  
	  tMaxMat  <- matrix(rep(tMax,times=maxRR),ncol=maxRR);
	  chkMax   <- ifelse(tMaxMat == tVals,1,0);
	  detct    <- tMax>t_alfa
	  detCorr  <- tMax>t_alfNow
	  truth    <- ifelse(rowSums(chkMax*myTruth)>0,1,0)  
	  
	  truePos <- ifelse(((truth==1)&(detCorr==1)&(h>=w)),1,0)
	  falsPos <- ifelse(((truth==0)&(detct==1)&(h>=w)),1,0)
	  
	  labs$alpha <- falsPos
	  labs$power <- truePos
	  labs$value <- labs$value + (truePos + falsPos)*V_novel
	  labs$age   <- labs$age + 1
	  labs$cTP   <- labs$cTP + truePos
	  labs$cFP   <- labs$cFP + falsPos 
	  
	  if(i>delayTime){
	  
	  indxDeath   <- which(deaths<=deathRate);
	  indxRR      <- which(rrs<=rrRate);
	  indxEff     <- which(effs<=effRate);
	  
	  mutList[i,1] <- length(indxDeath)
	  mutList[i,2] <- length(indxRR)
	  mutList[i,3] <- length(indxEff)
	  
	  L           <- length(indxRR);
	  L2          <- min(floor(LM*L),num_labs);
	 
	  bestLabs    <- sort.int(labs$value,index.return=TRUE,decreasing=TRUE)
	  
	  indxRep     <- sample(bestLabs$ix[1:L2],L,replace=FALSE)
	  
	  labs$replication_rate[indxRR] <- round(labs$replication_rate[indxRep] + rnorm(L,0,sigma_r))
	  labs$replication_rate         <- ifelse(labs$replication_rate>maxRR,maxRR,labs$replication_rate)
	  labs$replication_rate         <- ifelse(labs$replication_rate<1,1,labs$replication_rate)
	  
	  L           <- length(indxEff);
	  
	  L2          <- min(floor(LM*L),num_labs);
	  
	  indxRep     <- sample(bestLabs$ix[1:L2],L,replace=FALSE)
	  
	  labs$effort[indxEff]  <- round(labs$effort[indxRep] + rnorm(L,0,sigma_e))
	  labs$effort           <- ifelse(labs$effort>maxEff,maxEff,labs$effort)
	  labs$effort           <- ifelse(labs$effort<minEff,minEff,labs$effort) 
	  

	  L           <- length(indxDeath);
	  L2          <- min(floor(LM*L),num_labs);
	  
	  indxRep     <- sample(bestLabs$ix[1:L2],L,replace=FALSE)
	  
	  nDeath                 <- length(indxDeath)
	  labs[indxDeath,]       <- labs[indxRep,]
	  labs$age[indxDeath]    <-1
	  labs$runExp[indxDeath] <-0
	  labs$value[indxDeath]  <-0
	  labs$cFP[indxDeath]    <-0
	  labs$cTP[indxDeath]    <-0
	  
	  labs$id[indxDeath]     <- seq(max(labs$id)+1,max(labs$id)+L)
	  labs$effort[indxDeath] <- round(labs$effort[indxDeath] + rnorm(L,0,sigma_e))
	  labs$effort            <- ifelse(labs$effort>maxEff,maxEff,labs$effort)
	  labs$effort            <- ifelse(labs$effort<minEff,minEff,labs$effort)
	  
	  labs$replication_rate[indxDeath] <- round(labs$replication_rate[indxDeath] + rnorm(L,0,sigma_r))
	  labs$replication_rate            <- ifelse(labs$replication_rate>maxRR,maxRR,labs$replication_rate)
	  labs$replication_rate            <- ifelse(labs$replication_rate<1,1,labs$replication_rate)  
	  
	  # myDelta[indxDeath,]           <- myDelta[indxRep,] + matrix(rnorm((nDeath*maxRR),0,sigma_d),ncol=maxRR)
	  # myDelta                       <- ifelse(myDelta<1e-6,1e-6,myDelta)
	  }
	  
	  powerList[i,]  <- powerList[i,] + c(sum(labs$power)/num_labs-1/sqrt(num_labs),sum(labs$power)/num_labs,sum(labs$power)/num_labs+1/sqrt(num_labs))/nMC  
	  alfaList[i,]   <- alfaList[i,] + c(sum(labs$alpha)/num_labs-1/sqrt(num_labs),sum(labs$alpha)/num_labs,sum(labs$alpha)/num_labs+1/sqrt(num_labs))/nMC     
	  effortList[i,] <- effortList[i,] + quantile(labs$effort,probs = c(0.25,0.5,0.75))/nMC
	  valueList[i,]  <- valueList[i,] + quantile(labs$value,probs = c(0.25,0.5,0.75))/nMC
	  rrList[i,]     <- rrList[i,] + quantile(labs$replication_rate,probs = c(0.25,0.5,0.75))/nMC

	}
	medDelt         <- 0*labs$effort
	for (iii in 1:num_labs){
		 medDelt[iii]<- median(myDelta[iii,(1:labs$replication_rate[iii])])
	}
	finalRR[,iMC]   <- labs$replication_rate
	finalEff[,iMC]  <- labs$effort
	finalVal[,iMC]  <- labs$value
	finalFP[,iMC]   <- labs$alpha
	finalTP[,iMC]   <- labs$power
	finalAge[,iMC]  <- labs$age
	finalDelt[,iMC] <- medDelt
	finalcFP[,iMC]  <- labs$cFP
	finalcTP[,iMC]  <- labs$cTP
}

print(labs)

#print(library)
dev.new()
plot(powerList[,2],ylim=c(0,1),xlab="time step",ylab="True Discovery Rate")
lines(powerList[,1],col="red")
lines(powerList[,3],col="red")
dev.new()
plot(alfaList[,2],ylim=c(0,1),xlab="time step",ylab="False Discovery Rate")
lines(alfaList[,1],col="red")
lines(alfaList[,3],col="red")
dev.new()
plot(effortList[,2],xlab="time step",ylab="Effort")
lines(effortList[,1],col="red")
lines(effortList[,3],col="red")
dev.new()
plot(valueList[,2],xlab="time step",ylab="Value")
lines(valueList[,1],col="red")
lines(valueList[,3],col="red")
dev.new()
plot(rrList[,2],ylim=c(0,20),xlab="time step",ylab="re-try rate")
lines(rrList[,1],col="red")
lines(rrList[,3],col="red")
toc()
save.image(file=myFileName)
