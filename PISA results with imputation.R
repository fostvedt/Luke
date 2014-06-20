require(data.table)
require(quantreg)
require(reshape)
require(doBy)
require(mice)
require(lme4)
require(ggplot2)
setwd("/Users/lukefostvedt/Documents/PISA 2012/Luke QR/")
source("/Users/lukefostvedt/Documents/Journal Articles/quantile regression/Code from papers/Geraci MI code/Quantile MI functions.R")
load("USdata.Rdata")

scq <- read.csv("/Users/lukefostvedt/Documents/PISA 2012/CSV 2012 files/scq.csv")
schUS <- scq[which(scq$CNT=="United States of America"),]

#m <- match(paste("PV",1:5,"MATH", sep=""),names(stuUS))
#mCC <- match(paste("PV",1:5,"MACC", sep=""),names(stuUS))
#mCQ <- match(paste("PV",1:5,"MACQ", sep=""),names(stuUS))
#mCS <- match(paste("PV",1:5,"MACS", sep=""),names(stuUS))
#mCU <- match(paste("PV",1:5,"MACU", sep=""),names(stuUS))
#mPE <- match(paste("PV",1:5,"MAPE", sep=""),names(stuUS))
#mPF <- match(paste("PV",1:5,"MAPF", sep=""),names(stuUS))
#mPI <- match(paste("PV",1:5,"MAPI", sep=""),names(stuUS))
#s <- match(paste("PV",1:5,"SCIE", sep=""),names(stuUS))
#r <- match(paste("PV",1:5,"READ", sep=""),names(stuUS))
#
ind <- match(stuUS$SCHOOLID,schUS$SCHOOLID)
stuU <- cbind(stuUS,schUS[ind,])
cols <- c(1:7,405:550,631:644,900:935)
USdat <- stuU[,cols]




unique(USdat$FAMSTRUC)
apply(USdat, 2, function(x) sum(x=="M")) -> Mis
USdat2 <- USdat[,which(Mis < 3000)]

num <- c("HEDRES", "hisei", "HOMEPOS", "LMINS", "SMINS", "MMINS", "PARED", "WEALTH", "CULTPOS", "ESCS", "GRADE", "BFMJ2", "BMMJ1", names(schUS[,c(258:265, 267:275,278:290)]), names(stuU[,c(500:550,631:644)]),"W_FSTUWT")

fac <- c("CNT","STRATUM","SCHOOLID","TestLANG","LANGN","IMMIG","REPEAT","SCHLTYPE" )

for (j in num){
	USdat2[,j] <- as.numeric(USdat2[,j])
}

for (j in fac){
	ind <- which(USdat2[,j] == "M")
	if(length(ind>0)) USdat2[ind,j] <- NA
	USdat2[,j] <- as.factor(USdat2[,j])
}

summary(USdat2[,num])

USdat1 <- cbind(USdat2[,num], USdat2[,fac])
apply(USdat1, 2, function(x) sum(is.na(x))) -> Mis
m <- names(which(Mis>4000))
mm <- match(m,names(USdat1))
USdat1 <- USdat1[,-mm]
x <- USdat1


#USdat <- cbind(stuUS[,num], stuUS[,fac])

#var <- c("pvM","pvMCC","pvMCQ","pvMCS","pvMCU", "pvMPE", "pvMPF", "pvMPE", "pvMPI", "pvS","pvR","HEDRES", "HOMEPOS", "PARED", "WEALTH", "CULTPOS", "ESCS", "GRADE", "IMMIG", "REPEAT","SCHOOLID", names(schUS[,c(258:265,267:275,278:290)]))





hdat <- USdat1[,c("ESCS","HEDRES","CULTPOS","PARED", "WEALTH","SMINS","LMINS","MMINS","W_FSTUWT","W_FSCHWT",paste("PV",1:5,"MATH",sep=""),paste("PV",1:5,"SCIE",sep=""),paste("PV",1:5,"READ",sep=""))]

apply(hdat,1,function(s) sum(is.na(s))) -> miss
dim(hdat[as.numeric(which(miss>0)),])
dat <- hdat[-which(miss>0),]

model1 <- formula(ESCS~HEDRES+ CULTPOS+PARED+WEALTH)
model2S <- formula(y~Y + HEDRES  + CULTPOS+WEALTH  + Ynuhat)
model2M <- formula(y~Y + HEDRES  + CULTPOS+WEALTH  + Ynuhat)
model2R <- formula(y~Y + HEDRES  + CULTPOS+WEALTH  + Ynuhat)







taus <- c(.10,.25,.50,.75,.90)
tau <- taus[5]
npar = 6
vc1 <- array(0,dim=c(npar,npar,19))
betam<- array(0,dim=c(19,npar))
test <- "MATH"

for(mm in 1:5){ 
	data1 <- cbind(complete(imp,mm), hdat[,paste("PV",mm,test,sep="")], hdat[,c("W_FSCHWT","W_FSTUWT")])
	data1$y <- hdat[,paste("PV",mm,test,sep="")]
	data1$Y <- data1$ESCS
	m1 <- rq(model1, tau=tau, data=data1,weights=W_FSTUWT) 
	dat <- data1
	dat$Ystar <- predict(m1)
	dat$nuhat=dat$Y-dat$Ystar
	dat$Ynuhat=dat$Y*dat$nuhat
for(jj in 1:19){
mS <- rq(model2S, tau=jj/20, data=dat,weights=W_FSTUWT) 
betam[jj,] <- betam[jj,] + mS$coefficients
vc1[,,jj] <- vc1[,,jj] + summary.rq(mS,se="ker",covariance=T)$cov # produces a covariance matrix
}
}
vc1/5
betam/5


vc[[i]] <- vc1
}

vlistS <- list()
vlistM <- list()
vlistR <- list()
vcovS <- array(0,dim=c(19,5,5))
vcovM <- array(0,dim=c(19,5,5))
vcovR <- array(0,dim=c(19,5,5))



#var1 <- c("HEDRES", "HOMEPOS", "PARED", "WEALTH", "CULTPOS", "ESCS","CLSIZE")
#var1 <- c("HEDRES", "HOMEPOS", "PARED", "WEALTH", "CULTPOS", "ESCS","LMINS","MMINS","SMINS","CLSIZE")
#impdat <- USdat1[,var1]	
#method <- c("rq","rq","rq","rq","rq","rq","")
#imp <- mice(impdat,m=5,maxit=5,method=method)
#save(imp,file="m5imputation.Rdata")
#load(file="m5imputation.Rdata")


m <- 5
for(i in 1:m){	
var <- c("SCHOOLID", paste("PV",i,"SCIE",sep=""), paste("PV",i,"MATH",sep=""), paste("PV",i,"READ",sep=""),"CLSIZE","W_FSCHWT","W_FSTUWT")
dat <- cbind(USdat1[,var],complete(imp,1))
dat$ym <- dat[,paste("PV",i,"MATH",sep="")]
dat$ys <- dat[,paste("PV",i,"SCIE",sep="")]
dat$yr <- dat[,paste("PV",i,"READ",sep="")]

taus <- c(.10,.25,.50,.75,.90)

for(j in taus){
dat$Y <- dat$ESCS
m1 <- rq(model1, tau=j, data=dat,weights=W_FSCHWT) 
dat$Ystar <- predict(m1)
dat$nuhat=dat$Y-dat$Ystar
dat$Ynuhat=dat$Y*dat$nuhat

for(k in 1:19){
mS <- rq(model2S, tau=k/20, data=dat,weights=W_FSTUWT) 
mM <- rq(model2M, tau=k/20, data=dat,weights=W_FSTUWT) 
mR <- rq(model2R, tau=k/20, data=dat,weights=W_FSTUWT) 
vlistS[[k]] <- summary.rq(mS,se="ker",covariance=T)$cov
vlistM[[k]] <- summary.rq(mM,se="ker",covariance=T)$cov
vlistR[[k]] <- summary.rq(mR,se="ker",covariance=T)$cov

vcovS[[k]] <- vlistS[[k]] + vcovS[[k]]
vcovM[[k]] <- vlistM[[k]] + vcovM[[k]]
vcovR[[k]] <- vlistR[[k]] + vcovR[[k]]
}
}
}

#titleS <- paste('Sciencep=',taus*100,sep="")
#titleM <- paste('Mathp=',taus*100,sep="")
#titleR <- paste('Readingp=',taus*100,sep="")
#
#nameS <- paste(titleS,".pdf",sep="")
#nameM <- paste(titleM,".pdf",sep="")
#nameR <- paste(titleR,".pdf",sep="")

#pdf(nameS[i])
#plot(summary(mS),parm=2,main=titleS[i],cex=2,lwd=3,ylim=c(0,60));dev.off()
#pdf(nameM[i])
#plot(summary(mM),parm=2,main=titleM[i],cex=2,lwd=3,ylim=c(0,60));dev.off()
#pdf(nameR[i])
#plot(summary(mR),parm=2,main=titleR[i],cex=2,lwd=3,ylim=c(0,60));dev.off()

qplot(WEALTH,pvM,data=x)

summary(lmer(pvM~pvR+pvS+HOMEPOS+ ESCS+ GRADE + (1|SCHOOLID),data=x)->mlm)





