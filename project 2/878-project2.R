################################# problem 1 radon question
#read in data
radon <- read.csv(file="/Users/nannan/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/2022 Spring/STAT-878/Project/project 1/radon.csv")
head(radon)
radon <- radon$Radon

# 1.a create plot for acf and pacf
dev.new(width=8,height=6,pointsize=10)
par(mfcol=c(1,2))
acf(x=radon,lag.max = 20,ylim=c(-1,1),main="ACF of Radon time series")
pacf(x=radon,lag.max = 20,ylim=c(-1,1),main="PACF of Radon time series")

# 1.d plot true ACF and PACF plots of AR(1) with phi1=0.3013
dev.new(width=8,height=6,pointsize=10)
par(mfcol=c(1,2))
plot(y=ARMAacf(ar=0.3013,lag.max = 20),x=0:20,type = "h",ylim = c(-1,1),xlab = "h",ylab = expression(rho(h)),main=expression(paste("TRUE ACF of AR(1) with",phi1[1]==0.3013)))
abline(h=0)
plot(x=ARMAacf(ar=0.3013,lag.max = 20,pacf=TRUE),type = "h",ylim = c(-1,1),xlab = "h",ylab=expression(phi1[hh]),main=expression(paste("TRUE PACF of AR(1) with",phi1[1]==0.3013)))
abline(h=0)

# 1.e plot true ACF and PACF plots of MA(1) with theta1=0.2592
dev.new(width=8,height=6,pointsize=10)
par(mfcol=c(1,2))
plot(y=ARMAacf(ma=0.2592,lag.max = 20),x=0:20,type = "h",ylim = c(-1,1),xlab = "h",ylab = expression(rho(h)),main=expression(paste("TRUE ACF of MA(1) with",theta[1]==0.2592)))
abline(h=0)
plot(x=ARMAacf(ma=0.2592,lag.max = 20,pacf=TRUE),type = "h",ylim = c(-1,1),xlab = "h",ylab=expression(phi1[hh]),main=expression(paste("TRUE PACF of MA(1) with",theta[1]==0.2592)))
abline(h=0)

########################################## 2 #############################
# 2.c
ARMAtoMA <- ARMAtoMA(ar=0.3,ma=c(0.1,0.2),lag.max = 100)
ARMAtoMA
head(ARMAtoMA)
psiVector <- c(1,ARMAtoMA)
ARMAacf(ma=psiVector,lag.max = 20)
head(psiVector)
gamma.0 <- sum(psiVector*psiVector)
psiVector99 <- psiVector[1:99]
psiVector100 <- psiVector[2:100]
gamma.1 <- sum(psiVector99*psiVector100)

# 2. d
pho1 <- gamma.1/gamma.0
a1 <- a[1:99]
a2 <- a[2:100]
(sum(a1*a2)+0.4)/(1+sum(a*a))
head(ARMAacf(ar=0.3,ma=c(0.1,0.2),lag.max = 20))


# 2c
a.ts <- ts(a)
set <- ts.intersect(a.ts,lag(x=a.ts,k=1))
set1 <- as.matrix(set)
head(set1)
colnames(set1) <- c("original","lag1")
# must add 1*0.4 to get gamma(1)
sum(set1[,1]*set1[,2])+0.4
#2d
ARMAacf(ar=0.3,ma=c(0.1,0.2),lag.max = 20)

###################################### problem 3 #############################
# 3.a
dev.new(width=8,height=6,pointsize=10)
par(mfcol=c(1,2))
plot(y=ARMAacf(ar=0.3,ma=c(0.1,0.2),lag.max = 20),x=0:20,type = "h",ylim = c(-1,1),xlab = "h",ylab = expression(rho(h)))
abline(h=0)
plot(x=ARMAacf(ar=0.3,ma=c(0.1,0.2),lag.max = 20,pacf=TRUE),type = "h",ylim = c(-1,1),xlab = "h",ylab=expression(phi1[hh]))
abline(h=0)
TrueACF <- ARMAacf(ar=0.3,ma=c(0.1,0.2),lag.max = 20)
TrueACF[2]
c<- 0.3
0.3^3*0.2+0.3^4*0.1+0.3^5


#3.b (1)
# n=20
set.seed(1287)
n.20 <- arima.sim(model=list(order=c(1,0,2),ar=0.3,ma=c(0.1,0.2)),n=20,rand.gen = rnorm,sd=1)
head(n.20)
tail(n.20)
#n=100
set.seed(9198)
n.100 <- arima.sim(model=list(order=c(1,0,2),ar=0.3,ma=c(0.1,0.2)),n=100,rand.gen = rnorm,sd=1)
head(n.100)
tail(n.100)
#n=10000
set.seed(8712)
n.10000 <- arima.sim(model=list(order=c(1,0,2),ar=0.3,ma=c(0.1,0.2)),n=10000,rand.gen = rnorm,sd=1)
head(n.10000)
tail(n.10000)

#3.b (2) plot
# n=20
dev.new(width=8,height=6,pointsize=10)
par(mfcol=c(1,2))
acf(x=n.20,type = "correlation",lag.max = 20,ylim=c(-1,1),main="ACF for sample size 20")
pacf(x=n.20,lag.max = 20,ylim=c(-1,1),main="PACF for sample size 20")

# n=100
dev.new(width=8,height=6,pointsize=10)
par(mfcol=c(1,2))
acf(x=n.100,type = "correlation",lag.max = 20,ylim=c(-1,1),main="ACF for sample size 100")
pacf(x=n.100,lag.max = 20,ylim=c(-1,1),main="PACF for sample size 100")

# n=10000
dev.new(width=8,height=6,pointsize=10)
par(mfcol=c(1,2))
acf(x=n.10000,type = "correlation",lag.max = 20,ylim=c(-1,1),main="ACF for sample size 10000")
pacf(x=n.10000,lag.max = 20,ylim=c(-1,1),main="PACF for sample size 10000")

# 3.b (3) generate 1000 dataset for different sample size and find the mean of estimates of each sample size
# n=20
n20.summary <- matrix(data=NA,nrow = 1000,ncol=3)
set.seed(1287)
for (i in 1:1000){
  x <- arima.sim(model=list(order=c(1,0,2),ar=0.3,ma=c(0.1,0.2)),n=20,rand.gen = rnorm,sd=1)
  acfs <- acf(x=x,lag.max = 20,type = "correlation")
  acf1 <- acfs$acf[2]
  z <- acf1*sqrt(1000)
  if (z>qnorm(p=0.975) | z< -qnorm(p=0.975)){
    reject=1
  }
  else{
    reject=0
  }
  n20.summary[i,] <- c(acf1,reject,z)
}
mean(n20.summary[,1])
mean(n20.summary[,2])

# n=100
n100.summary <- matrix(data=NA,nrow = 1000,ncol=3)
set.seed(9198)
for (i in 1:1000){
  x <- arima.sim(model=list(order=c(1,0,2),ar=0.3,ma=c(0.1,0.2)),n=100,rand.gen = rnorm,sd=1)
  acfs <- acf(x=x,lag.max = 20,type = "correlation")
  acf1 <- acfs$acf[2]
  z <- acf1*sqrt(100)
  if (z>qnorm(p=0.975) | z< -qnorm(p=0.975)){
    reject=1
  }
  else{
    reject=0
  }
  n100.summary[i,] <- c(acf1,reject,z)
}
mean(n100.summary[,1])
mean(n100.summary[,2])

# n=10000 8712 
n10000.summary <- matrix(data=NA,nrow = 1000,ncol=3)
set.seed(8712)
for (i in 1:1000){
  x <- arima.sim(model=list(order=c(1,0,2),ar=0.3,ma=c(0.1,0.2)),n=10000,rand.gen = rnorm,sd=1)
  acfs <- acf(x=x,lag.max = 20,type = "correlation")
  acf1 <- acfs$acf[2]
  z <- acf1*sqrt(10000)
  if (z>qnorm(p=0.975) | z< -qnorm(p=0.975)){
    reject=1
  }
  else{
    reject=0
  }
  n10000.summary[i,] <- c(acf1,reject,z)
}
mean(n10000.summary[,1])
mean(n10000.summary[,2])

ACF.true <- ARMAacf(ar=0.3,ma=c(0.1,0.2),lag.max = 20)
ACF.true[2]






























