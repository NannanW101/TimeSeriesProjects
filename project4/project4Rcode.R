# read in Radon dataset 
radon <- read.csv(file="/Users/nannan/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/2022 Spring/STAT-878/project data/radon.csv",header = T)
head(radon)
x <- radon$Radon
# define diagnostic function examine.mod()
examine.mod <- function(mod.fit.obj, mod.name, max.lag = 20) {
  
  z <- mod.fit.obj$coef/sqrt(diag(mod.fit.obj$var.coef))
  p.value <- 2*(1-pnorm(q = abs(z), mean = 0, sd = 1))
  
  dev.new(width = 8, height = 6, pointsize = 10)
  tsdiag(object = mod.fit.obj, gof.lag = max.lag)
  mtext(text = paste("Model is", mod.name), side = 1, line = 3,
        at = 1, col = "red", adj = 0)
  
  dev.new(width = 8, height = 6, pointsize = 10)
  par(mfrow = c(2,2)) 
  pacf(x = as.numeric(mod.fit.obj$residuals), lag.max = max.lag,
       main = paste("Estimated PACF for residuals", mod.name),
       xlim = c(1,max.lag), ylim = c(-1,1))
  hist(x = mod.fit.obj$residuals, main = "Histogram of residuals",
       xlab = "Residuals", freq = FALSE, col = NA)
  curve(expr = dnorm(x, mean = mean(mod.fit.obj$residuals),
                     sd = sd(mod.fit.obj$residuals)), col = "red", add = TRUE)
  qqnorm(y = mod.fit.obj$residuals, ylab = "Residuals",
         panel.first = grid(col = "gray", lty = "dotted"))
  qqline(y = mod.fit.obj$residuals, col = "red")
  par(mfrow = c(1,1)) 
  
  list(z = z, p.value = p.value)
}

# step 1: plot the time series data 
plot(x=x,type = 'l',col="red")
points(x=x,pch=20,col="blue")

acf(x=x,lag.max = 20,type = "correlation")
# step 2: acf and pacf plots
par(mfrow=c(1,2))
acf(x=x,type = "correlation",lag.max = 20)
pacf(x=x,lag.max = 20)

# step 3: estimate the model
model.101 <- arima(x=x,order = c(1,0,1),include.mean = TRUE)
model.101
examine.mod(model.101,mod.name = "ARIMA(1,0,1)")
model.100 <- arima(x=x, order =c(1,0,0),include.mean = TRUE)
model.100
examine.mod(model.100,mod.name = "AR(1)")
model.001 <- arima(x=x,order = c(0,0,1),include.mean = T)
model.001
fore.mod.100 <- predict(object=model.100, n.ahead=12, se.fit = TRUE)
fore.mod.100

# step 6: calculate 95% confidence interval
low <- fore.mod.100$pred-qnorm(p=0.975,mean=0,sd=1)*fore.mod.100$se
up <- fore.mod.100$pred+qnorm(p=0.975,mean=0,sd=1)*fore.mod.100$se
data.frame(low,up)

# step 6: forecast plots
plot(x = x, ylab = expression(x[t]), xlab = "t", type = 
         "o", col = "red", lwd = 1, pch = 20, 
       panel.first = grid(col = "gray", lty = "dotted"), xlim 
       = c(1, 100))
lines(x = c(x - model.100$residuals, fore.mod.100$pred), lwd 
        = 1, col = "black", type = "o", pch = 17) 
lines(y = low, x = 89:100, lwd = 1, col = "darkgreen", 
        lty = "dashed") 
lines(y = up, x = 89:100, lwd = 1, col = "darkgreen", 
        lty = "dashed") 
legend(locator(1), legend = c("Observed", "Forecast", 
                                "95% C.I."), lty = c("solid", "solid", "dashed"),
         col = c("red", "black", "darkgreen"), pch = c(20, 
                                                       17, NA), bty = "n")
# step 6: zoom in plot
plot(x = x, ylab = expression(x[t]), xlab = "t", type = 
       "o", col = "red", lwd = 1, pch = 20, 
     panel.first = grid(col = "gray", lty = "dotted"), xlim 
     = c(80, 100))
lines(x = c(x - model.100$residuals, fore.mod.100$pred), lwd 
      = 1, col = "black", type = "o", pch = 17) 
lines(y = low, x = 89:100, lwd = 1, col = "darkgreen", 
      lty = "dashed") 
lines(y = up, x = 89:100, lwd = 1, col = "darkgreen", 
      lty = "dashed") 
legend(locator(1), legend = c("Observed", "Forecast", 
                              "95% C.I."), lty = c("solid", "solid", "dashed"),
       col = c("red", "black", "darkgreen"), pch = c(20, 
                                                     17, NA), bty = "n")

# call examine.mod function for three candidate models
examine.mod(mod.fit.obj = model.101,mod.name = "ARIMA(1,0,1)", max.lag = 20)
examine.mod(mod.fit.obj = model.100,mod.name = "ARIMA(1,0,0)", max.lag = 20)
examine.mod(mod.fit.obj = model.001,mod.name = "ARIMA(0,0,1)", max.lag = 20)

############################# problem 2
# read in TB dataset
tb <- read.csv(file = "/Users/nannan/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/2022 Spring/STAT-878/project data/tb.csv",header = TRUE)
tb.2007 <- tb[1:96,]
x.2007 <- tb.2007$count
# 2.a (1) plot the TB time series data
plot(x=x.2007,type = "l",col="red")
points(x=x.2007,pch=20,col="blue")
abline(v=c(12,24,36,48,60,72,84,96),lty="dotted")
axis(side=1,at=c(12,24,36,48,60,72,84,96),tck=-.01) 
# plot observed acf plot
acf(x=x.2007,lag.max = 90)
axis(side=1,at=c(12,24,36,48,60,72,84,96),tck=-.01) 
abline(v=c(12,24,36,48,60,72,84,96),lty="dotted")
# examine when seasonal parameter is 12
x.s12 <- diff(x=x.2007,lag=12,differences = 1)
x.d1D1 <- diff(x.s12,lag = 1,differences = 1)
plot(x.d1D1,type = "o",col="blue")
plot(x=x.s12,type = "o",col="red")
points(x=x.s12,pch=20,col="blue")

# 2.a (2) examine CDC model
model.CDC <- arima(x=x.2007,order = c(0,1,1),seasonal = list(order=c(0,1,1), period=12))
model.CDC
examine.mod(model.CDC,mod.name = "CDC model",max.lag = 50)

# 2.b model building
# step 1: Construct plots of the time series
tb <- read.csv(file = "/Users/nannan/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/2022 Spring/STAT-878/project data/tb.csv",header = TRUE)
tb.2007 <- tb[1:96,]
x.2007 <- tb.2007$count
x.s12 <- diff(x=x.2007,lag=12,differences = 1)
plot(x=x.s12,type = "o",col="red")
points(x=x.s12,col="blue",pch=20)
# step 2: Construct plots of the estimated ACF and PACF
model.000.010 <- arima(x=x.2007,order = c(0,0,0),seasonal = list(order=c(0,1,0), period=12))
model.000.010
par(mfrow=c(1,2))
acf(x=model.000.010$residuals,lag.max = 100)
abline(v=c(12,24,36,48,60,72,84,96),lty="dotted")
axis(side=1,at=c(12,24,36,48,60,72,84,96),tck=-.01) 
pacf(x=model.000.010$residuals,lag.max = 100)
abline(v=c(12,24,36,48,60,72,84,96),lty="dotted")
# step 3: Find the estimated model using MLE
model.000.111 <- arima(x=x.2007,order = c(0,0,0),seasonal = list(order=c(1,1,1), period=12))
model.000.111
# step 4: Investigate the diagnostic measures
examine.mod(model.000.111,mod.name = "model.000.111",max.lag = 80)
# propose and examine another model
model.202.111 <- arima(x=x.2007,order = c(2,0,2),seasonal = list(order=c(1,1,1), period=12))
model.202.111

# 2.c 
library(tidyverse)
# 2008 and 2009 data
tb2008 <- tb[97:108,]
tb2009 <- tb[109:120,]
mod.fit3 <- arima(x = x.2007, order = c(2, 0, 2), seasonal = list(order = c(1, 1, 1), period = 12))
mod.fit3
# forcasts based on mod.fit3
tb.2007 <- tb[1:96,]
tb.2008 <- tb[97:108,]
tb.2009 <- tb[109:120,]
x.2007 <- tb.2007$count
x.2008 <- tb.2008$count
x.2009 <- tb.2009$count
fore.mod <- predict(object = mod.fit3, n.ahead = 24, se.fit = TRUE)
fore.mod
# 95% confidence interval for 2008 observations
low.2008=fore.mod$pred[1:12]-qnorm(p=0.975)*fore.mod$se[1:12]
up.2008=fore.mod$pred[1:12]+qnorm(p=0.975)*fore.mod$se[1:12]
data.frame(x.2008,fore.mod$pred[1:12],low.2008,up.2008)
# 95% confidence interval for 2009 observations
low.2009=fore.mod$pred[13:24]-qnorm(p=0.975)*fore.mod$se[13:24]
up.2009=fore.mod$pred[13:24]+qnorm(p=0.975)*fore.mod$se[13:24]
data.frame(x.2009,fore.mod$pred[13:24],low.2009,up.2009)
# MSE for 2008 and 2009
mse.2008 <- sum((fore.mod$pred[1:12]-x.2008)^2)/12
mse.2008
mse.2009 <- sum((fore.mod$pred[13:24]-x.2008)^2)/12
mse.2009

# 2.d plot observations and forecasts
# lower bound of 95% confidence interval
low=c(low.2008,low.2009)
up=c(up.2008,up.2009)
# upper bound of 95% confidence interval
plot(x = tb$count, ylab = expression(x[t]), xlab = "t", type = "o", col = "red", lwd = 1, pch = 20, main ="Forecasted TB cases for 2008-2009",panel.first=grid(col = "gray", lty = "dotted"), xlim = c(1,120))
lines(x = c(x.2007 - mod.fit3$residuals, fore.mod$pred), lwd = 1, col = "black", type = "o", pch = 17)
lines(y = low, x =97:120, lwd = 1, col = "darkgreen", lty = "dashed") 
lines(y = up, x = 97:120, lwd = 1, col = "darkgreen", lty = "dashed")
legend(locator(1), legend = c("Observed", "Forecast","95% C.I."), lty = c("solid", "solid", "dashed"),col = c("red", "black", "darkgreen"), pch = c(20,+17, NA), bty = "n")
# zoom in forecast plot
plot(x = tb$count, ylab = expression(x[t]), xlab = "t", type = "o", col = "red", lwd = 1, pch = 20, main ="Forecasted TB cases for 2008-2009",panel.first=grid(col = "gray", lty = "dotted"), xlim = c(96,120))
lines(x = c(x.2007 - mod.fit3$residuals, fore.mod$pred), lwd = 1, col = "black", type = "o", pch = 17)
lines(y = low, x =97:120, lwd = 1, col = "darkgreen", lty = "dashed") 
lines(y = up, x = 97:120, lwd = 1, col = "darkgreen", lty = "dashed")
legend(locator(1), legend = c("Observed", "Forecast","95% C.I."), lty = c("solid", "solid", "dashed"),col = c("red", "black", "darkgreen"), pch = c(20, +17, NA), bty = "n")




