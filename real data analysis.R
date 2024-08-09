rm(list = ls(all = TRUE))
library(survival)
data(cancer, package = "survival")

gbsg1 <- gbsg$rfstime[which(gbsg$age < 46)]

gbsg2 <- gbsg$rfstime[which(46 <= gbsg$age & gbsg$age <=55 )]

gbsg3 <- gbsg$rfstime[which(56 <= gbsg$age & gbsg$age <=65 )]

gbsg4 <- gbsg$rfstime[which(gbsg$age > 65)]
###################### acf ################
par(mfrow=c(2,2))
acf(gbsg1, lag.max = NULL, type = "correlation", plot = TRUE, na.action = na.fail, demean = TRUE, main ="")
title("(a)  Age <= 46")
acf(gbsg2, lag.max = NULL, type = "correlation", plot = TRUE, na.action = na.fail, demean = TRUE, main ="")
title("(b)  46 <= Age <= 55")
acf(gbsg3, lag.max = NULL, type = "correlation", plot = TRUE, na.action = na.fail, demean = TRUE, main ="")
title("(c) 56<= Age <= 65")
acf(gbsg4, lag.max = NULL, type = "correlation", plot = TRUE, na.action = na.fail, demean = TRUE, main ="")
title("(d) Age > 65")


F1 <- function(h,t,e){
  res <- mean(pnorm((t-e)/h, mean = 0, sd = 1))
  return(res)
}

F2 <- function(t,e)
{
  res <- sum(as.numeric(e <= t))/length(e)
  return(res)
}

F3 <- function(t,e){
  sigmax <- sd(e)
  n <- length(e)
  i <- seq(1,n,by=1)
  h <- 2*sigmax*i^(-1/4)
  res <- mean(pnorm((t-e)/h, mean = 0, sd = 1))
}

f1 <- function(t,e){
  n <- length(e)
  i <- seq(1,n,by=1)
  hi <- (log(i)/i)^(1/6)
  res <- mean(qnorm((t-e)/hi,mean=0,sd=1)/hi)
  return(res)
}

f2 <- function(t,e){
  res <- sum(as.numeric(e <= t))/length(e)
}

f3 <- function(h,t,e){
  res = 1/(sqrt(2*pi)*h)*mean(exp(-(t-e)^2/(2*h^2))) 
  return(res)
}

f4 <- function(h, t, e){
  res = length(which( (t-h) < e & e<= (t+h)))/(2*length(e)*h)
}


gbsg <- gbsg1
n <- length(gbsg)
hn = sd(gbsg)*(log(n)/n)^(1/4) # bandwidths
# hn = n1^(-1/5)
m <- 400
t = seq(0, 2500, length = m)

fhat <- fn <- matrix(NA, m, 1)
for (i in 1:m){
  fhat[i] <- f3(h=hn, t=t[i], e = gbsg)
}

for (i in 1:m){
  fn[i] <- f4(h=hn, t=t[i], e= gbsg)
}

Fn <- Fhat <- matrix(NA, m, 1)
for (i in 1:m){
  Fn[i] <- F2(t=t[i], e = gbsg)
}
for (i in 1:m){
  Fhat[i] <- F1(h=hn, t=t[i], e = gbsg)
}

### 1 - F(x) distribution function ###
par(mai=c(0.7,.7,.4,.4),cex=0.8)
plot(t,(1-Fhat), xlab="x", ylab= "survivor functions", type = "n", ylim = c(0, 1.05),lwd=2)
lines(t, 1-Fn, col="blue", lty=3, lwd=2)
lines(t, (1-Fhat), col="red",lty=4,lwd=2)
title("(d)  Age > 65")
legend("topright", c( expression(1-F[n](x)), expression(1-hat(F)[n](x))), 
       col = c("blue","red"), lty = c(3,4), lwd=2)


### hazard function ########
rhat = rn = matrix(NA, m, 1)
rn <- fn/(1 - Fn)
rhat <- fhat/(1 - Fhat)
par(mai=c(0.7,0.7,0.5,0.5),cex=0.8)
plot(t, rn, type = "n", xlab = "x", ylab = "harzard functions")
# plot(rn)
lines(t,rn,col="orange",lty=2,lwd=2)
lines(t,rhat,col="green",lty=5,lwd=2)
title("(a)  Age <= 46")
legend("topleft", c(expression(r[n](x)), expression(hat(r)[n](x))), 
       col = c("orange","green"), lty = c(2,5), lwd=2)


