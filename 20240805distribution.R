rm(list = ls(all = TRUE))
n = 50
N = 1
alpha.n = (1/seq(1:n))^2
XN = matrix(NA,n,N)
for (j in 1:N) {
set.seed(j)
  y <- rnorm((n+1), mean = 0, sd = 1)
  x <- c()
  for (i in 1:n) {
    x[i] <- (alpha.n[i]*y[i+1] + y[i])/((1 + alpha.n[i]^2)^{1/2})
  }
  XN[,j] <- x # born-in first 50
}

# plot()

f1 <- function(h,t,e){
  res <- mean(pnorm((t-e)/h, mean = 0, sd = 1))
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

hn <- c(n^{-1/3}, n^{-1/4}*(log(n))^{1/4}, n^{-1/6}) #*(log(n))^{1/4}
# hn = (log(n)/n)^(1/4) #设置窗宽
# QU <- quantile(eN[,1], probs = c(0.25,0.75), na.rm = FALSE)
# IQR = QU[2]-QU[1]
# hn = IQR*n^(-1)
# hn = kCDF(eN[,1])$bw #设置窗宽 install.packages("sROC") library(sROC)
lh <- length(hn)
t = seq(-3, 3, by = 0.01)
m = length(t)
Fhatm = matrix(NA, m, lh)
for (l in 1:lh) {
Fhat = matrix(NA, m, N)
for (j in 1:N){
  Fhat.temp = c()
  for (i in 1:m){
    Fhat.temp[i] <- f1(h=hn[l], t=t[i], e = XN[,1])
  }
  Fhat[,j] = Fhat.temp
}
Fhatm[,l] <- apply(Fhat, 1, mean)
}

Fn = matrix(NA, m, N)
for (j in 1:N){
  Ft.temp = c()
  for (i in 1:m){
    Ft.temp[i] <- f2(t=t[i], e= XN[,j])
  }
  Fn[,j]= Ft.temp
}
Fnm <- apply(Fn, 1, mean)
# par(mar=c(1,1,1,1))   #设置边距参数

# hist(XN[,1],freq = F)
# curve(dnorm(x, mean = 0, sd = 1), add = T, lwd=3, col= "black", lty=1)

plot(t, Fhat[,1], xlab = "x", ylab = "distribution function", ylim=c(0,1), type = "n")
lines(t, Fhatm[,1], lwd=2, col="limegreen",lty=4)
lines(t, Fhatm[,2], lwd=2, col="blue",lty=2)
lines(t, Fhatm[,3], lwd=2, col="purple",lty=3)
lines(t, Fnm, lwd=2, col="red",lty=5)
curve(pnorm(x, mean = 0, sd = 1), add = T, lwd=3, col= "black", lty=1)
legend("bottomright", c(expression(hat(F)[n](x)(h[n]==n^{-1/3})), expression(hat(F)[n](x)(h[n]==n^{-1/4}*(log(n))^{1/4})), 
                        expression(hat(F)[n](x)(h[n]==n^{-1/6})), expression(F[n](x)), "F(x)"), 
       col = c("limegreen","blue","purple","red","black"), lty = c(4,2,3,5,1), lwd=c(2,2,2,2,3))

