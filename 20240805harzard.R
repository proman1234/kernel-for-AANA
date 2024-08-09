rm(list = ls(all = TRUE))
n = 600 # sample size
N = 1
alpha.n = 1/seq(1:n)
XN = matrix(NA,n,N)
for (j in 1:N) {
  y <- rnorm((n+1), mean = 0, sd = 1)
  x <- c()
  for (i in 1:n) {
    x[i] <- (alpha.n[i]*y[i+1] + y[i])/((1 + alpha.n[i]^2)^{1/2})
  }
  XN[,j] <- x 
}


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

hn <- c(n^{-1/2}, n^{-1/3}, n^{-1/5}) # bandwidth
t = seq(-3, 3, by = 0.01)
m = length(t)
Fhat = matrix(NA, m, 3)
for (l in 1:3){
  Fhat.temp = c()
  for (i in 1:m){
    Fhat.temp[i] <- f1(h=hn[l], t=t[i], e= XN[,j])
  }
  Fhat[,l] = Fhat.temp
}

Fn = matrix(NA, m, 1)
for (j in 1:N){
  Ft.temp = c()
  for (i in 1:m){
    Ft.temp[i] <- f2(t=t[i], e= XN[,j])
  }
  Fn[,j]= Ft.temp
}
# par(mar=c(1,1,1,1))   # 设置边距参数

fhat = matrix(NA, m, 3)
for (l in 1:3){
  fhat.temp = c()
  for (i in 1:m){
    fhat.temp[i] <- f3(h=hn[l], t=t[i], e= XN[,j])
  }
  fhat[,l] = fhat.temp
}
fn <- matrix(NA, m, 1)
for (j in 1:N){
  ft.temp = c()
  for (i in 1:m){
    ft.temp[i] <- f4(h=hn[2], t=t[i], e= XN[,j])
  }
  fn[,j]= ft.temp
}
# par(mar=c(1,1,1,1))  # 设置边距参数


rn <- fn[,1]/(1 - Fn[,1])
rhat <- matrix(NA, m, 3)
for (l in 1:3){
  rhat[,l] <- fhat[,l]/(1 - Fhat[,l])
}


####  plot  ####
rf <- function(t){
#  t=0
  out <- dnorm(t, mean = 0, sd = 1)/(1 - pnorm(t, mean = 0, sd = 1))
}
length(rn)
length(t)
length(fn)
length(Fn[,1])

plot(t, rhat[,1], xlab = "x", ylab = "harzard function", ylim = c(0, 4), type = "n")
lines(t, rhat[,1], lwd=2, col="limegreen",lty=4)
lines(t, rhat[,2], lwd=2, col="blue",lty=2)
lines(t, rhat[,3], lwd=3, col="purple",lty=3)
lines(t, rn, lwd=2, col="red",lty=5)
curve(rf, add = T, lwd=3, col= "black", lty=1)
legend("topleft", c(expression(h[n]==n^{-1/2}), expression(h[n]==n^{-1/3}), 
                     expression(h[n]==n^{-1/5}), expression(r[n](x)), "r(x)"), 
       col = c("limegreen","blue","purple","red","black"), lty = c(4,2,3,5,1), lwd=c(2,2,3,2,3))
