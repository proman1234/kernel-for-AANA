rm(list = ls(all = TRUE))
n = 200 # sample size
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


t = seq(-3, 3, by = 0.01)
m = length(t)
hn <- c(n^{-1/3}, n^{-1/4}*(log(n))^{1/4}, n^{-1/6}) #*(log(n))^{1/4}
lh <- length(hn)
fhatm = matrix(NA, m, lh)
for (l in 1:lh) {
  fhat = matrix(NA, m, N)
  for (j in 1:N){
    fhat.temp = c()
    for (i in 1:m){
      fhat.temp[i] <- f3(h=hn[l], t=t[i], e = XN[,j])
    }
    fhat[,j] = fhat.temp
  }
  fhatm[,l] <- apply(fhat, 1, mean)
}


fn = matrix(NA, m, N)
for (j in 1:N){
  ft.temp = c()
  for (i in 1:m){
    ft.temp[i] <- f4(h=hn[2], t=t[i], e= XN[,j])
  }
  fn[,j]= ft.temp
}
fnm <- apply(fn, 1, mean)

plot(t, fhatm[,1], xlab = "x", ylab = "density function", ylim=c(0,0.55), type = "n")
lines(t, fhatm[,1], lwd=2, col="limegreen",lty=4)
lines(t, fhatm[,2], lwd=2, col="blue",lty=2)
lines(t, fhatm[,3], lwd=2, col="purple",lty=3)
lines(t, fnm, lwd=2, col="red",lty=5)
curve(dnorm(x, mean = 0, sd = 1), add = T, lwd=3, col= "black", lty=1)
legend("topright", c(expression(hat(f)[n](x)(h[n]==n^{-1/3})), expression(hat(f)[n](x)(h[n]==n^{-1/4}*(log(n))^{1/4})), 
                     expression(hat(f)[n](x)(h[n]==n^{-1/6})), expression(f[n](x)(h[n]==n^{-1/4}*(log(n))^{1/4})), "f(x)"), 
       col = c("limegreen","blue","purple","red","black"), lty = c(4,2,3,5,1), lwd=c(2,2,2,2,3))

