rm(list = ls(all = TRUE)) ## remove (almost) everything in the working environment.
n = 50 # sample size
N = 1
alpha.n = (1/seq(1:n))^2
XN = matrix(NA,n,N)
for (j in 1:N) {
set.seed(j+3)
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

hn <- c(n^{-1/3}, n^{-1/4}*(log(n))^{1/4}, n^{-1/6}) # bandwidth
t = seq(-3, 3, by = 0.01)
m = length(t)
lh <- length(hn)
Fhatm = matrix(NA, m, lh)
for (l in 1:lh) {
Fhat = matrix(NA, m, N)
for (j in 1:N){
  Fhat.temp = c()
  for (i in 1:m){
    Fhat.temp[i] <- f1(h=hn[l], t=t[i], e= XN[,j])
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
# par(mar=c(1,1,1,1))  


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

# par(mar=c(1,1,1,1)) 
rnm <- fnm/(1 - Fnm)
rhatm <- matrix(NA, m, lh)
for (l in 1:lh){
  rhatm[,l] <- fhatm[,l]/(1 - Fhatm[,l])
}

####  plot  ####
rf <- function(t){
#  t=0
  out <- dnorm(t, mean = 0, sd = 1)/(1 - pnorm(t, mean = 0, sd = 1))
}


plot(t, rhatm[,1], xlab = "x", ylab = "harzard function", ylim = c(0, 4), type = "n")
lines(t, rhatm[,1], lwd=2, col="limegreen",lty=4)
lines(t, rhatm[,2], lwd=2, col="blue",lty=2)
lines(t, rhatm[,3], lwd=3, col="purple",lty=3)
lines(t, rnm, lwd=2, col="red",lty=5)
curve(rf, add = T, lwd=3, col= "black", lty=1)
legend("topleft", c(expression(hat(r)[n](x)(h[n]==n^{-1/3})), expression(hat(r)[n](x)(h[n]==n^{-1/4}*(log(n))^{1/4})), 
                       expression(hat(r)[n](x)(h[n]==n^{-1/6})), expression(r[n](x)(h[n]==n^{-1/4}*(log(n))^{1/4})), "r(x)"), 
       col = c("limegreen","blue","purple","red","black"), lty = c(4,2,3,5,1), lwd=c(2,2,3,2,3))
