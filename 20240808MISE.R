rm(list = ls(all = TRUE))
#################################################################
############# compute index MISE for Table 2 ####################
#################################################################
########## generate data  
n = 400  ###    sample size 
N = 1000     ####### replication 
# alpha.n = 1/seq(1:n)
 alpha.n = rep(0,n)
XN = matrix(NA,n,N)
for (j in 1:N) {
  y <- rnorm((n+1), mean = 0, sd = 1)
  x <- c()
  for (i in 1:n) {
    x[i] <- (alpha.n[i]*y[i+1] + y[i])/((1 + alpha.n[i]^2)^{1/2})
  }
  XN[,j] <- x # 
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
####### hazard function #########
rf <- function(t){
  out <- dnorm(t, mean = 0, sd = 1)/(1-pnorm(t, mean = 0, sd = 1))
}

hn = (log(n)/n)^(1/4) # set bandwidth 

m <- 400
t = seq(-3, 3, length = m)
######## distribution estimator ############
Fnhat = Fn = matrix(NA, m, N)
for (j in 1:N)
{ ###  kernel distribution estimator
  Fnhat.temp = c()
  for (i in 1:m){
    Fnhat.temp[i] <- f1(h=hn, t=t[i], e = XN[,j])
  }
  Fnhat[,j] = Fnhat.temp
  ### empirical distribution estimator
  Fnt.temp = c()
  for (i in 1:m){
    Fnt.temp[i] <- f2(t=t[i], e= XN[,j])
  }
  Fn[,j]= Fnt.temp
}

######## density estimator ##########
fnhat = fn = matrix(NA, m, N)
for (j in 1:N){
  fnhat.temp = c()
  for (i in 1:m){
    fnhat.temp[i] <- f3(h=hn, t=t[i], e= XN[,j])
  }
  fnhat[,j] = fnhat.temp
  fn.temp = c()
  for (i in 1:m){
    fn.temp[i] <- f4(h=hn, t=t[i], e= XN[,j])
  }
  fn[,j]= fn.temp
}

### hazard estimator
rnhat = rn = matrix(NA, m, N)
for (j in 1:N){
  rn[,j] <- fn[,j]/(1-Fn[,j])
  rnhat[,j] <- fnhat[,j]/(1 - Fnhat[,j])
}

#### compute index MISE for F(x) ######
ISE1 = NULL
for (i in 1:N){
  ISE1[i] <- sum(( Fn[,i] - pnorm(t, mean = 0, sd = 1))^2)*0.05 }
MISE.F <- mean(ISE1)

ISE2 = NULL
for (i in 1:N){
  ISE2[i] <- sum(( Fnhat[,i]- pnorm(t, mean = 0, sd = 1) )^2)*0.05 }
MISE.Fhat <- mean(ISE2)
Mratio = MISE.F/MISE.Fhat


#### compute index MISE for f(x) ######
ise1 = NULL
for (i in 1:N){
  ise1[i] <- sum(( fn[,i] - dnorm(t, mean = 0, sd = 1))^2)*0.01 }
MISE.f <- mean(ise1)

ise2 = NULL
for (i in 1:N){
  ise2[i] <- sum(( fnhat[,i]- dnorm(t, mean = 0, sd = 1))^2)*0.01 }
MISE.fhat <- mean(ise2)
Mfratio = MISE.f/MISE.fhat

###############################################
########## hazard function ########################
m <- 400
t = seq(-3, 1, length = m)
######## distribution estimator ############
Fnhat = Fn = matrix(NA, m, N)
for (j in 1:N)
{ ###  kernel distribution estimator
  Fnhat.temp = c()
  for (i in 1:m){
    Fnhat.temp[i] <- f1(h=hn, t=t[i], e = XN[,j])
  }
  Fnhat[,j] = Fnhat.temp
  ### empirical distribution estimator
  Fnt.temp = c()
  for (i in 1:m){
    Fnt.temp[i] <- f2(t=t[i], e= XN[,j])
  }
  Fn[,j]= Fnt.temp
}

######## density estimator ##########
fnhat = fn = matrix(NA, m, N)
for (j in 1:N){
  fnhat.temp = c()
  for (i in 1:m){
    fnhat.temp[i] <- f3(h=hn, t=t[i], e= XN[,j])
  }
  fnhat[,j] = fnhat.temp
  fn.temp = c()
  for (i in 1:m){
    fn.temp[i] <- f4(h=hn, t=t[i], e= XN[,j])
  }
  fn[,j]= fn.temp
}

### hazard estimator
rnhat = rn = matrix(NA, m, N)
for (j in 1:N){
  rn[,j] <- fn[,j]/(1-Fn[,j])
  rnhat[,j] <- fnhat[,j]/(1 - Fnhat[,j])
}

#### compute index MISE ######
rISE1 = NULL
for (j in 1:N){ rISE1[j] <- sum(( rn[,j] - rf(t) )^2)*0.01  }
MISE.r <- mean(rISE1)
rnhat.ISE2 = NULL
for (j in 1:N){ rnhat.ISE2[j] <- sum(( rnhat[,j]- rf(t) )^2)*0.01 }
MISE.rhat <- mean(rnhat.ISE2)
Mr.ratio = MISE.r/MISE.rhat
########################################################
################ Output the result ####################
index = c(MISE.Fhat, Mratio, MISE.fhat, Mfratio, MISE.rhat, Mr.ratio)
index <- round(index, 4)
index
