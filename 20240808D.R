rm(list = ls(all = TRUE))  ## remove (almost) everything in the working environment.
#################################################################
############# compute index D(F) for Table 1 ####################
#################################################################
########## generate data  
n = 100  ###    sample size 
N = 1000     ####### replication 
 alpha.n = (1/seq(1:n))^2
# alpha.n = rep(0,n)
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

#### compute index D(F) ######
DF = DFhat = NULL
DF.temp =  DFhat.temp = matrix(NA,N,1)
for (i in 1:N) {
  DF.temp[i] <- max(Fn[,i]- pnorm(t, mean = 0, sd = 1))
  DFhat.temp[i] <- max(Fnhat[,i] - pnorm(t, mean = 0, sd= 1))
}
DF <- mean(DF.temp)
DFhat <- mean(DFhat.temp)
Dratio = DF/DFhat


#### compute index D(f) ######
Df = Dfhat = NULL
Df.temp=  Dfhat.temp= matrix(NA,N,1)
for (i in 1:N) {
  Df.temp[i] <- max(fn[,i] - dnorm(t, mean = 0, sd = 1))
  Dfhat.temp[i] <- max(fnhat[,i] - dnorm(t, mean = 0, sd = 1))
}
Df <- mean(Df.temp)
Dfhat <- mean(Dfhat.temp)
Dfratio = Df/Dfhat

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

#### compute index D(r) ######
rf <- function(t){
  out <- dnorm(t, mean = 0, sd = 1)/(1-pnorm(t, mean = 0, sd = 1))
}
Dr = Drhat = NULL
Dr.temp =  Drhat.temp = matrix(NA, N, 1)
for (j in 1:N) {
  Dr.temp[j] <- max(rn[,j] - rf(t) )
  Drhat.temp[j] <- max(rnhat[,j] - rf(t) )
}
Dr <- mean(Dr.temp)
Drhat <- mean(Drhat.temp)
Dr.ratio = Dr/Drhat

index = c(DFhat,Dratio, Dfhat, Dfratio, Drhat, Dr.ratio)
index <- round(index, 4)
index
