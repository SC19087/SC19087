## -----------------------------------------------------------------------------
peBS <- function(x, y, h){ # probability estimator of BS density

  m <- length(x)
  n <- length(y)
  ye <- rep(0, m)
  for (i in 1:n){
    ye <- ye + as.numeric((x >= y[i] - h) & (x < y[i] + h))
  }
  ye <- ye / (2*h*n)
  return(ye)
}

## -----------------------------------------------------------------------------
rBS <- function(n){ 
  u <- runif(n)
  y <- u
  ind <- which(u > 0.5) #index for those generated from N(0,1)
  y[ind] <- rnorm(length(ind), 0, 1)
  for (j in 0:4) {
    ind <- which(u > j * 0.1 & u <= (j+1) * 0.1)
    #index for those generated from N(j/2-1,1/10^2)
    y[ind] <- rnorm(length(ind), j/2 -1, 1/10)
  }
  return(y)
}

## -----------------------------------------------------------------------------
data<-rnorm(100,mean=0,sd=1)
shapiro.test(data)

## -----------------------------------------------------------------------------
library(vcd)
knitr::kable(head(Arthritis))

## -----------------------------------------------------------------------------
qqnorm(data)

## -----------------------------------------------------------------------------

Rayleigh<-function(sigma)
{ U<-runif(10000)
  X<-1:10000
  for(i in 1:10000){
  X[i]<-sqrt((-2*(sigma)^2)*log(1-U[i]))}
  X
}  
checkfunc<-function(sigma){
  x<-Rayleigh(sigma)
  t=1/(sigma*sqrt(exp(1)))
  hist(x,prob=TRUE,main="check the pdf",ylim=c(0,t))
  y<-seq(0,5,.01)
  lines(y,y*exp(-y^2/(2*sigma^2))/(sigma)^2)
}

## -----------------------------------------------------------------------------
head(Rayleigh(1))

## -----------------------------------------------------------------------------
checkfunc(1)

## -----------------------------------------------------------------------------
head(Rayleigh(0.5))

## -----------------------------------------------------------------------------
checkfunc(0.5)  

## -----------------------------------------------------------------------------
head(Rayleigh(0.75))

## -----------------------------------------------------------------------------
checkfunc(0.75)  

## -----------------------------------------------------------------------------
head(Rayleigh(0.25))

## -----------------------------------------------------------------------------
checkfunc(0.25)  

## -----------------------------------------------------------------------------
normalmix<-function(p){
  U<-runif(1000)
  X<-1:1000
  for(i in 1:1000){
    X[i]=(U[i]<=p)*rnorm(n=1,mean=0,sd=1)+(U[i]>p)*rnorm(n=1,mean=3,sd=1)
  }
  X
}
refunc<-function(p){
  x<-normalmix(p)
  cat("The head of sample for p1=",p,"is:",head(x),"\n")
  hist(x,prob=TRUE,ylim=c(0,0.8))
  y<-seq(-10,10,0.01)
  lines(y,p*dnorm(x=y,mean=0,sd=1)+(1-p)*dnorm(x=y,mean=3,sd=1))
}

## -----------------------------------------------------------------------------
refunc(0.75)

## -----------------------------------------------------------------------------
refunc(0.25)

## -----------------------------------------------------------------------------
refunc(0.5)

## -----------------------------------------------------------------------------
refunc(0.2)

## -----------------------------------------------------------------------------
refunc(0.3)

## -----------------------------------------------------------------------------
refunc(0.7)

## -----------------------------------------------------------------------------
refunc(0.8)

## -----------------------------------------------------------------------------

wishart<-function(n) #n:sample size
{
  arr<-array(1:(9*n),dim=c(3,3,n))  #save the T
  matrix<-array(1:(9*n),dim=c(3,3,n))  #save the LAL^T
  t<-rnorm(3*n) #tij:i>j
  s1<-rchisq(n,df=5)
  s2<-rchisq(n,df=4)
  s3<-rchisq(n,df=3)
  t1<-sqrt(s1) #t11
  t2<-sqrt(s2) #t22
  t3<-sqrt(s3) #t33
  S<-matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),ncol=3)
  l<-chol(S)  #matrix l is a upper triangular,t(l) is what we need
  for(i in 1:n){
    arr[,,i]<-matrix(c(t1[i],t[3*i-2],t[3*i-1],0,t2[i],t[3*i],0,0,t3[i]),ncol=3)
    matrix[,,i]<-t(l)%*%arr[,,i]%*%t(arr[,,i])%*%l
  }
  matrix
}


## -----------------------------------------------------------------------------

wishart(10)


## -----------------------------------------------------------------------------

m<-1e4
T<-runif(m,min=0,max=pi/3)  
theta.hat<-mean(sin(T)*pi/3)  
print(c(theta.hat,-cos(pi/3)+cos(0)))  

## -----------------------------------------------------------------------------
m<-1e4
MC1<-1:1000
MC2<-1:1000
set.seed(123)
for(i in 1:1000){

u<-runif(m/2)
v<-runif(m/2)
x1<-c(u,1-u)       ##sample with antithetic reduction
x2<-c(u,v) ##sample without variance reduction
MC1[i]<-mean(exp(-x1)/(1+x1^2))
MC2[i]<-mean(exp(-x2)/(1+x2^2))
}
print(c(mean(MC1),mean(MC2),1-var(MC1)/var(MC2)))

## -----------------------------------------------------------------------------
a<-c(0.2,0.4,0.6,0.8)
p<--log(1-(1-exp(-1))*a)
u1<-runif(2000)
u2<-runif(2000)
u3<-runif(2000)
u4<-runif(2000)
u5<-runif(2000)
x1<--log(1-(1-exp(-1))*u1/5)     ## xi: the random sample in ith subinterval
x2<--log(exp(-p[1])-(1-exp(-1))*u2/5)
x3<--log(exp(-p[2])-(1-exp(-1))*u3/5)
x4<--log(exp(-p[3])-(1-exp(-1))*u4/5)
x5<--log(exp(-p[4])-(1-exp(-1))*u5/5)
y1<-exp(-x1)/(1+x1^2)/(5*exp(-x1)/(1-exp(-1)))
y2<-exp(-x2)/(1+x2^2)/(5*exp(-x2)/(1-exp(-1)))
y3<-exp(-x3)/(1+x3^2)/(5*exp(-x3)/(1-exp(-1)))
y4<-exp(-x4)/(1+x4^2)/(5*exp(-x4)/(1-exp(-1)))
y5<-exp(-x5)/(1+x5^2)/(5*exp(-x5)/(1-exp(-1)))
fg<-y1+y2+y3+y4+y5
mean(fg)
sd(fg)

## -----------------------------------------------------------------------------
n<-20
alpha<-0.05
m=1e3
UCL1<-numeric(m)
LCL1<-numeric(m)
judge1<-numeric(m)
UCL2<-numeric(m)
LCL2<-numeric(m)
judge2<-numeric(m)
for(i in 1:m){
x<-rchisq(n,df=2)
UCL1[i]<-mean(x)+qt(1-alpha/2,n-2)*sd(x)/sqrt(n)
LCL1[i]<-mean(x)-qt(1-alpha/2,n-2)*sd(x)/sqrt(n)
UCL2[i]<-(n-1)*var(x)/qchisq(alpha,df=n-1) 
judge1[i]<-(UCL1[i]>2&&LCL1[i]<2)+0  ##t-interval
judge2[i]<-(UCL2[i]>4)+0  ##interval for variance
}  ## for chisq distribution with df=n,the mean of r.v. is n,the var is 2n.X~chisq(2),so mean(X)=2,var(X)=4


## -----------------------------------------------------------------------------
mean(judge1)

## -----------------------------------------------------------------------------
mean(judge2)

## -----------------------------------------------------------------------------
n<-1e3
m<-1e4
sk<-numeric(n)
for(i in 1:n){
  x<-rnorm(m)  ## generate samples from N(0,1)
  xbar<-mean(x)
  m3<-mean((x-xbar)^3)
  m2<-mean((x-xbar)^2)
  sk[i]<-m3/m2^1.5
}
sk<-sort(sk)
xq<-sk[c(floor(0.025*n),floor(0.05*n),ceiling(0.95*n),ceiling(0.975*n))]
print(xq)##the estimated quantiles  

## -----------------------------------------------------------------------------
q<-c(0.025,0.05,0.95,0.975)
var_xq<-q*(1-q)/(n*dnorm(xq,mean=0,sd=sqrt(6*(n-2)/(n+1)/(n+3)))^2)
sd_xq<-sqrt(var_xq)
print(sd_xq)

## -----------------------------------------------------------------------------
print(xq)  
print(qnorm(c(0.025,0.05,0.95,0.975),mean=0,sd=sqrt(6/n)))

## -----------------------------------------------------------------------------

n<-1e4
m<-1e4
sk<-numeric(n)
for(i in 1:n){
  x<-rnorm(m)  ## generate samples from N(0,1)
  xbar<-mean(x)
  m3<-mean((x-xbar)^3)
  m2<-mean((x-xbar)^2)
  sk[i]<-m3/m2^1.5
}
sk<-sort(sk)
xq<-sk[c(floor(0.025*n),floor(0.05*n),ceiling(0.95*n),ceiling(0.975*n))]
print(xq)
print(qnorm(c(0.025,0.05,0.95,0.975),mean=0,sd=sqrt(6/n)))

## -----------------------------------------------------------------------------
q<-c(0.025,0.05,0.95,0.975)
var_xq<-q*(1-q)/(n*dnorm(xq,mean=0,sd=sqrt(6*(n-2)/(n+1)/(n+3)))^2)
sd_xq<-sqrt(var_xq)
print(sd_xq)

## -----------------------------------------------------------------------------
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
alpha<-0.1
n<-30
m<-2500
cv<-qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
sktests1<-sktests2<-numeric(m)
for(i in 1:m){
  x1<-rbeta(n,alpha,alpha)
  x2<-rt(n,df=2)
  sktests1[i]<-as.integer(abs(sk(x1)>=cv))
  sktests2[i]<-as.integer(abs(sk(x2)>=cv))
}
pwj1<-mean(sktests1)
pwj2<-mean(sktests2)
print(c(pwj1,pwj2))


## -----------------------------------------------------------------------------
n<-20
alpha=0.05
mu0<-1#these three distribution have the same mean

m<-10000
p1<-p2<-p3<-numeric(m)
for(j in 1:m){
  x1<-rchisq(n,df=1)##X~chisq(1)
  ttest1<-t.test(x1,alternative="two.sided",mu=mu0)
  p1[j]<-ttest1$p.value
  x2<-runif(n,min=0,max=2)##X~U(0,2)
  ttest2<-t.test(x2,alternative="two.sided",mu=mu0)
  p2[j]<-ttest2$p.value
  x3<-rexp(n,rate=1)##X~exp(1)
  ttest3<-t.test(x3,alternative="two.sided",mu=mu0)
  p3[j]<-ttest3$p.value
}
p1.hat<-mean(p1<alpha)
p2.hat<-mean(p2<alpha)
p3.hat<-mean(p3<alpha)
print(c(p1.hat,p2.hat,p3.hat))

## -----------------------------------------------------------------------------
library(bootstrap)
di<-function(x,y){
  plot(scor[,x],scor[,y],xlab=colnames(scor)[x],ylab=colnames(scor)[y])
}

for(i in 1:4){
  for(j in (i+1):5){
    di(i,j)
  }
}

## -----------------------------------------------------------------------------
cor(scor)  

## -----------------------------------------------------------------------------
B<-200
n<-nrow(scor)
rho12<-rho34<-rho35<-rho45<-numeric(B)
for(j in 1:B){
  i<-sample(1:n,size=n,replace=TRUE)# a sample
  mec<-scor$mec[i]          
  vec<-scor$vec[i]
  alg<-scor$alg[i]
  ana<-scor$ana[i]
  sta<-scor$sta[i]
  rho12[j]<-cor(mec,vec)   # the correlation coefficient of the sample
  rho34[j]<-cor(alg,ana)
  rho35[j]<-cor(alg,sta)
  rho45[j]<-cor(ana,sta)
}
print(c(sd(rho12),sd(rho34),sd(rho35),sd(rho45)))

## -----------------------------------------------------------------------------
n<-20
B<-200
N<-1e3
alpha<-0.05
sknormal<-skchisq<-numeric(B)# save sample skewness
judge11<-judge12<-judge13<-numeric(N)#three judge variables for normal population
judge21<-judge22<-judge23<-numeric(N)#three judge variables for chisq population
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
for(i in 1:N){
  x<-rnorm(n)
  y<-rchisq(n,df=5)
  for(j in 1:B){
    b<-sample(1:n,size=n,replace=TRUE)
    samp1<-x[b]
    samp2<-y[b]
    sknormal[j]<-sk(samp1)
    skchisq[j]<-sk(samp2)
  }
  sknormal<-sort(sknormal)
  hat1<-mean(sknormal)
  se1<-sd(sknormal)
  skchisq<-sort(skchisq)
  hat2<-mean(skchisq)
  se2<-sd(skchisq)
  judge11[i]<-(0>hat1-qnorm(1-alpha/2)*se1)+(0>hat1+qnorm(1-alpha/2)*se1)+0
  judge12[i]<-(0>hat1*2-sknormal[B*(1-alpha/2)])+(0>hat1*2-sknormal[B*alpha/2])+0
  judge13[i]<-(0>sknormal[B*alpha/2])+(0>sknormal[B*(1-alpha/2)])
  judge21[i]<-((4/sqrt(10))>hat2-qnorm(1-alpha/2)*se2)+((4/sqrt(10))>hat2+qnorm(1-alpha/2)*se2)+0
  judge22[i]<-((4/sqrt(10))>hat2*2-skchisq[B*(1-alpha/2)])+((4/sqrt(10))>hat2*2-skchisq[B*alpha/2])+0
  judge23[i]<-((4/sqrt(10))>skchisq[B*alpha/2])+((4/sqrt(10))>skchisq[B*(1-alpha/2)]) 
}
##For normal population
##standard normal bootstrap confidence interval
cat("The coverage probability,the proportion of miss on the left,the proportion of miss on the right is",c(sum(judge11==1),sum(judge11==0),sum(judge11==2))/N)
##basic bootstrap confidence interval
cat("The coverage probability,the proportion of miss on the left,the proportion of miss on the right is",c(sum(judge12==1),sum(judge12==0),sum(judge12==2))/N)
##percentile confidence interval
cat("The coverage probability,the proportion of miss on the left,the proportion of miss on the right is",c(sum(judge13==1),sum(judge13==0),sum(judge13==2))/N)
##For chisq population
##standard normal bootstrap confidence interval
cat("The coverage probability,the proportion of miss on the left,the proportion of miss on the right is",c(sum(judge21==1),sum(judge21==0),sum(judge21==2))/N)
##basic bootstrap confidence interval
cat("The coverage probability,the proportion of miss on the left,the proportion of miss on the right is",c(sum(judge22==1),sum(judge22==0),sum(judge22==2))/N)
##percentile confidence interval
cat("The coverage probability,the proportion of miss on the left,the proportion of miss on the right is",c(sum(judge23==1),sum(judge23==0),sum(judge23==2))/N)

## -----------------------------------------------------------------------------
library(bootstrap)
n<-nrow(scor)
lambda_hat <- eigen(cov(scor))$values
theta_hat <- lambda_hat[1] / sum(lambda_hat)
theta_j <- rep(0,n)
for (i in 1:n) {
x<-scor [-i,]
lambda<-eigen(cov(x))$values
theta_j[i]<-lambda[1]/sum(lambda)
}
bias_jack<-(n-1)*(mean(theta_j)-theta_hat)
# the estimated bias of theta_hat
se_jack<-(n-1)*sqrt(var(theta_j)/n)
# the estimated se of theta_hat
print(c(bias_jack,se_jack))

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n<-length(magnetic)
e1<-e2<-e3<-e4<-numeric(n)
for(k in 1:n){
  y<-magnetic[-k]
  x<-chemical[-k]
  
  L1<-lm(y~x)
  yhat1<-L1$coef[1]+L1$coef[2]*chemical[k]
  e1[k]<-magnetic[k]-yhat1
  
  L2<-lm(y~x+I(x^2))
  yhat2<-L2$coef[1]+L2$coef[2]*chemical[k]+L2$coef[3]*chemical[k]^2
  e2[k]<-magnetic[k]-yhat2
  
  L3<-lm(log(y)~x)
  logyhat3<-L3$coef[1]+L3$coef[2]*chemical[k]
  yhat3<-exp(logyhat3)
  e3[k]<-magnetic[k]-yhat3
  
  L4<-lm(y~x+I(x^2)+I(x^3))
  yhat4<-L4$coef[1]+L4$coef[2]*chemical[k]+L4$coef[3]*chemical[k]^2+L4$coef[4]*chemical[k]^3
  e4[k]<-magnetic[k]-yhat4
  
}
c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))

## -----------------------------------------------------------------------------
LM2<-lm(magnetic~chemical+I(chemical^2))
cat("Y=",LM2$coef[1],LM2$coef[2],"X","+",LM2$coef[3],"X^2")

## -----------------------------------------------------------------------------

LM1<-lm(magnetic~chemical)
LM2<-lm(magnetic~chemical+I(chemical^2))
LM3<-lm(log(magnetic)~chemical)
LM4<-lm(magnetic~chemical+I(chemical^2)+I(chemical^3))
c(summary(LM1)$adj.r.squared,summary(LM2)$adj.r.squared,summary(LM3)$adj.r.squared,summary(LM4)$adj.r.squared)

## -----------------------------------------------------------------------------

cat("Y=",LM2$coef[1],LM2$coef[2],"X","+",LM2$coef[3],"X^2")

## -----------------------------------------------------------------------------
R<-999
n1<-20
n2<-30
mu1<-mu2<-0
sigma1<-sigma2<-1
x<-rnorm(n1,mean=mu1,sd=sigma1)
y<-rnorm(n2,mean=mu2,sd=sigma2)

z<-c(x,y)
K<-1:(n1+n2)
D<-numeric(R)
options(warn=-1)
D0<-ks.test(x,y,exact=FALSE)$statistic
for(i in 1:R){
  k<-sample(K,size=n1,replace=FALSE)
  x1<-z[k]
  y1<-z[-k]
  D[i]<-ks.test(x1,y1,exact=FALSE)$statistic
}
p<-mean(c(D0,D)>=D0)
options(warn=0)
print(p)

## -----------------------------------------------------------------------------
hist(D, main = "", freq = FALSE, breaks = "scott")
points(D0, 0, cex = 1, pch = 16)  


## -----------------------------------------------------------------------------
library(Ball)
library(boot)
library(mvtnorm)
alpha<-0.1
n0<-c(10,20,30,40,50,60,70,80,90,100)
dCov <- function(x, y) {
x <- as.matrix(x)
y <- as.matrix(y)
n <- nrow(x)
m <- nrow(y)
if (n != m || n < 2) stop("Sample sizes must agree")
if (! (all(is.finite(c(x, y)))))
stop("Data contains missing or infinite values")
Akl <- function(x) {
d <- as.matrix(dist(x))
m <- rowMeans(d)
M <- mean(d)
a <- sweep(d, 1, m)
b <- sweep(a, 2, m)
return(b + M)
}
A <- Akl(x)
B <- Akl(y)
dCov <- sqrt(mean(A * B))
dCov
}
pd1<-pd2<-pb1<-pb2<-numeric(100)
pod1<-pod2<-pob1<-pob2<-numeric(10)
ndCov2 <- function(z, ix, dims) {
  #dims contains dimensions of x and y 
  p <- dims[1]
  q <- dims[2]
  d <- p + q
  x <- z[ , 1:p] #leave x as is
  y <- z[ix, -(1:p)] #permute rows of y 
  return(nrow(z) * dCov(x, y)^2)
}
for(i in 1:10){
  for(j in 1:100){
  I2<-matrix(c(1,0,0,1),nrow=2)
  x<-rmvnorm(n0[i],rep(0,2),I2)
  e<-rmvnorm(n0[i],rep(0,2),I2)
  y1<-x/4+e
 
  y2<-x/4*e
  z1<-cbind(x,y1)
  z2<-cbind(x,y2)
  boot.obj1<-boot(data=z1,statistic=ndCov2,R=99,sim="permutation",dims=c(2,2))
  boot.obj2<-boot(data=z2,statistic=ndCov2,R=99,sim="permutation",dims=c(2,2))
  tb1<-c(boot.obj1$t0,boot.obj1$t)
  tb2<-c(boot.obj2$t0,boot.obj2$t)
  pd1[j]<-mean(tb1>=tb1[1])
  pd2[j]<-mean(tb2>=tb2[1])
  pb1[j]<-bcov.test(z1[,1:2],z1[,3:4],R=100,seed=j*12345)$p.value
  pb2[j]<-bcov.test(z2[,1:2],z2[,3:4],R=100,seed=j*23456)$p.value
  }
  pod1[i]<-mean(pd1<alpha)
  pod2[i]<-mean(pd2<alpha)
  pob1[i]<-mean(pb1<alpha)
  pob2[i]<-mean(pb2<alpha)
}
## Model 1
plot(n0,pod1,type="l",xlab='sample size',ylab='power',ylim=c(0,1),main='Model 1')
lines(n0,pob1,type="l")
points(n0,pod1)
points(n0,pob1)
## Model 2
plot(n0,pod2,type="l",xlab='sample size',ylab='power',ylim=c(0,1),main='Model 2')
lines(n0,pob2,type="l")
points(n0,pod2)
points(n0,pob2)

## -----------------------------------------------------------------------------
sigma<-c(0.05,0.5,2,16)
dl<-function(x){  ## pdf of standard Laplace
  return(exp(-abs(x))/2)
}
rw.M<-function(sigma,N){   ## generate sample for different sigma
  x0<-rnorm(1,mean=0,sd=sigma)
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for(i in 2:N){
    y<-rnorm(1,mean=x[i-1],sd=sigma)
    if(u[i]<=dl(y)/dl(x[i-1])) {
      x[i]<-y
      k<-k+1}
    else x[i]<-x[i-1]
  }
  return(list(x=x,k=k))
}

N<-2000
rw1<-rw.M(sigma[1],N)  ## four samples
rw2<-rw.M(sigma[2],N)
rw3<-rw.M(sigma[3],N)
rw4<-rw.M(sigma[4],N)

plot(1:2000,rw1$x,main=expression(sigma==0.05),type="l",xlab=' ',ylab='X')
plot(1:2000,rw2$x,main=expression(sigma==0.5),type="l",xlab=' ',ylab='X')
plot(1:2000,rw3$x,main=expression(sigma==2),type="l",xlab=' ',ylab='X')
plot(1:2000,rw4$x,main=expression(sigma==16),type="l",xlab=' ',ylab='X')
cat('the acceptance rate is',c(rw1$k,rw2$k,rw3$k,rw4$k)/N)

## -----------------------------------------------------------------------------
x<-1:100
y<-log(exp(x))
z<-exp(log(x))
print(c(sum((x-y>0)),sum(y-z>0),sum(x-z>0)))
cat(all.equal(x,y),all.equal(y,z),all.equal(x,z))

## -----------------------------------------------------------------------------
A114<-A115<-numeric(25)
k<-c(4:25,100,500,1000)
S<-function(k,x){
  return(1-pt(sqrt(x^2*k/(k+1-x^2)),df=k))
         }
f114<-function(k,x){
  return(S(k,x)-S(k-1,x))
}
for(i in 1:25){
  f<-function(x){
    return(f114(k[i],x))
  }
  b<-sqrt(k[i]-1)
  A114[i]<-uniroot(f,c(0.00001,b-0.00001))$root
}
ck<-function(k,a){
  sqrt(a^2*k/(k+1-a^2))
  }
S2<-function(k,a){

  
  log(2)+lgamma((k+1)/2)-0.5*log(pi*k)-lgamma(k/2)+log(integrate(function(u){(1+u^2/k)^(-(k+1)/2)},0,ck(k,a))$value)-(log(2)+lgamma((k)/2)-0.5*log(pi*(k-1))-lgamma((k-1)/2)+log(integrate(function(u){(1+u^2/(k-1))^(-k/2)},0,ck(k-1,a))$value))

}


for(i in 1:22){
  b<-A114[i]
  A115[i]<-uniroot(S2,interval=c(b-0.1,b+0.1),k=k[i])$root
}

## -----------------------------------------------------------------------------
print(A114[1:22])
print(A115[1:22])
print(A114[1:22]-A115[1:22])

## -----------------------------------------------------------------------------
iter<-function(p,q){
  r<-1-p-q
  a<-2*28*(p^2)/(p^2+2*p*r)+2*p*r/(p^2+2*p*r)*28+70
  b<-2*24*(q^2)/(q^2+2*q*r)+2*q*r/(q^2+2*q*r)*24+70
  c<-2*41+2*p*r/(p^2+2*p*r)*28+2*q*r/(q^2+2*q*r)*24
  pnew<-a/(a+b+c)
  qnew<-b/(a+b+c)

  return(c(pnew,qnew))
}
p<-0.1
q<-0.1
E<-NULL
for(i in 1:10000){
  p<-c(p,iter(p[i],q[i])[1])
  q<-c(q,iter(p[i],q[i])[2])
  if((p[i+1]-p[i]<10^(-5))&(q[i+1]-q[i]<10^(-5)))
    break
}
n<-1:length(p)
cat('the MLE of p,q is',c(p[length(p)],q[length(q)]))
plot(n,p,type="l",col="red")
lines(n,q,col="blue")

## -----------------------------------------------------------------------------
attach(mtcars)
formulas <- list(  
  mpg ~ disp,  
  mpg ~ I(1/disp),  
  mpg ~ disp + wt,  
  mpg ~ I(1 / disp) + wt
)

lapply(formulas,lm,data=mtcars)
for(i in 1:4){
  print(lm(formulas[[i]],data=mtcars))
}

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {  
rows <- sample(1:nrow(mtcars), rep = TRUE)  
mtcars[rows, ]  
})
lapply(bootstraps,lm,formula=mpg~disp)
for(i in 1:10){
  print(lm(mpg ~ disp, data = bootstraps[[i]]) )
}

## -----------------------------------------------------------------------------
attach(mtcars)
rsq <- function(mod) summary(mod)$r.squared
formulas <- list(  
  mpg ~ disp,  
  mpg ~ I(1/disp),  
  mpg ~ disp + wt,  
  mpg ~ I(1 / disp) + wt
)
bootstraps <- lapply(1:10, function(i) {  
rows <- sample(1:nrow(mtcars), rep = TRUE)  
mtcars[rows, ]  
})
list1<-lapply(formulas,lm)
list2<-lapply(bootstraps,lm,formula=mpg~disp)
lapply(list1,rsq)
lapply(list2,rsq)

## -----------------------------------------------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
sapply(trials, '[[', 3) 

## -----------------------------------------------------------------------------

mcsapply<-function(x,f,n){
    f<-match.fun(f)
    answer <- mclapply(x,f,mc.cores=n)
    if (USE.NAMES && is.character(x) && is.null(names(answer))) 
        names(answer) <- x
    if (!isFALSE(simplify) && length(answer)) 
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
##Rcpp method
cppFunction('List rw_Mc(double x0, double sigma, int N){
NumericVector x(N);
x[0]=x0;
int k=0;
for(int i=1;i<N;i++){

double y=as<double>(rnorm(1,x[i-1], sigma));
double u=as<double>(runif(1));

if (u<=exp(-abs(y))/exp(-abs(x[i-1]))) {
x[i]=y;
k=k+1;  
}
else x[i]=x[i-1];
}

return(List::create(Named("x")=x,Named("k")=k));
}')
## R method
rw.M<-function(x0,sigma,N){   ## generate sample for different sigma
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for(i in 2:N){
    y<-rnorm(1,mean=x[i-1],sd=sigma)
    if(u[i]<=exp(-abs(y))/exp(-abs(x[i-1]))) {
      x[i]<-y
      k<-k+1}
    else x[i]<-x[i-1]
  }
  return(list(x=x,k=k))
}
sigma<-c(0.05,0.5,2,16)
x0<-25
N<-2000
rw1<-rw.M(x0,sigma[1],N)  ## four samples
rw2<-rw.M(x0,sigma[2],N)
rw3<-rw.M(x0,sigma[3],N)
rw4<-rw.M(x0,sigma[4],N)
rwC1<-rw_Mc(x0,sigma[1],N)  ## four samples
rwC2<-rw_Mc(x0,sigma[2],N)
rwC3<-rw_Mc(x0,sigma[3],N)
rwC4<-rw_Mc(x0,sigma[4],N)

##comparison
qqplot(rw1$x,rwC1$x)
qqplot(rw2$x,rwC2$x)
qqplot(rw3$x,rwC4$x)
qqplot(rw4$x,rwC4$x)
##computation time
ts1<-microbenchmark(rw.M(x0,sigma[1],N),rw_Mc(x0,sigma[1],N))
ts2<-microbenchmark(rw.M(x0,sigma[2],N),rw_Mc(x0,sigma[2],N))
ts3<-microbenchmark(rw.M(x0,sigma[3],N),rw_Mc(x0,sigma[3],N))
ts4<-microbenchmark(rw.M(x0,sigma[4],N),rw_Mc(x0,sigma[4],N))
summary(ts1)[,c(1,3,5,6)]
summary(ts2)[,c(1,3,5,6)]
summary(ts3)[,c(1,3,5,6)]
summary(ts4)[,c(1,3,5,6)]

