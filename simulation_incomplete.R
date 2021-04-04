.libPaths("/home/grad/rondai/RLibrary")

suppressMessages({
  if (!require(doParallel)) install.packages('doParallel', repos =
                                               'https://cran.revolutionanalytics.com/')
  
  if (!require(rqPen)) install.packages('rqPen', repos =
                                          'https://cran.revolutionanalytics.com/')
  
  if (!require(mvtnorm)) install.packages('mvtnorm', repos =
                                            'https://cran.revolutionanalytics.com/')
  
  if (!require(writexl)) install.packages('writexl', repos =
                                            'https://cran.revolutionanalytics.com/')
  
  
  library(MASS)
  library(rqPen)
  library(Matrix)
  library(doParallel)
  library(mvtnorm)
  library(writexl)
})


source("/home/grad/rondai/di/simulation_code/di_functions.R")


# Incomplete grouping

registerDoParallel(detectCores())

K=2 # number of experiments
n=100 # sample size
p=100 # number of parameters
#p=200
s=100 # number of repetitions
tau=seq(1/6,5/6,1/6)
m=length(tau)
indc1=1
indc2=3
coe=0.7
bict=3

# covariance matrix of X_{k..}
cov.x=matrix(0,p,p)
for(i in 1:p)
{
  for (j in 1:p) cov.x[i,j]=0.5^(abs(i-j))
}

# covariance matrix of Y
cov.epsilon=matrix(0.7,K,K)
diag(cov.epsilon)=1

# generate theta and store it in a p*K-dimensional vector
theta0=list(rep(0,p),rep(0,p))

ind1=c(4,6,9,12,15,20)
ind2=c(1,6,12,15,20,25)
ind=sort(unique(c(ind1,ind2)))
set.seed(999)
theta0[[1]][ind1]=runif(length(ind1),0.05,1)
theta0[[2]][ind2]=runif(length(ind2),0.05,1)

inda1=sort(c(ind1,indc1))
inda2=sort(c(ind2,indc2))
inda=sort(unique(c(inda1,inda2)))

q=length(inda) # number of non-zero parameters 


#1(1) and 2(3) are correlated with the error
thetatrue=numeric(K*m*p)
for(i in 1:K) 
{
  for(j in 1:m)
  {
    thetatrue[((i-1)*m*p+(j-1)*p+1):((i-1)*m*p+j*p)]=theta0[[i]]
  }
}
for(j in 1:m) 
{
  thetatrue[(j-1)*p+indc1]=qnorm(tau[j])*coe
  thetatrue[m*p+(j-1)*p+indc2]=qnorm(tau[j])*coe
}

thetatrue1=thetatrue[c((p+1):(2*p),((m+1)*p+1):((m+2)*p))]
thetatrue2=thetatrue[c((2*p+1):(3*p),((m+2)*p+1):((m+3)*p))]

r=length(inda1)+length(inda2)

out=numeric(9)

lambda.set=exp(seq(-3.5,-1.5,length=30))


error=foreach(i=1:s, .packages = c("MASS", "rqPen", "Matrix", "mvtnorm")) %dopar%
{
  # generate X and store it in a n*p*K array
  set.seed(i*777)
  
  x=replicate(K,mvrnorm(n,rep(0,p),cov.x),simplify=F)
  x[[1]][,indc1]=pnorm(x[[1]][,indc1])
  x[[2]][,indc2]=pnorm(x[[2]][,indc2])
  
  # error distribution
  epsilon=mvrnorm(n,rep(0,K),cov.epsilon)
  # epsilon=rmvt(n,sigma=cov.epsilon,df=3)
  
  epsilon[,1]=epsilon[,1]*x[[1]][,indc1]*coe
  epsilon[,2]=epsilon[,2]*x[[2]][,indc2]*coe
  y=mapply(function(x,y) x%*%y,x=x,y=theta0)+epsilon
  y=lapply(1:K,function(i,y) y[,i],y=y)
  
  f1=final.group1.hetero(y=y,x=x,tau=tau,lambdaset=lambda.set,penalty="SCAD",bict=bict,intercept=F)
  
  mod1=f1$coefficients
  
  mod2=list()
  for(t in 1:K)
  {
    for(j in 1:m)
    {
      f2=final.group1.hetero(y=list(y[[t]]),x=list(x[[t]]),tau=tau[j],
                             lambdaset=lambda.set,penalty="SCAD",bict=bict,intercept=F) 
      mod2=c(mod2,f2$coefficients)
    }
  }
  
  ca1=mod2[c(2,2+m)]
  ca2=mod2[c(3,3+m)]
  
  out[1]=sum(Reduce('+',lapply(mod1,abs))[inda]!=0)/q
  out[2]=sum(Reduce('+',lapply(ca1,abs))[inda]!=0)/q
  out[3]=sum(Reduce('+',lapply(ca2,abs))[inda]!=0)/q
  
  
  out[4]=sum(Reduce('+',lapply(mod1,abs))[-inda]!=0)/(p-q)
  out[5]=sum(Reduce('+',lapply(ca1,abs))[-inda]!=0)/(p-q)
  out[6]=sum(Reduce('+',lapply(ca2,abs))[-inda]!=0)/(p-q)
  
  out[7]=sum(abs(unlist(mod1)-thetatrue))/(m*K)
  out[8]=sum(abs(unlist(ca1)-thetatrue1))/(K)
  out[9]=sum(abs(unlist(ca2)-thetatrue2))/(K)
  
  
  out
}

res=t(matrix(unlist(error),nrow=9,ncol=s))

res[, 1 : 6] = res[, 1 : 6] * 100

colnames(res)=c("PSR of DI","PSR of CA1","PSR of CA2",
                "FDR of DI","FDR of CA1","FDR of CA2",
                "SE of DI","SE of CA1","SE of CA2")
ave = apply(res,2,mean)
std = apply(res,2,sd)

ta = matrix(0, 3, 6)

for(i in 1 : 3)
{
  
  ta[, (i - 1) * 2 + 1] = ave[(1 : 3) + 3 * (i - 1)]
  ta[, (i - 1) * 2 + 2] = std[(1 : 3) + 3 * (i - 1)]
  
}

tabular = list("tabular" = data.frame(ta))

write_xlsx(tabular, path = "/home/grad/rondai/di/incomplete.xlsx")

ave

std





