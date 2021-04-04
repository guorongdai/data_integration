.libPaths("/home/grad/rondai/RLibrary")

suppressMessages({
  
  if (!require(doParallel)) install.packages("doParallel", 
                                             repos="http://cran.r-project.org")
  if (!require(PMA)) install.packages('PMA', repos =
                                        'https://cran.revolutionanalytics.com/')
  
  if (!require(grpreg)) install.packages('grpreg', repos =
                                           'https://cran.revolutionanalytics.com/')
  
  if (!require(writexl)) install.packages('writexl', repos =
                                            'https://cran.revolutionanalytics.com/')
  
  if (!require(magic)) install.packages('magic', repos =
                                          'https://cran.revolutionanalytics.com/')
  
  if (!require(rqPen)) install.packages("rqPen", 
                                        repos="http://cran.r-project.org")
  if (!require(FusionLearn)) install.packages("FusionLearn", 
                                              repos="http://cran.r-project.org")
  if (!require(quantreg)) install.packages("quantreg", 
                                           repos="http://cran.r-project.org")
  
  
  library(MASS)
  library(PMA)
  library(Matrix)
  library(doParallel)
  library(grpreg)
  library(writexl)
  library(magic)
  library(rqPen)
  library(FusionLearn)
  library(quantreg)
  
})

# Functions
########################### 
final.group1.hetero=function(y, x, tau=.5, lambdaset, penalty="SCAD", 
                             initial_beta=NULL, maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                             a=3.7,intercept=T)
{
  if(intercept==1) 
  {
    models=sapply(lambdaset,myQICD.group1.hetero,y=y,x=x,tau=tau,penalty=penalty)
  } else
  {
    models=sapply(lambdaset,myQICD.group1.hetero.no.inter,y=y,x=x,tau=tau,penalty=penalty)
  }
  res=models[1,which.min(models[2,])][[1]]
  lambda.final=lambdaset[which.min(models[2,])]
  cqbic=cbind(lambdaset,models[2,])
  return(list("coefficients"=res,
              "lambda"=lambda.final,"cqbic"=cqbic))
}


myQICD.group1.hetero = function(y, x, tau=.5, lambda, penalty="SCAD", 
                                initial_beta=NULL, maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                                a=3.7, ...)
{
  #cleanInputs(y, x, lambda, initial_beta, penalty, a)
  dyn.load("cqpen.so")
  require(rqPen)
  
  # Set penalty function
  if(penalty == "SCAD"){
    pentype = as.integer(0)  
  } else if (penalty == "MCP"){
    pentype = as.integer(1)
  } else{
    pentype = as.integer(2)
  }
  
  K         = as.integer(length(y))
  p         = as.integer( ncol(x[[1]]) )
  m         = as.integer(length(tau))
  n         = as.integer(length(y[[1]]))
  
  beta=numeric(K*p*m)
  intval=numeric(m*K)
  residuals=numeric(n*m*K)
  
  if( is.null(initial_beta) ){
    for(i in 1:K)
    {
      for(j in 1:m)
      {
        coef=LASSO.fit(y=y[[i]],x=x[[i]],tau=tau[j],lambda=lambda,
                       coef.cutoff=coef.cutoff,intercept=1)
        beta[(p*m*(i-1)+p*(j-1)+1):(p*m*(i-1)+p*j)]=coef[-1]
        intval[m*(i-1)+j]=coef[1]
        residuals[(n*m*(i-1)+n*(j-1)+1):(n*m*(i-1)+n*j)]=y[[i]]-x[[i]]%*%coef[-1]-coef[1]
      }
    }  
  }
  
  
  
  
  y         = lapply(y,as.double)
  xdouble     = matrix(as.double(unlist(x)),n,p*K)
  tau       = as.double(tau)
  a         = as.double(a)
  eps       = as.double(eps)
  maxin     = as.integer(maxin)
  maxout=as.integer(maxout)
  betadouble=as.double(beta)
  intvaldouble=as.double(unlist(intval))
  residualsdouble=matrix(as.double(residuals),n,m*K)
  lambda    = as.double(lambda)
  
  # groupl1 = rep(0, p)
  
  out=.Fortran("multicqpenhetero",xdouble,betadouble,intvaldouble,residualsdouble,
               n,p,K,tau,m,a,eps,maxin,maxout,lambda,pentype)
  
  out[[2]][abs(out[[2]])<coef.cutoff]=0
  
  coefficients=list()
  
  for (i in 1:(m*K)) coefficients[[i]]=c(out[[3]][i],out[[2]][((i-1)*p+1):(i*p)])
  
  rho=numeric(m*K*n)
  
  for (i in 1:K)
  {
    for (j in 1:m) 
    {
      residuals=y[[i]]-x[[i]]%*%coefficients[[(i-1)*m+j]][-1]-
        coefficients[[(i-1)*m+j]][1]
      rho[(m*n*(i-1)+n*(j-1)+1):(m*n*(i-1)+n*(j-1)+n)]=residuals*(tau[j]-(residuals<0))
    }
    
  }
  
  cqbic=log(sum(rho))+sum(Reduce('+',lapply(coefficients,function(x) abs(x[-1])))!=0)*log(n)*log(p)/(2*n)
  
  return( list(coefficients,cqbic) )
}
########################### 

data("stockindexVIX")
data("stockindexGSPC")
data("stockindexDJI")
data("validVIX")
data("validGSPC")
data("validDJI")

valVIX=validDJI
valDJI=validVIX
valGSPC=validGSPC

# Function generating model size and prediction error
###############################
rda=function(coef,n,m, tau)
{
  num=sum(Reduce('+',lapply(coef,abs))[-1]!=0)
  
  residVIX=sapply(coef[1:m],function(y,x,beta) y-x%*%beta[-1]-beta[1],
                  y=valVIX[,1],x=valVIX[,-1])
  residDJI=sapply(coef[(m+1):(2*m)],function(y,x,beta) y-x%*%beta[-1]-beta[1],
                  y=valDJI[,1],x=valDJI[,-1])
  residGSPC=sapply(coef[(2*m+1):(3*m)],function(y,x,beta) y-x%*%beta[-1]-beta[1],
                   y=valGSPC[,1],x=valGSPC[,-1])
  
  
  
  seVIX=matrix(0,n,m)
  seDJI=matrix(0,n,m)
  seGSPC=matrix(0,n,m)
  for (i in 1:m)
  {
    seVIX[,i]=residVIX[,i]*(tau[i]-(residVIX[,i]<0))
    seDJI[,i]=residDJI[,i]*(tau[i]-(residDJI[,i]<0))
    seGSPC[,i]=residGSPC[,i]*(tau[i]-(residGSPC[,i]<0))
  }
  sseVIX=sum(seVIX)
  sseDJI=sum(seDJI)
  sseGSPC=sum(seGSPC)
  
  final=list("Model size"=num,"sseVIX"=sseVIX,"sseDJI"=sseDJI,"sseGSPC"=sseGSPC)
  
  return(final)
  
}

rda.median=function(coef,n,m, tau, factor)
{
  num=sum(Reduce('+',lapply(coef,abs))[-1]!=0)
  
  residVIX=sapply(coef[1:m],function(y,x,beta) y-x%*%beta[-1]-beta[1],
                  y=valVIX[,1],x=valVIX[,-1])
  residDJI=sapply(coef[(m+1):(2*m)],function(y,x,beta) y-x%*%beta[-1]-beta[1],
                  y=valDJI[,1],x=valDJI[,-1])
  residGSPC=sapply(coef[(2*m+1):(3*m)],function(y,x,beta) y-x%*%beta[-1]-beta[1],
                   y=valGSPC[,1],x=valGSPC[,-1])
  
  
  
  seVIX=matrix(0,n,m)
  seDJI=matrix(0,n,m)
  seGSPC=matrix(0,n,m)
  for (i in 1:m)
  {
    seVIX[,i]=residVIX[,i]*(tau[i]-(residVIX[,i]<0))
    seDJI[,i]=residDJI[,i]*(tau[i]-(residDJI[,i]<0))
    seGSPC[,i]=residGSPC[,i]*(tau[i]-(residGSPC[,i]<0))
  }
  sseVIX=sum(seVIX)
  sseDJI=sum(seDJI)
  sseGSPC=sum(seGSPC)
  
  final=list("Model size"=num,"sseVIX"=sseVIX * factor,"sseDJI"=sseDJI * factor,
             "sseGSPC"=sseGSPC * factor)
  
  return(final)
  
}
###############################

y=list("VIX"=stockindexVIX[,1],"GSPC"=stockindexGSPC[,1],"DJI"=stockindexDJI[,1])
x=list(stockindexVIX[,-1],stockindexGSPC[,-1],stockindexDJI[,-1])
tau=seq(1/20,19/20,by=1/20)
m=length(tau)
n=nrow(valVIX)
p = ncol(valVIX) - 1
n0 = nrow(stockindexDJI)
K = 3

lambda.set=exp(seq(-2.9,-0.8,length=40))

f1=final.group1.hetero(y,x,tau=tau,lambdaset=lambda.set,penalty="SCAD")
coef1=f1$coefficients

f2=final.group1.hetero(y,x,tau=tau,lambdaset=lambda.set,penalty="MCP")
coef2=f2$coefficients

coef3=list()
for(i in 1:3)
{
  for(j in 1:m)
  {
    f3=final.group1.hetero(y[i],x[i],tau=tau[j],lambdaset=lambda.set,penalty="SCAD")
    coef3=c(coef3,f3$coefficients)
  }
}

coef4=list()
for(i in 1:3)
{
  for(j in 1:m)
  {
    f4=final.group1.hetero(y[i],x[i],tau=tau[j],lambdaset=lambda.set,penalty="MCP")
    coef4=c(coef4,f4$coefficients)
  }
}

f5=final.group1.hetero(y,x,tau=tau,lambdaset=0,penalty="SCAD")
coef5=f5$coefficients

xmat = cbind(x[[1]], x[[2]], x[[3]])
ymat = cbind(y[[1]], y[[2]], y[[3]])

set.seed(1)
f6.cv = CCA.permute(x = xmat, z = ymat, typex = "standard",
                    typez = "standard", trace = F)
f6 = CCA(x = xmat, z = ymat, typex = "standard", typez = "standard",
         trace = F) $ u
ind = which((f6[1 : p] != 0) | (f6[(p + 1) : (2 * p)]) | (f6[(1 + 2 * p) : (3 * p)]))
coef6 = rep( list(numeric(p + 1)), K )
for(i in 1 : K) coef6[[i]][c(1, ind + 1)] = rq(y[[i]] ~ x[[i]][, ind]) $ coefficients



xx = adiag( cbind(rep(1,n0), x[[1]]), cbind(rep(1, n0), x[[2]]), cbind(rep(1, n0), x[[3]]) )
yy = c(y[[1]], y[[2]], y[[3]])
group = as.integer(c(0, 1 : p, 0, 1 : p, 0, 1 : p))
lambda = cv.grpreg(X = xx, y = yy, group = group, penalty = "grSCAD",
                   family = "gaussian", seed = 1) $ lambda.min
f7 = (grpreg(X = xx, y = yy, group = group,
             penalty = "grSCAD", family = "gaussian", lambda = lambda) $ beta)[-1]
coef7 = list(f7[1 : (p + 1)], f7[(p + 2) : (2 + 2 * p)], f7[(3 + 2 * p) : (3 + 3 * p)])


print("##### DI-SCAD and DI-MCP ########")
sum((Reduce('+',lapply(coef1,abs))[-1]!=0) & (Reduce('+',lapply(coef2,abs))[-1]!=0)) # DI-SCAD and DI-MCP
print("##### DI-SCAD and CA-SCAD ########")
sum((Reduce('+',lapply(coef1,abs))[-1]!=0) & (Reduce('+',lapply(coef3,abs))[-1]!=0)) # DI-SCAD and CA-SCAD
print("##### DI-SCAD and CA-MCP ########")
sum((Reduce('+',lapply(coef1,abs))[-1]!=0) & (Reduce('+',lapply(coef4,abs))[-1]!=0)) # DI-SCAD and CA-MCP
print("##### DI-SCAD and CCA ########")
sum((Reduce('+',lapply(coef1,abs))[-1]!=0) & (Reduce('+',lapply(coef6,abs))[-1]!=0)) # DI-SCAD and CCA
print("##### DI-SCAD and LS ########")
sum((Reduce('+',lapply(coef1,abs))[-1]!=0) & (Reduce('+',lapply(coef7,abs))[-1]!=0)) # DI-SCAD and LS

rda(coef1,n,m, tau)
rda(coef2,n,m, tau)
rda(coef3,n,m, tau)
rda(coef4,n,m, tau)
rda(coef5,n,m, tau)
rda.median(coef6,n,1, 0.5, m)
rda.median(coef7,n,1, 0.5, m)





