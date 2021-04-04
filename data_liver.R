.libPaths("/home/grad/rondai/RLibrary")
setwd("/home/grad/rondai/di/data_code")

suppressMessages({
  
  if (!require(doParallel)) install.packages("doParallel", 
                                             repos="http://cran.r-project.org")
  
  library(doParallel)
  
})


source("di_functions.R")

# Function to generate model size and prediction error
###############################
rda=function(coef,y,x,tau)
{
  k=length(y)
  m=length(tau)
  
  num=sum(Reduce('+',lapply(coef,abs))[-1]!=0)
  
  pe=0
  for(i in 1:k)
  {
    for(j in 1:m)
      resid=y[[i]]-x[[i]]%*%coef[[m*(i-1)+j]][-1]-coef[[1]][1]
    err=resid*(tau[j]-(resid<0))
    pe=pe+sum(err)
  }
  
  
  final=list("Model size"=num,"Prediction error"=pe)
  
  return(final)
  
}
###############################

registerDoParallel(detectCores())

load("liver.RData") # data
x=list(dat1[,-1],dat2[,-1])
y=list(dat1[,1],dat2[,1])

tau=seq(1/10,9/10,1/10)
m=length(tau)
n=nrow(x[[1]])

bict = 6

lambda.set=exp(seq(-4,-2,length=20))

f1=final.group1.hetero(y,x,tau=tau,lambdaset=lambda.set,penalty="SCAD",bict=bict)
coef1=f1$coefficients

f2=final.group1.hetero(y,x,tau=tau,lambdaset=lambda.set,penalty="MCP",bict=bict)
coef2=f2$coefficients


coef3=list()
for(i in 1:2)
{
  for(j in 1:m)
  {
    f3=final.group1.hetero(y[i],x[i],tau=tau[j],lambdaset=lambda.set,penalty="SCAD",bict=bict)
    coef3=c(coef3,f3$coefficients)
  }
}


coef4=list()
for(i in 1:2)
{
  for(j in 1:m)
  {
    f4=final.group1.hetero(y[i],x[i],tau=tau[j],lambdaset=lambda.set,penalty="MCP",bict=bict)
    coef4=c(coef4,f4$coefficients)
  }
}

print("##### Size of DI-SCAD ########")
sum(Reduce('+',lapply(coef1,abs))[-1]!=0) # Size of DI-SCAD
print("##### Size of DI-MCP ########")
sum(Reduce('+',lapply(coef2,abs))[-1]!=0) # Size of DI-MCP
print("##### Size of CA-SCAD ########")
sum(Reduce('+',lapply(coef3,abs))[-1]!=0) # Size of CA-SCAD
print("##### Size of CA-MCP ########")
sum(Reduce('+',lapply(coef4,abs))[-1]!=0) # Size of CA-MCP

print("##### DI-SCAD and DI-MCP ########")
sum((Reduce('+',lapply(coef1,abs))[-1]!=0) & (Reduce('+',lapply(coef2,abs))[-1]!=0)) # DI-SCAD and DI-MCP
print("##### CA-SCAD and CA-MCP ########")
sum((Reduce('+',lapply(coef3,abs))[-1]!=0) & (Reduce('+',lapply(coef4,abs))[-1]!=0)) # CA-SCAD and CA-MCP
print("##### DI-SCAD and CA-SCAD ########")
sum((Reduce('+',lapply(coef1,abs))[-1]!=0) & (Reduce('+',lapply(coef3,abs))[-1]!=0)) # DI-SCAD and CA-SCAD
print("##### DI-SCAD and CA-MCP ########")
sum((Reduce('+',lapply(coef1,abs))[-1]!=0) & (Reduce('+',lapply(coef4,abs))[-1]!=0)) # DI-SCAD and CA-MCP
print("##### DI-MCP and CA-SCAD ########")
sum((Reduce('+',lapply(coef2,abs))[-1]!=0) & (Reduce('+',lapply(coef3,abs))[-1]!=0)) # DI-MCP and CA-SCAD
print("##### DI-MCP and CA-MCP ########")
sum((Reduce('+',lapply(coef2,abs))[-1]!=0) & (Reduce('+',lapply(coef4,abs))[-1]!=0)) # DI-MCP and CA-MCP


out=numeric(8)
error=foreach(t=1:50) %dopar%
  {
    set.seed(t)
    
    ind=sample(1:n,24)
    
    x.t=list(dat1[ind,-1],dat2[ind,-1])
    y.t=list(dat1[ind,1],dat2[ind,1])
    
    x.v=list(dat1[-ind,-1],dat2[-ind,-1])
    y.v=list(dat1[-ind,1],dat2[-ind,1])
    
    f11=final.group1.hetero(y.t,x.t,tau=tau,lambdaset=lambda.set,penalty="SCAD")
    coef11=f11$coefficients
    
    f22=final.group1.hetero(y.t,x.t,tau=tau,lambdaset=lambda.set,penalty="MCP")
    coef22=f22$coefficients
    
    
    coef33=list()
    for(i in 1:2)
    {
      for(j in 1:m)
      {
        f33=final.group1.hetero(y.t[i],x.t[i],tau=tau[j],lambdaset=lambda.set,penalty="SCAD")
        coef33=c(coef33,f33$coefficients)
      }
    }
    
    
    coef44=list()
    for(i in 1:2)
    {
      for(j in 1:m)
      {
        f44=final.group1.hetero(y.t[i],x.t[i],tau=tau[j],lambdaset=lambda.set,penalty="MCP")
        coef44=c(coef44,f44$coefficients)
      }
    }
    
    out[1:2]=rda(coef11,y.v,x.v,tau)
    out[3:4]=rda(coef22,y.v,x.v,tau)
    out[5:6]=rda(coef33,y.v,x.v,tau)
    out[7:8]=rda(coef44,y.v,x.v,tau)
    
    out
  }

res=t(matrix(unlist(error),nrow=8,ncol=length(error)))
colnames(res)=c("Size of DI-SCAD","PE of DI-SCAD","Size of DI-MCP","PE of DI-MCP",
                "Size of CA-SCAD","PE of CA-SCAD","Size of CA-MCP","PE of CA-MCP")

print("#########################Mean#########################")
apply(res,2,mean)
print("#########################SD#########################")
apply(res,2,sd)






