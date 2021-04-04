########################### 
# heterogeneous

final.group1.hetero=function(y, x, tau=.5, lambdaset, penalty="SCAD", bict = 6,
                             maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                             a=3.7,intercept=T)
{
  
  dyn.load("cqpen.so")
  
  # Set penalty function
  if(penalty == "SCAD"){
    pentype = as.integer(0)  
  } else if (penalty == "MCP"){
    pentype = as.integer(1)
  } else{
    pentype = as.integer(2)
  }
  
  K = as.integer(length(y))
  p = as.integer( ncol(x[[1]]) )
  m = as.integer(length(tau))
  n = as.integer(length(y[[1]]))
  ll = as.integer(length(lambdaset))
  
  residualsx=matrix(0, n * ll, m * K)
  
  
  betax = matrix(as.double(0), ll, p * K * m)
  
  for(i in 1 : K)
  {
    
    residualsx[, (m * (i - 1) + 1) : (m * i)] = matrix(as.double(rep(y[[i]], ll * m)), n * ll, m)
    
  }
  
  y         = lapply(y,as.double)
  xdouble     = matrix(as.double(unlist(x)),n,p*K)
  tau       = as.double(tau)
  a         = as.double(a)
  eps       = as.double(eps)
  maxin     = as.integer(maxin)
  maxout=as.integer(maxout)
  lambdaset    = as.double(lambdaset)
  
  if(intercept == F)
  {
    
    out=.Fortran("pmqmlam",xdouble,betax,residualsx,
                 n,p,K,tau,m,a,eps,maxin,maxout,lambdaset,pentype,ll)
    
    
    betaset = out[[2]]
    
    betaset[abs(betaset)<coef.cutoff]=0
    
    cqbic = numeric(ll)
    
    for(tt in 1 : ll)
    {
      
      coefficients=list()
      
      for (i in 1:(m*K)) coefficients[[i]]=betaset[tt, ((i-1)*p+1):(i*p)]
      
      rho=numeric(m*K*n)
      
      for (i in 1:K)
      {
        for (j in 1:m) 
        {
          residuals=y[[i]]-x[[i]]%*%coefficients[[(i-1)*m+j]]
          rho[(m*n*(i-1)+n*(j-1)+1):(m*n*(i-1)+n*(j-1)+n)]=residuals*(tau[j]-(residuals<0))
        }
        
      }
      
      card = sum(Reduce('+',lapply(coefficients,function(x) abs(x[-1])))!=0)
      cqbic[tt]=log(sum(rho))+card*log(n)*log(p)/(2*bict*K*n)
      
    }
    
  } else {
    
    intvalx = matrix(as.double(0), ll, m * K) 
    
    
    out=.Fortran("pmqmlamint",xdouble,betax,intvalx,residualsx,
                 n,p,K,tau,m,a,eps,maxin,maxout,lambdaset,pentype,ll)
    
    
    
    betaset = out[[2]]
    intset = out[[3]]
    
    betaset[abs(betaset)<coef.cutoff]=0
    
    cqbic = numeric(ll)
    
    for(tt in 1 : ll)
    {
      
      coefficients=list()
      
      for (i in 1:(m*K)) coefficients[[i]]=c(intset[tt, i], betaset[tt, ((i-1)*p+1):(i*p)])
      
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
      
      card = sum(Reduce('+',lapply(coefficients,function(x) abs(x[-1])))!=0)
      cqbic[tt]=log(sum(rho))+card*log(n)*log(p)/(2*bict*K*n)
      
    }
    
  }
  
  
  bicvalue = cbind(lambdaset,cqbic)
  loc = which.min(cqbic)
  
  lambda.final = lambdaset[loc]
  
  res=list()
  
  for(i in 1:(m*K))
  {
    
    
    if(intercept == T) 
    {
      
      res[[i]] = c( intset[loc, i], betaset[loc, ((i-1)*p+1):(i*p)] )
      
    } else {
      
      res[[i]] = betaset[loc, ((i-1)*p+1):(i*p)]
      
    }
    
  }

  return(list("coefficients"=res, "lambda"=lambda.final, "cqbic"=bicvalue))
  
}

