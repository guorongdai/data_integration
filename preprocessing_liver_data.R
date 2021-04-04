#install.packages("BiocManager")
#BiocManager::install("mixOmics")

library(mixOmics)

data("liver.toxicity") 
dat=liver.toxicity
ind1=which(dat$treatment$Dose.Group<200)
ind2=which(dat$treatment$Dose.Group>200)
x1=dat$gene[ind1,]
x2=dat$gene[ind2,]
y1=log(dat$clinic$Cholesterol.mg.dL.[ind1],10)
y2=log(dat$clinic$Cholesterol.mg.dL.[ind2],10)


p=ncol(x1)

corr1=numeric(p)
corr2=numeric(p)
for(i in 1:p)
{
  corr1[i]=cor(y1,x1[,i])
  corr2[i]=cor(y2,x2[,i])
}

ind=unique(c(order(corr1,decreasing=T)[1:200],order(corr2,decreasing=T)[1:200]))
xx1=x1[,ind]
xx2=x2[,ind]

dat1=cbind(y1,as.matrix(xx1))
dat2=cbind(y2,as.matrix(xx2))

save(dat1,dat2,file="liver.RData")


