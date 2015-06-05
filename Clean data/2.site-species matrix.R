## setwd("C:/Users/Hercules/McGill Stuff/NetworkStuff")
dat<-read.csv("finally_hooray.csv",header=TRUE, as.is=TRUE)

## create site-species matrix
ss.mat<-matrix(nrow=length(unique(dat$SITE_ID)),
               ncol=length(unique(dat$TAXANAME)))
ss.mat<-as.data.frame(ss.mat)
rownames(ss.mat)<-unique(dat$SITE_ID)
colnames(ss.mat)<-unique(dat$TAXANAME)
for(i in 1:nrow(ss.mat)){
  site<-rownames(ss.mat)[i]
  submat<-dat[which(dat$SITE_ID==site),]
  for(j in 1:nrow(submat)){
    spec<-as.character(submat[j,]$TAXANAME)
    ss.mat[i,which(colnames(ss.mat)==spec)]<-submat[j,]$abund_ml
  }
}


##  New code
library(reshape2)
dat<-read.csv("finally_hooray.csv",header=TRUE, as.is=TRUE)
aux = dcast(dat, SITE_ID ~ TAXANAME, value.var= 'abund_ml', fill =0, fun=mean)

## Analysis
## The dimensions do not match because a TAXANAME exist with ("")
dim(aux)
dim(ss.mat)
aux[1:5,1:5]
ss.mat[1:5,'Bosmina']
aux[1:5, 'Bosmina']
## There exist an TAXANAME that is empty (""), see below
unique(dat$TAXANAME)[which.min(nchar(unique(dat$TAXANAME)))]
