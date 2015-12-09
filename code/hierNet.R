##################################
# R source code file for trying to trick the hierNet package
# into doing GxE interactions
# Created by Sahir, Dec 7, 2015
# Updated 
# hosted on Github repo 'sahirbhatnagar/interactions'
# NOTE: this script will automatically install packages not
# currently available locally. source 'packages.R' before running
##################################

rm(list = ls())
source("packages.R")
"%ni%" <- Negate("%in%")

set.seed(12)
p = 100
n = 200
x=matrix(rnorm(n*p),ncol=p)
e=matrix(rbinom(n,1,0.5), ncol=1)



X <- hierNet::compute.interactions.c(cbind2(x,e)) 
(columns.to.keep <- c(paste0(1:p, ":",1:p), X %>% colnames() %>% grep((p+1) %>% as.character(),., value = T)))
X[,colnames(X) %ni% columns.to.keep ] <- 0

head(X)

#X=scale(X[,columns.to.keep],TRUE,TRUE)

y=X[,"1:1"]+ X[,"5:5"] + X[,"101:101"] + 1.5*X[,"1:101"] + 1.5*X[,"5:101"]+1.5*rnorm(n)
fit=hierNet.path(cbind2(x,e),y, zz = X)
fitcv=hierNet.cv(fit,cbind2(x,e),y)
print(fitcv)
plot(fitcv)


fitcv

(d <- fit$th[,,which(fitcv$lamhat.1se==fit$lamlist)])

fit$bp[,which(fitcv$lamhat.1se==fit$lamlist), drop = F] - fit$bn[,which(fitcv$lamhat.1se==fit$lamlist), drop = F]


