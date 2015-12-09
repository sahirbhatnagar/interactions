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
#source("packages.R")

library(hierNet)
library(magrittr)

"%ni%" <- Negate("%in%")

set.seed(12)

# number of predictors
p = 10

# number of subjects
n = 200

# predictors
x = matrix(rnorm(n*p),ncol=p)

# environment
e = matrix(rbinom(n,1,0.5), ncol=1)

# create the zz matrix
X <- hierNet::compute.interactions.c(cbind2(x,e)) 

# find columnames corresponding to main effects and interactions with 'e'
(columns.to.keep <- c(paste0(1:p, ":",1:p), X %>% 
                        colnames() %>% 
                        grep((p+1) %>% as.character(),., value = T)))

# set all other pairwise interactions to 0
X[,colnames(X) %ni% columns.to.keep ] <- 0

head(X)

# generate response
y = X[,"1:1"]+ X[,"5:5"] + X[,"11:11"] + 1.5*X[,"1:11"] + 2*X[,"5:11"]+3*rnorm(n)

fit = hierNet.path(cbind2(x,e),y, zz = X)
fitcv = hierNet.cv(fit,cbind2(x,e),y)
print(fitcv)
plot(fitcv)

# matrix of estimated interaction coefficients 
(d <- fit$th[,,which(fitcv$lamhat.1se==fit$lamlist)])

# main effect estimates
fit$bp[,which(fitcv$lamhat.1se==fit$lamlist), drop = F] - fit$bn[,which(fitcv$lamhat.1se==fit$lamlist), drop = F]


