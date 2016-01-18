##################################
# R source code file to analyse the simulated data
# Created by Sahir, Dec 7, 2015
# Updated 
# hosted on Github repo 'sahirbhatnagar/interactions'
# NOTE: the shim function implements the Choi et al. algorithm
# for one pair of tuning parameters. I still need to work on the
# cross validation scheme for choosing the tuning parameters
##################################

dev.off()
rm(list = ls())
options(scipen=999, digits = 10)

source("packages.R")
source("data.R")
source("https://raw.githubusercontent.com/noamross/noamtools/master/R/proftable.R")
source("functions.R")
"%ni%" <- Negate("%in%")

system.time(res <- shim_multiple(x = X, y = Y, main.effect.names = main_effect_names,
                     interaction.names = interaction_names,
                     lambda.beta = NULL , lambda.gamma = NULL,
                     threshold = 1e-5 , max.iter = 100 , initialization.type = "ridge",
                     nlambda.gamma = 5, nlambda.beta = 10, cores = 3))

# 1 core is faster than more ...
system.time(res2 <- shim_multiple_faster(x = X, y = Y, main.effect.names = main_effect_names,
                     interaction.names = interaction_names,
                     lambda.beta = NULL , lambda.gamma = NULL,
                     threshold = 1e-5 , max.iter = 100 , initialization.type = "ridge",
                     nlambda.gamma = 5, nlambda.beta = 10, cores = 1))

# user defined lambda sequence
system.time(res <- shim_multiple(x = X, y = Y, main.effect.names = main_effect_names,
                     interaction.names = interaction_names,
                     lambda.beta = seq(0.1, 10, length.out = 5) , 
                     lambda.gamma = seq(0.1, 100, length.out = 5),
                     threshold = 1e-5 , max.iter = 100 , initialization.type = "ridge",
                     nlambda.gamma = 5, nlambda.beta = 5, cores = 8))

betas <- matrix(unlist(res2$beta), ncol = length(res2$beta), byrow = TRUE)
dim(betas)

gammas <- matrix(unlist(res2$gamma), ncol = length(res2$gamma), byrow = TRUE)
dim(gammas)

matplot(t(betas), type="l")
matplot(t(gammas), type="l")
length(res$beta)
matplot(res2$Q , type="l")

res2[[1]]


system.time(res2 <- parallel::mcmapply(shim, lambda.beta = seq(0.1, 10, length.out = 25), 
                          lambda.gamma = seq(0.1, 100, length.out = 25), 
                          MoreArgs = list(x = X, y = Y, main.effect.names = main_effect_names, 
                                          interaction.names = interaction_names,
                                          threshold = 1e-5, max.iter = 500, 
                                          initialization.type = "ridge"), SIMPLIFY = F, mc.cores = 8))

lapply(res2, function(i) i$m)

res2[[10]]



# Trying to estimate both betas and gammas --------------------------------

true.betas.and.alphas <- matrix(rep(0,55),nrow = 55, ncol=1) %>% 
  magrittr::set_rownames(colnames(X))
true.betas.and.alphas[names(beta5),] <- beta5
true.betas.and.gammas <- convert(true.betas.and.alphas, main_effect_names, interaction_names)

res <- shim(x = X, y = Y, main.effect.names = main_effect_names, 
            interaction.names = interaction_names,
            lambda.beta = 1, lambda.gamma = 0.5, threshold = 1e-5, max.iter = 500, 
            initialization.type = "ridge")


library(doParallel)
library(foreach)
library(parallel)
options("mc.cores" = 10L)
getOption("mc.cores", 5L)

parallel::detectCores()

res <- parallel::mcmapply(shim, lambda.beta = seq(0,10,1), lambda.gamma = seq(0,10,1), 
                          MoreArgs = list(x = X, y = Y, main.effect.names = main_effect_names, 
                                          interaction.names = interaction_names,
                                          threshold = 1e-5, max.iter = 500, 
                                          initialization.type = "ridge"), SIMPLIFY = F)


res2 <- shim(lambda.beta = 1.5, lambda.gamma = 2, 
x = X, y = Y, main.effect.names = main_effect_names, 
                                          interaction.names = interaction_names,
                                          threshold = 1e-5, max.iter = 500, 
                                          initialization.type = "ridge")

res$Q[complete.cases(res$Q),]

res2$beta[complete.cases(res2$beta)]


# plot of cofficients at each iteration
matplot(res2$beta[,1:res2$m] %>% t, type = "l", ylab="")
matplot(res2$gamma[,1:res2$m] %>% t, type = "l", ylab="")

cbind2(round(res$beta[,1:res$m],2), true.betas.and.gammas[main_effect_names,,drop=F])
cbind2(round(res$gamma[,1:res$m],2), true.betas.and.gammas[interaction_names,,drop=F])




# Trying to estimate gammas with fixed beta -------------------------------



res <- shim_fix_betas(x = X, y = Y, main.effect.names = main_effect_names, 
                      interaction.names = interaction_names,
                      lambda.beta = 0.5, lambda.gamma = 0.5, threshold = 1e-8, max.iter = 500,
                      fixed.beta = true.betas.and.gammas[main_effect_names,, drop = F])

# plot of cofficients at each iteration
matplot(res$beta[,1:res$m] %>% t, type = "l", ylab="")
matplot(res$gamma[,1:res$m] %>% t, type = "l", ylab="")
cbind2(res$gamma[,1:res$m], true.betas.and.gammas[interaction_names,,drop=F])


# Trying to estimate betas with fixed gamma -------------------------------

res <- shim_fix_gamma(x = X, y = Y, main.effect.names = main_effect_names, 
                      interaction.names = interaction_names,
                      lambda.beta = 0.5, lambda.gamma = 0.5, threshold = 1e-8, 
                      max.iter = 500, initialization.type = "ridge",
                      fixed.gamma = true.betas.and.gammas[interaction_names, , drop=F])

# plot of cofficients at each iteration
matplot(res$beta[,1:res$m] %>% t, type = "l", ylab="")
matplot(res$gamma[,1:res$m] %>% t, type = "l", ylab="")

cbind2(res$beta[,1:res$m], true.betas.and.gammas[main_effect_names,,drop=F])

plot(res$Q[1:res$m,], type="l")




