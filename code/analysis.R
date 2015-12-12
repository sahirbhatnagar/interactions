##################################
# R source code file to analyse the simulated data
# Created by Sahir, Dec 7, 2015
# Updated 
# hosted on Github repo 'sahirbhatnagar/interactions'
# NOTE: the shim function implements the Choi et al. algorithm
# for one pair of tuning parameters. I still need to work on the
# cross validation scheme for choosing the tuning parameters
##################################

rm(list = ls())
options(scipen=999, digits = 10)

source("packages.R")
source("data.R")
source("functions.R")



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

res[[1]]$Q[complete.cases(res[[1]]$Q),]

res[[1]]$beta[complete.cases(res[[1]]$beta)]


# plot of cofficients at each iteration
matplot(res$beta[,1:res$m] %>% t, type = "l", ylab="")
matplot(res$gamma[,1:res$m] %>% t, type = "l", ylab="")

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




