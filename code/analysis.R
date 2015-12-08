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

source("packages.R")
source("data.R")
source("functions.R")



true.betas.and.alphas <- matrix(rep(0,55),nrow = 55, ncol=1) %>% magrittr::set_rownames(colnames(X))
true.betas.and.alphas[names(beta4),] <- beta4
true.betas.and.alphas

data.frame(uni_fun(colnames(X), X, Y, type = "univariate"), 
           uni_fun(colnames(X), X, Y, type = "ridge"),
           true.betas.and.alphas)







res <- shim(x = X, y = Y, main.effect.names = main_effect_names, 
            interaction.names = interaction_names,
            lambda.beta = 2.5, lambda.gamma = 2.5, threshold = 1e-8, max.iter = 500)






# plot of cofficients at each iteration
matplot(res$beta[,1:res$m] %>% t, type = "l", ylab="")
matplot(res$gamma[,1:res$m] %>% t, type = "l", ylab="")

cbind2(res$beta[,1:res$m], true.betas.and.gammas[main_effect_names,,drop=F])
cbind2(res$gamma[,1:res$m], true.betas.and.gammas[interaction_names,,drop=F])





true.betas.and.gammas <- convert(true.betas.and.alphas, main_effect_names, interaction_names)

res <- shim_fix_betas(x = X, y = Y, main.effect.names = main_effect_names, 
            interaction.names = interaction_names,
            lambda.beta = 0.5, lambda.gamma = 0.5, threshold = 1e-8, max.iter = 500)

# plot of cofficients at each iteration
matplot(res$beta[,1:res$m] %>% t, type = "l", ylab="")
matplot(res$gamma[,1:res$m] %>% t, type = "l", ylab="")
cbind2(res$gamma[,1:res$m], true.betas.and.gammas[interaction_names,,drop=F])





res <- shim_fix_gamma(x = X, y = Y, main.effect.names = main_effect_names, 
                      interaction.names = interaction_names,
                      lambda.beta = 1.5, lambda.gamma = 1.5, threshold = 1e-8, max.iter = 500)

# plot of cofficients at each iteration
matplot(res$beta[,1:res$m] %>% t, type = "l", ylab="")
matplot(res$gamma[,1:res$m] %>% t, type = "l", ylab="")
cbind2(res$beta[,1:res$m], true.betas.and.gammas[main_effect_names,,drop=F])





