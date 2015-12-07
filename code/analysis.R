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

res <- shim(x = X, y = Y, main.effect.names = main_effect_names, 
            interaction.names = interaction_names,
            lambda.beta = 1, lambda.gamma = 1, threshold = 1e-4, max.iter = 500)

# plot of cofficients at each iteration
matplot(res$beta[,1:res$m] %>% t, type = "l", ylab="")
matplot(res$gamma[,1:res$m] %>% t, type = "l", ylab="")



