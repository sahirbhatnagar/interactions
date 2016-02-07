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
options(scipen=999, digits = 9)

source("packages.R")
source("data.R")
source("https://raw.githubusercontent.com/noamross/noamtools/master/R/proftable.R")
#source("~/Dropbox/Winter 2014/MATH 783/Assignments/A3/multiplot.R")
source("functions.R")
library(ggrepel)
# library(latex2exp)
# library(gridExtra)
library(cowplot)



# crossvalidation ---------------------------------------------------------

require(doMC)
registerDoMC(cores=10)

system.time(cvfit <- cv.shim(x = X, y = Y, main.effect.names = main_effect_names,
                   interaction.names = interaction_names,
                   lambda.beta = NULL , lambda.gamma = NULL,
                   threshold = 1e-4 , max.iter = 200 , initialization.type = "ridge",
                   intercept = TRUE, normalize = TRUE,
                   nlambda.gamma = 5, nlambda.beta = 20, cores = 1,
                   type.measure = c("mse", "deviance", "class", "auc", "mae"), 
                   nfolds = 10, grouped = TRUE, keep = FALSE, parallel = TRUE))


plot(cvfit)
coef(cvfit)
cvfit$lambda.min.beta


cvfit$glmnet.fit$converged
cvfit$glmnet.fit$beta[,which(cvfit$glmnet.fit$converged==1)]
cvfit$glmnet.fit$alpha[,which(cvfit$glmnet.fit$converged==1)]
cvfit$glmnet.fit$b0[which(cvfit$glmnet.fit$converged==1)]
cvfit$glmnet.fit$nlambda.gamma

cvfit$lambda.beta[which(cvfit$glmnet.fit$converged==1)] %>% unique() %>% length
cvfit$lambda.gamma[which(cvfit$glmnet.fit$converged==1)] %>% unique %>% length()

which(cvfit$lambda.beta==cvfit$lambda.1se.beta)
which(cvfit$lambda.gamma==cvfit$lambda.1se.gamma)[1]

(cvfit$converged!=10) %>% sum

coef(cvfit, s.beta = "s76.beta", s.gamma = "s76.gamma")
coef(cvfit, s.beta = "lambda.min.beta", s.gamma = "lambda.min.gamma")
coef(cvfit, s.beta = "lambda.1se.beta", s.gamma = "lambda.1se.gamma")
names(cvfit)
plot(cvfit$cvm[cvfit$cvm<1e4])
length(cvfit$cvm)
d <- coef(cvfit$glmnet.fit) %>% as.matrix()
apply(coef(cvfit$glmnet.fit), 2, nonzero)

# plot for Celia  ---------------------------------------------------------

# 1 core is faster than more ...
system.time(res2 <- shim_multiple_faster(x = X, y = Y, main.effect.names = main_effect_names,
                     interaction.names = interaction_names,
                     lambda.beta = NULL , lambda.gamma = NULL,
                     threshold = 1e-4 , max.iter = 500 , initialization.type = "ridge",
                     nlambda.gamma = 20, nlambda.beta = 10, cores = 10))

DT_norm <- as.data.table(t(rbind(as.matrix(res2$beta), as.matrix(res2$alpha))))
colnames(DT_norm)
DT_norm[, lambda:=paste0("s",1:nrow(DT_norm))]
DT_norm[,index:= 1:nrow(DT_norm)]
DT_norm[, lambda.beta:=res2$lambda.beta]
DT_norm[, lambda.gamma:=res2$lambda.gamma]
setkey(DT_norm, lambda)
DT_norm[, L1beta := rowSums(abs(.SD)), .SDcols = main_effect_names]
DT_norm[, L1alpha := rowSums(abs(.SD)), .SDcols = interaction_names]
#DT_norm[, L1gamma := rowSums(abs(.SD)), .SDcols = interaction_names]
DT_norm_melt <- melt(DT_norm, id.vars = c("lambda", "L1beta", "L1alpha","index", "lambda.beta", "lambda.gamma"))

DT_norm_melt[variable %in% c("x1","x2","x1:x2"), group:="x1:x2"] 
DT_norm_melt[variable %in% c("x3","x4","x3:x4"), group:="x3:x4"]
DT_norm_melt[, table(group, exclude = F)]


a <- ggplot(DT_norm_melt[!is.na(group)][variable %in% c("x1","x2", "x3","x4")][lambda %in% paste0("s",41:60)], aes(x = index, y = value, group = variable, color = group)) + 
    geom_line(size=1) + 
    coord_cartesian(xlim = c(40, 60 + 2)) +
    geom_text_repel(data = dplyr::distinct(DT_norm_melt[!is.na(group)][variable %in% c("x1","x2", "x3","x4")][lambda %in% paste0("s",41:60)][index==60], variable),
                                        aes(label = variable),
                                        size = 4,
                                        nudge_x = 0,
                                        nudge_y = 0.0,
                                        segment.color = NA) + 
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "tuning parameter index", y = "main effect")+ 
    scale_x_continuous(breaks = seq(40,60,5), labels = seq(0,20,5))

b <- ggplot(DT_norm_melt[!is.na(group)][variable %in% c("x1:x2","x3:x4")][lambda %in% paste0("s",41:60)], aes(x = index, y = value, group = variable, color = group)) + 
    geom_line(size=1) + 
    coord_cartesian(xlim = c(40, 60 + 2)) +
    geom_text_repel(data = dplyr::distinct(DT_norm_melt[!is.na(group)][variable %in% c("x1:x2", "x3:x4")][lambda %in% paste0("s",41:60)][index==60], variable),
                    aes(label = variable),
                    size = 4,
                    nudge_x = 0,
                    nudge_y = 0.0,
                    segment.color = NA) + 
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "tuning parameter index", y = "interaction effect") + 
    scale_x_continuous(breaks = seq(40,60,5), labels = seq(0,20,5))

final_tuning_plot <- cowplot::plot_grid(a,b, labels = c("A", "B"), align = "v", nrow =  2)
ggsave("heredity_plot.png")

# old graphs stuff --------------------------------------------------------


matplot(t(res2$beta)[41:60,true_var_names[1:4]], type = "l", lty = 1)
matplot(t(res2$gamma)[41:60,c("x1:x2","x3:x4")], type = "l", lty = 1)
dev.off()
par(mfrow=c(2,1))
index1 = apply(abs(res2$beta), 2, sum) 
index2 = apply(abs(res2$alpha), 2, sum) 
index3 = apply(abs(res2$gamma), 2, sum) 
matplot(index1[order(index1)], t(res2$beta)[order(index1),true_var_names[1:4]] , lty = 1, xlab = "L1 norm", ylab = "Coefficients", 
        type = "l", xlim = c(0,40))
matplot(index2[order(index2)], t(res2$alpha)[order(index2),c("x1:x2","x3:x4")], lty = 1, xlab = "L1 norm", ylab = "Coefficients", 
        type = "l", xlim = c(0,50))
matplot(index3[order(index3)], t(res2$gamma)[order(index3),true_var_names[5:10]], lty = 1, xlab = "L1 norm", ylab = "Coefficients", 
        type = "l", xlim = c(0,4))

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

# 1 core is faster than more ...
lam <- lambda_sequence(X,Y, ridge_weights(X,Y, main_effect_names, interaction_names), nlambda = 5)
lam
system.time(res2 <- shim_multiple_faster(x = X, y = Y, main.effect.names = main_effect_names,
                                         interaction.names = interaction_names,
                                         lambda.beta = lam ,
                                         lambda.gamma = 50,
                                         threshold = 1e-4 , max.iter = 200 , initialization.type = "ridge",
                                         nlambda.gamma = 5, nlambda.beta = 10, cores = 1))
res2$dfbeta
res2$dfalpha



dev.off()
par(mfrow=c(2,1))
res2$beta %>% t %>% matplot(type="l")
res2$alpha %>% t %>% matplot(type="l")

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




