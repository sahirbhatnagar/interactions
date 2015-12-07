##################################
# R source code file used to create simulated data
# from Choi et al 2009 JASA
# Created by Sahir, November 5, 2015
# Updated Dec 6th, 2015
# hosted on Github repo 'sahirbhatnagar/interactions'
# NOTE: need to load packages.R before running this file
##################################

source("packages.R")

## ---- data ----

set.seed(12345)

# number of predictors
p = 10 

# number of test subjects
n = 200 

# correlation between X's
rho = 0.5

# signal to noise ratio
signal_to_noise_ratio = 4

# names of the main effects, this will be used in many of the functions
main_effect_names <- paste0("x",1:p) 

# names of the active set
true_var_names <- c("x1","x2","x3","x4","x1:x2", "x1:x3", "x1:x4", "x2:x3", "x2:x4", "x3:x4")

# different true coefficient vectors as in Table 1 of Choi et al. 
beta1 <- c(7,2,1,1,0,0,0,0,0,0) %>% magrittr::set_names(true_var_names)
beta2 <- c(7,2,1,1,1,0,0,0.5,0.4,0.1) %>% magrittr::set_names(true_var_names)
beta3 <- c(7,2,1,1,7,7,7,2,2,1) %>% magrittr::set_names(true_var_names)
beta4 <- c(7,2,1,1,14,14,14,4,4,2) %>% magrittr::set_names(true_var_names)
beta5 <- c(0,0,0,0,7,7,7,2,2,1) %>% magrittr::set_names(true_var_names)

# simulate Toeplitz like correlation structure between X's
H <- abs(outer(1:p, 1:p, "-"))
cor <- rho^H

# generate X's from multivariate normal and label the matrix
DT <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = cor) %>% 
    magrittr::set_colnames(paste0("x",1:p)) %>% 
    set_rownames(paste0("Subject",1:n))

# create X matrix which contains all main effects and interactions
# but not the intercept
# each column is standardized to mean 0 and sd 1
X <- model.matrix(
    as.formula(paste0("~(",paste0(main_effect_names, collapse = "+"),")^2-1")), 
    data = DT %>% as.data.frame()) %>% scale

# check that means of columns are 0 and sd 1
colMeans(X) %>%  sum
apply(X, 2, sd) %>% sum

# generate response with user defined signal to noise ratio and center 
# the response
y.star <- X[,names(beta4)] %*% beta4
error <- rnorm(n)
k <- sqrt(var(y.star)/(signal_to_noise_ratio*var(error))) 
Y <- (y.star + k*error) %>% scale(center = TRUE, scale = FALSE) 
colnames(Y) <- "Y"

# record mean of response before centering
(b0 <- mean(y.star + k*error))

# names of interaction variables assuming interaction terms contain a ":"
# this will be used in many of the functions
interaction_names <- colnames(X) %>% grep(":",., value = T)