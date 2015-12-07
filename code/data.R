## ---- data ----

set.seed(123456)
p = 10 # number of predictors
n = 200 # number of test subjects
m = 200 # number of validation subjects
rho = 0.5
signal_to_noise_ratio = 4

main_effect_names <- paste0("x",1:p) # names of the main effects
# names of the active set
true_var_names <- c("x1","x2","x3","x4","x1:x2", "x1:x3", "x1:x4", "x2:x3", "x2:x4", "x3:x4")
beta1 <- c(7,2,1,1,0,0,0,0,0,0) %>% magrittr::set_names(true_var_names)
beta2 <- c(7,2,1,1,1,0,0,0.5,0.4,0.1) %>% magrittr::set_names(true_var_names)
beta3 <- c(7,2,1,1,7,7,7,2,2,1) %>% magrittr::set_names(true_var_names)
beta4 <- c(7,2,1,1,14,14,14,4,4,2) %>% magrittr::set_names(true_var_names)
beta5 <- c(0,0,0,0,7,7,7,2,2,1) %>% magrittr::set_names(true_var_names)


# simulate Toeplitz like structure
H <- abs(outer(1:p, 1:p, "-"))
cor <- rho^H

DT <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = cor) %>% 
    magrittr::set_colnames(paste0("x",1:p)) %>% 
    set_rownames(paste0("Subject",1:n))
head(DT)

# create X matrix which contains all main effects and interactions
# each column is standardized to mean 0 and sd 1
X <- model.matrix(as.formula(paste0("~(",paste0(main_effect_names, collapse = "+"),")^2-1")), DT %>% as.data.frame()) %>% scale
crossprod(X) %>% diag

apply(X, 2, mean) %>%  sum
apply(X, 2, sd) %>% sum

# # alternative, according to Bien et al 2013, they first standardize the main effects, 
# # then take all pairwise products for interactions, then center those products
# # first scale the main effects
# X_temp <- DT %>% scale
# apply(X_temp, 2, mean) %>%  sum
# apply(X_temp, 2, sd) %>% sum
# 
# # then create all pairwise interactions
# X_temp_2 <- model.matrix(as.formula(paste0("~(",paste0(main_effect_names, collapse = "+"),")^2-1")), X %>% as.data.frame())
# apply(X_temp_2[,(p+1):ncol(X_temp_2)], 2, sd) %>% hist
# apply(X_temp_2[,(p+1):ncol(X_temp_2)], 2, mean) %>% hist
# 
# # then center the product terms and combine with standardized main effects
# X <- cbind(X_temp_2[,1:p], scale(X_temp_2[,(p+1):ncol(X_temp_2)], center = T, scale = F))
# dim(X)
# str(X)
# X %>% crossprod() %>%  diag

y.star <- X[,names(beta4)] %*% beta4
error <- rnorm(n)
k <- sqrt(var(y.star)/(signal_to_noise_ratio*var(error)))
Y <- (y.star + k*error) %>% scale(center = TRUE, scale = FALSE)
colnames(Y) <- "Y"

(b0 <- mean(y.star + k*error))

# names of interaction variables assuming interaction terms contain a ":"
interaction_names <- colnames(X) %>% grep(":",., value = T)