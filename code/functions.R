## ---- functions ----

#' Univariate regressions 
#' 
#' @description function used to create initial estimates in fitting 
#' algorithm
#' @param variables character vector of variable names for which you want 
#' the univariate regression estimate. Must be contained in the 
#' column names of x
#' @param x design matrix of dimension n x q, where n is the number of 
#' subjects and q is the total number of variables. This must include all main 
#' effects and interactions as well, with column names corresponding to the 
#' names of the variables (e.g. x1, x2, ...) and their 
#' interactions (e.g. x1:x2, x1:x3, ...). All columns should be scaled to have 
#' mean 0 and variance 1
#' @param y response (matrix form) of dimension n x 1
#' @param include.intercept logical if intercept should be fitted. default is 
#' FALSE. Should be set to TRUE if y is not centered
#' @return OLS coefficients as a q x 1 data.frame
#' @note to stay consistent with the notation of Choi et al., p is defined as
#' the number of main effects. I have introduced q as being the total number
#' of variables (e.g. the number of columns in the design matrix). For their 
#' specific setting (i.e. all pairwise interactions) q = p + p*(p-1)/2

uni_fun <- function(variables, x, y, include.intercept = F) {
    
    res <- plyr::ldply(variables, function(i) {
            # dont need to add intercept because y has been centered
            fit <- if (include.intercept) {
                lm.fit(x = cbind2(rep(1, nrow(x)),x[,i, drop = F]), y = y )
                fit$coefficients[2] 
                } else {
                    lm.fit(x = x[,i, drop = F], y = y )
                    fit$coefficients[1] 
                } 
            }) %>% 
        magrittr::set_rownames(variables) %>% 
        magrittr::set_colnames("univariate_beta") %>% 
        as.matrix
    
    return(res)
}

#' Convert alphas to gammas
#' 
#' @description function that takes a vector of betas (which are the 
#' main effects) and alphas (which are the interaction effects) and converts 
#' the alphas to gammas.
#' note that 
#' \deqn{y = \beta_0 + \beta_1 x_1 + \cdots + \beta_p x_p 
#' + \alpha_{12} x_1 x_2 + \cdots + \alpha_{p-1,p} x_p x_{p-1} }
#' and
#' \deqn{\alpha_{ij} = \gamma_{ij} * \beta_i*\beta_j , i < j}
#' this function is used because the fitting algorithm estimates the gammas,
#' and furthermore, the L1 penalty is placed on the gammas. It is used only
#' in the initialization step.
#' @param betas.and.alphas q x 1 data.frame or matrix of main effects and 
#' interaction estimates. For example the output from the \code{uni_fun} 
#' function. The rownames must be appropriately labelled because these labels 
#' will be used in other functions
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#' separated by a ':' (e.g. x1:x2)
#' @param epsilon threshold to avoid division by a very small beta e.g. if 
#' any of the main effects are less than epsilon, set gamma to zero. This 
#' should not really be an important parameter because this function is only
#' used in the initialization step, where the intial estimates are from OLS or
#' ridge regression and therefor should not be very close to 0
#' @return a labelled q x 1 data.frame of betas and gammas 

convert <- function(betas.and.alphas, main.effect.names, interaction.names, 
                    epsilon = 1e-5) {
    
    betas_and_gammas <- matrix(nrow = nrow(betas.and.alphas)) %>% 
                        magrittr::set_rownames(rownames(betas.and.alphas))
    
    for (k in interaction.names) {
        
        # get names of main effects corresponding to interaction
        main <- betas.and.alphas[k, , drop = F] %>% 
                rownames %>% stringr::str_split(":") %>% 
                unlist
        
        # convert alpha to gamma BUT NEED TO CHECK IF BETAS ARE 0
        betas_and_gammas[k,] <- if (any(abs(betas.and.alphas[main,]) < epsilon )) 0 else 
                                betas.and.alphas[k,]/prod(betas.and.alphas[main,]) 
    }
    
    # add back the main effects which dont need to be transformed
    for (j in main.effect.names) {
        betas_and_gammas[j,] <- betas.and.alphas[j,]
    }
    
    return(betas_and_gammas)
}

#' Convert gammas to alphas
#' 
#' @description function that takes a vector of betas (which are the 
#' main effects) and gammas and converts the alphas to gammas.
#' This function is used to calculate the linear predictor of the likelihood
#' function (the Q function in the fitting algorithm)
#' @param betas.and.gammas q x 1 data.frame or matrix of betas and 
#' gamma estimates. For example the output from the \code{convert} 
#' function. The rownames must be appropriately labelled because these labels 
#' will be used in other functions
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#' separated by a ':' (e.g. x1:x2)
#' @return a labelled q x 1 data.frame of betas and alphas 

convert2 <- function(betas.and.gammas, main.effect.names, interaction.names) {

    # create output matrix
    betas.and.alphas <- matrix(nrow = nrow(betas.and.gammas)) %>% 
                        magrittr::set_rownames(rownames(betas.and.gammas))
    
    for (k in interaction.names) {
        
        # get names of main effects corresponding to interaction
        main <- betas.and.gammas[k, , drop = F] %>% 
                rownames %>% 
                stringr::str_split(":") %>% 
                unlist
        
        # convert alpha to gamma
        betas.and.alphas[k,] <- betas.and.gammas[k,]*prod(betas.and.gammas[main,]) 
    }
    
    # add back the main effects which dont need to be transformed
    for (j in main.effect.names) {
        betas.and.alphas[j,] <- betas.and.gammas[j,]
    }
    
    return(betas.and.alphas)
    
}

#' Calculate working X's 
#' 
#' @description function used to calculate working X's (xtilde) in 
#' step 3 of algorithm
#' @param interaction.names character vector of interaction names. must be 
#' separated by a ':' (e.g. x1:x2)
#' @param data.main.effects data frame or matrix containing the main effects 
#' data
#' @param beta.main.effects data frame or matrix containing the coefficients 
#' of main effects
#' @return matrix of working X's (xtilde) of dimension n x (p*(p-1)/2)
#' @note Only the main effects are used in this step, 
#' however you can provide this function, either the betas.and.alphas or 
#' betas.and.gammas because only the betas (main effect) parameters 
#' are used in the calculation of xtilde
xtilde <- function(interaction.names, data.main.effects, beta.main.effects){
    
    # create output matrix
    xtildas <- matrix(ncol = length(interaction.names), 
                      nrow = nrow(data.main.effects)) %>% 
               magrittr::set_colnames(interaction.names)  
    
    for (k in interaction.names) {
        
        # get names of main effects corresponding to interaction
        main <- k %>% 
                stringr::str_split(":") %>% 
                unlist
        
        # step 3 to calculate x tilda
        xtildas[,k] <- prod(beta.main.effects[main,]) * 
                       data.main.effects[,main[1],drop = F] * 
                       data.main.effects[,main[2],drop = F]
    }
    
    return(xtildas)
}


#' Calculate Adaptive Weights 
#' 
#' @description uses ridge regression from glmnet package to calculate 
#' the weights used in the fitting algorithm. 
#' @param x design matrix of dimension n x q, where n is the number of 
#' subjects and q is the total number of variables. This must include all main 
#' effects and interactions as well, with column names corresponding to the 
#' names of the variables (e.g. x1, x2, ...) and their 
#' interactions (e.g. x1:x2, x1:x3, ...). All columns should be scaled to have 
#' mean 0 and variance 1
#' @param y response (matrix form) of dimension n x 1
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#' separated by a ':' (e.g. x1:x2)
#' @param include.intercept logical if intercept should be fitted. default is 
#' FALSE. Should be set to TRUE if y is not centered
#' @return q x 1 matrix of weights

ridge_weights <- function(x, y, main.effect.names, interaction.names, 
                          include.intercept = F) {

    # fit the ridge to get betas and alphas
    fit <- glmnet::cv.glmnet(x = x, y = y, alpha = 0, 
                              standardize = F, 
                              intercept = include.intercept)
    
    n <- length(y)
    
    # remove intercept (even if include.intercept is FALSE, coef.glmnet returns
    # an intercept set to 0)
    betas.and.alphas <- coef(fit, s = "lambda.1se") %>% 
                        as.matrix() %>% 
                        magrittr::extract(-1, ,drop = F)
    
    # create output matrix
    weights <- matrix(nrow = nrow(betas.and.alphas)) %>% 
               magrittr::set_rownames(rownames(betas.and.alphas))
    
    # main effects weights
    for (j in main.effect.names) {
        weights[j,] <- abs(1/betas.and.alphas[j,]) * log(n)/n
    }
    
    for (k in interaction.names) {
        
        # get names of main effects corresponding to interaction
        main <- betas.and.alphas[k, , drop = F] %>% 
                rownames %>% stringr::str_split(":") %>% 
                unlist
        
        weights[k,] <- abs(prod(betas.and.alphas[main,])/betas.and.alphas[k,]) * log(n)/n 
    }
    
    return(weights)
}

#' Likelihood function
#' 
#' @description calculates likelihood function. Used to assess convergence
#' of fitting algorithm. This corresponds to the Q(theta) function in the 
#' paper
#' @param x design matrix of dimension n x q, where n is the number of 
#' subjects and q is the total number of variables. This must include all main 
#' effects and interactions as well, with column names corresponding to the 
#' names of the variables (e.g. x1, x2, ...) and their 
#' interactions (e.g. x1:x2, x1:x3, ...). All columns should be scaled to have 
#' mean 0 and variance 1
#' @param y response (matrix form) of dimension n x 1
#' @param beta p x 1 matrix of main effect estimates
#' @param gamma p*(p-1)/2 x 1 matrix of gamma estimates
#' @param weights adaptive weights calculated by \code{ridge_weights} function
#' with rownames corresponding to column names of x
#' @param lambda.beta tuning parameter for main effects
#' @param lambda.gamma tuning parameter for gammas
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#' separated by a ':' (e.g. x1:x2)
#' @return value of likelihood function

Q_theta <- function(x, y, beta, gamma, weights, 
                    lambda.beta, lambda.gamma, main.effect.names, 
                    interaction.names){
    
    # first convert gammas to alphas which will be used to calculate
    # the linear predictor
    betas.and.alphas <- convert2(betas.and.gammas = rbind2(beta,gamma), 
                                 main.effect.names = main.effect.names, 
                                 interaction.names = interaction.names)
    
    crossprod(y - x %*% betas.and.alphas) +  
        lambda.beta * (crossprod(weights[main.effect.names,], abs(beta))) + 
        lambda.gamma * (crossprod(weights[interaction.names,], abs(gamma)))
}



soft <- function(x, y, beta, lambda, weight) {
    # user must supply x AND y, or beta.. but not both
    # i set it up this way because to get the sequence of lambdas, I use the beta argument
    # so that I only compute this once. I use the x, y argument for the CV folds
    # lambda can be a vector and this functions will return each thresholded beta for each 
    # lambda e.g. soft(0.25, lambda = seq(0.001,0.65,length.out = 100), 1.5)
    
    if (missing(x) & missing(y) & missing(beta)) stop("user must supply x AND y, or beta but not both")
    if (missing(x) & missing(y)) return(list("beta" = sign(beta) * pmax(0, abs(beta) - lambda*weight)))
    if (missing(beta)) {
        #(beta <- lm.fit(x = cbind2(rep(1, length(y)),x), y = y) %>% coef %>% magrittr::extract(2))
        (beta <- lm.fit(x = x[,1,drop=F], y = y) %>% coef %>% magrittr::extract(1))
        
        #lm.fit(x = cbind2(rep(1, length(y_tilde_2)),x_tilde_2), y = y_tilde_2) %>% coef %>% magrittr::extract(2)
        b_lasso <- sign(beta)* pmax(0, abs(beta) - lambda*weight)
        #return(list("beta" = b_lasso, "df" = nonzero(b_lasso)))
        return(b_lasso)
    }
    
}

