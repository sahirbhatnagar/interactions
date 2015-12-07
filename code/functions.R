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
        if (include.intercept) {
            fit <- lm.fit(x = cbind2(rep(1, nrow(x)),x[,i, drop = F]), y = y )
            fit$coefficients[2] 
        } else {
            fit <- lm.fit(x = x[,i, drop = F], y = y)
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


#' Calculate Adaptive Weights based on Ridge Regression
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
        weights[j,] <- abs(1/betas.and.alphas[j,]) 
    }
    
    for (k in interaction.names) {
        
        # get names of main effects corresponding to interaction
        main <- betas.and.alphas[k, , drop = F] %>% 
                rownames %>% stringr::str_split(":") %>% 
                unlist
        
        weights[k,] <- abs(prod(betas.and.alphas[main,])/betas.and.alphas[k,]) 
    }
    
    return(weights)
}


#' Update Weights based on betas and gammas. Currently not being used.
#' 
#' @description uses betas and gammas to update weights. this is used to update
#' the weights at each iteration of the fitting algorithm in the \code{shim}
#' function
#' @param betas.and.gammas q x 1 data.frame or matrix of betas and 
#' gamma estimates. The rownames must be appropriately labelled because 
#' these labels are be used in this function and must match those in the arguments
#' \code{main.effect.names} and \code{interaction.names}
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#' separated by a ':' (e.g. x1:x2)
#' @param include.intercept logical if intercept should be fitted. default is 
#' FALSE. Should be set to TRUE if y is not centered
#' @note Currently this is not being used in the shim function i.e. we are not
#' updating the weights at each iteration
#' @return q x 1 matrix of weights

update_weights <- function(betas.and.gammas, 
                           main.effect.names, 
                           interaction.names,
                           epsilon = 1e-5) {
    
    # create output matrix
    weights <- matrix(nrow = nrow(betas.and.gammas)) %>% 
        magrittr::set_rownames(rownames(betas.and.gammas))
    
    # main effects weights
    for (j in main.effect.names) {
        
        weights[j,] <- if (betas.and.gammas[j,] < epsilon) 1e7  else 
                        abs(1/betas.and.gammas[j,]) 
    }
    
    for (k in interaction.names) {

        weights[k,] <- if (betas.and.gammas[k,]<epsilon) 1e7 else 
                        abs(1/betas.and.gammas[k,]) 
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
#' @param lambda.beta a single tuning parameter for main effects
#' @param lambda.gamma a single tuning parameter for gammas
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



#' Fit the Strong Heredity Interactions Model
#' 
#' @description This is the main workhorse function that fits the 
#' Strong Heredity Interactions Model (SHIM) of Choi et al 2009 (JASA) for 
#' a given pair of tuning paramters
#' @param x design matrix of dimension n x q, where n is the number of 
#' subjects and q is the total number of variables. This must include all main 
#' effects and interactions as well, with column names corresponding to the 
#' names of the variables (e.g. x1, x2, ...) and their 
#' interactions (e.g. x1:x2, x1:x3, ...). All columns should be scaled to have 
#' mean 0 and variance 1
#' @param y response (matrix form) of dimension n x 1 (should be centered)
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#' separated by a ':' (e.g. x1:x2)
#' @param lambda.beta a single tuning parameter for main effects
#' @param lambda.gamma a single tuning parameter for gammas
#' @param threshold this corresponds to delta in step 5 of the algorithm. It
#' is the threhold for the relative difference between likelihoods at successive
#' iterations (i.e. the algortihm will stop once the relative difference is less
#' than threshold)
#' @param max.iter the maximum number of iterations. If algorithm hasn't 
#' converged, try increasing this number.
#' @return A list containing the following 
#' \enumerate{
#'   \item beta p x niter matrix of beta coefficients at each iteration
#'   \item gamma p*(p-1)/2 x niter matrix of gamma coefficients at each iteration
#'   \item Q a matrix with the iteration number and corresponding value
#'   of the likelihood function
#'   \item m the number of iterations
#' }

shim <- function(x, y, main.effect.names, interaction.names, 
                 lambda.beta, lambda.gamma, threshold, max.iter) {
    
    adaptive.weights <- ridge_weights(x = x, y = y, 
                                       main.effect.names = main.effect.names, 
                                       interaction.names = interaction.names)
    
    # initialization
    betas_and_alphas <- uni_fun(variables = colnames(x), x = x, y = y, include.intercept = F)
    
    # this converts the alphas to gammas
    uni_start <- betas_and_alphas %>% 
                 convert(., main.effect.names = main.effect.names, 
                         interaction.names = interaction.names)
    
    # initialize beta_hat_next also because we are updating each beta individually
    #beta_hat_previous <- beta_hat_next <- uni_start[main.effect.names, , drop = F]
    
    beta_hat_previous <- uni_start[main.effect.names, , drop = F]
    gamma_hat_previous <- uni_start[interaction.names, , drop = F]
    
    m = 1 # iteration counter
    delta = 1 # threshold initialization
    
    # store likelihood values at each iteration in a matrix Q
    Q <- matrix(c(seq(0,max.iter), rep(NA,max.iter+1)), nrow = max.iter+1, ncol = 2) %>% 
         magrittr::set_colnames(c("iteration", "Q.theta"))
    
    # matrix to store betas and gammas of every iteration
    betas <- matrix(nrow = length(main.effect.names), ncol = max.iter + 1)
    gammas <- matrix(nrow = length(interaction.names), ncol = max.iter + 1)
    
    betas[,1] <- beta_hat_previous
    gammas[,1] <- gamma_hat_previous
    
    # store the value of the likelihood at the 0th iteration
    Q[1,2] <- Q_theta(x = x, y = y, beta = beta_hat_previous, 
                      gamma = gamma_hat_previous, weights = adaptive.weights, 
                      lambda.beta = lambda.beta, lambda.gamma = lambda.gamma,
                      main.effect.names = main.effect.names, 
                      interaction.names = interaction.names)
    
    print(paste("Iteration: 0, Q(theta):",Q[1,2]))
    
    while (threshold < delta && m < max.iter){
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # update gamma (interaction parameter)
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y_tilde <- y - x[,main.effect.names] %*% beta_hat_previous
        
        
        x_tilde <- xtilde(interaction.names = interaction.names, 
                           data.main.effects = x[,main.effect.names],
                           beta.main.effects = beta_hat_previous)
        
        # update the gammas using glmnet
        # x_tilde only has the interaction columns, therefore, penalty.factor 
        # must also only include the weights for the interaction terms
        fit_gamma_hat_glmnet <- glmnet::glmnet(x = x_tilde, 
                                               y = y_tilde, 
                                               nlambda = 1, 
                                               lambda = lambda.gamma, 
                                               penalty.factor = adaptive.weights[colnames(x_tilde),],
                                               standardize = T, intercept = T)
        
        # get gamma coefficients and remove intercept
        gamma_hat_next <- coef(fit_gamma_hat_glmnet, s = lambda.gamma) %>% 
                          as.matrix %>% 
                          magrittr::extract(-1, ,drop=F)
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # update beta (main effect parameter) step 4 of algortihm in Choi et al
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        beta_hat_next <- beta_hat_previous
        
        for (j in main.effect.names) {
            
            # determine the main effects not in j
            j_prime_not_in_j <- dplyr::setdiff(main.effect.names,j)
            
            y_tilde_2 <- y - x[,j_prime_not_in_j] %*% beta_hat_next[j_prime_not_in_j,] - 
                (xtilde(interaction.names = interaction.names[-grep(j, interaction.names)],
                        data.main.effects = x[,j_prime_not_in_j],
                        beta.main.effects = beta_hat_next[j_prime_not_in_j,,drop=F]) %>% 
                     rowSums() %>% as.matrix(ncol = 1))
            
            # index data.frame to figure out which j < j'
            index <- data.frame(main.effect.names, seq_along(main.effect.names), 
                                stringsAsFactors = F) %>% 
                        magrittr::set_colnames(c("main.effect.names","index"))
            
            # j' less than j
            j.prime.less <- index[which(index[,"index"] < index[which(index$main.effect.names == j),2]),
                                  "main.effect.names"]
            
            # the if conditions in term1 and term2 are to check if there are 
            # any variables greater or less than j            
            term_1 <- if (length(j.prime.less) != 0) { 
                x[,paste(j.prime.less,j,sep = ":")] %*% 
                (gamma_hat_next[paste(j.prime.less,j,sep = ":"),, drop = F] * 
                     beta_hat_next[j.prime.less,,drop = F])
                } else 0
            
            # j' greater than j
            j.prime.greater <- index[which(index[,"index"] > index[which(index$main.effect.names == j),2]),
                                     "main.effect.names"]
            
            term_2 <- if (length(j.prime.greater) != 0) { 
                x[,paste(j,j.prime.greater,sep = ":")] %*% 
                (gamma_hat_next[paste(j, j.prime.greater,sep = ":"),, drop = F] * 
                     beta_hat_next[j.prime.greater,,drop = F]) 
                } else 0
            
            x_tilde_2 <- x[,j, drop = F] + term_1 + term_2
            
            # need to add a column of zeros to the design matrix, because
            # glmnet returns an error if the design matrix only has one 
            # column. also need to give this column a weight of 0 
            beta_hat_next[j,] <- glmnet(x = cbind2(x_tilde_2, rep(0,nrow(x_tilde_2))), 
                                        y = y_tilde_2, 
                                        nlambda = 1, intercept = T,
                                        lambda = lambda.beta,
                                        standardize = T,
                                        penalty.factor = c(adaptive.weights[colnames(x_tilde_2),],0)) %>%  
                                coef(., s = lambda.beta) %>% 
                                as.matrix %>% 
                                magrittr::extract(colnames(x_tilde_2), ,drop=F)
        }
        
        Q[m+1,2] <- Q_theta(x = x, y = y, beta = beta_hat_next, 
                            gamma = gamma_hat_next, weights = adaptive.weights, 
                            lambda.beta = lambda.beta, lambda.gamma = lambda.gamma,
                            main.effect.names = main.effect.names, 
                            interaction.names = interaction.names)
        
        betas[,m+1] <- beta_hat_next
        gammas[,m+1] <- gamma_hat_next
        
        delta <- abs(Q[m,2] - Q[m+1,2])/abs(Q[m,2])
        
        print(paste("Iteration:",m, ", Q(theta):",Q[m+1,2]))
        
        m = m+1
        
        beta_hat_previous <- beta_hat_next
        
    }
    
    return(list(beta = betas, gamma = gammas, Q = Q, m = m))
    
}



#' My soft thresholding function. Currently not being used

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

