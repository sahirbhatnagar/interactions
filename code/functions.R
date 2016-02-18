##################################
# R source code file for functions
# Created by Sahir, Dec 7, 2015
# Updated 
# hosted on Github repo 'sahirbhatnagar/interactions'
# NOTE: the main function is shim. Most of the other functions
# are being called by shim
##################################


## ---- functions ----

#' Univariate regressions
#' 
#' @description function used to create initial estimates in fitting algorithm
#' @param variables character vector of variable names for which you want the
#'   univariate regression estimate. Must be contained in the column names of x
#' @param x design matrix of dimension n x q, where n is the number of subjects
#'   and q is the total number of variables. This must include all main effects
#'   and interactions as well, with column names corresponding to the names of
#'   the variables (e.g. x1, x2, ...) and their interactions (e.g. x1:x2, x1:x3,
#'   ...). All columns should be scaled to have mean 0 and variance 1
#' @param y response (matrix form) of dimension n x 1
#' @param include.intercept logical if intercept should be fitted. default is 
#'   FALSE. Should be set to TRUE if y is not centered
#' @param The procedure used to estimate betas and alphas. must be the character
#'   "ridge" or "univariate".
#' @return OLS coefficients as a q x 1 data.frame
#' @note to stay consistent with the notation of Choi et al., p is defined as 
#'   the number of main effects. I have introduced q as being the total number 
#'   of variables (e.g. the number of columns in the design matrix). For their 
#'   specific setting (i.e. all pairwise interactions) q = p + p*(p-1)/2

uni_fun <- function(variables, x, y, include.intercept = F, type="ridge") {
  
  res <- switch(type,
                univariate = {plyr::ldply(variables, function(i) {
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
                    as.matrix},
                ridge = {
                  # fit the ridge to get betas and alphas
                  glmnet::cv.glmnet(x = x, y = y, alpha = 0, 
                                    standardize = F, 
                                    intercept = include.intercept) %>% 
                    # remove intercept (even if include.intercept is FALSE, coef.glmnet returns
                    # an intercept set to 0)
                    coef(., s = "lambda.min") %>% 
                    as.matrix() %>% 
                    magrittr::extract(-1, ,drop = F)
                })
  
  return(res)
}

#' Convert alphas to gammas
#' 
#' @description function that takes a vector of betas (which are the main
#'   effects) and alphas (which are the interaction effects) and converts the
#'   alphas to gammas. note that \deqn{y = \beta_0 + \beta_1 x_1 + \cdots +
#'   \beta_p x_p + \alpha_{12} x_1 x_2 + \cdots + \alpha_{p-1,p} x_p x_{p-1} } 
#'   and \deqn{\alpha_{ij} = \gamma_{ij} * \beta_i*\beta_j , i < j} this
#'   function is used because the fitting algorithm estimates the gammas, and
#'   furthermore, the L1 penalty is placed on the gammas. It is used only in the
#'   initialization step.
#' @param betas.and.alphas q x 1 data.frame or matrix of main effects and 
#'   interaction estimates. For example the output from the \code{uni_fun} 
#'   function. The rownames must be appropriately labelled because these labels 
#'   will be used in other functions
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#'   separated by a ':' (e.g. x1:x2)
#' @param epsilon threshold to avoid division by a very small beta e.g. if any
#'   of the main effects are less than epsilon, set gamma to zero. This should
#'   not really be an important parameter because this function is only used in
#'   the initialization step, where the intial estimates are from OLS or ridge
#'   regression and therefor should not be very close to 0
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
#' @description function that takes a vector of betas (which are the main
#'   effects) and gammas and converts the alphas to gammas. This function is
#'   used to calculate the linear predictor of the likelihood function (the Q
#'   function in the fitting algorithm)
#' @param betas.and.gammas q x 1 data.frame or matrix of betas and gamma
#'   estimates. For example the output from the \code{convert} function. The
#'   rownames must be appropriately labelled because these labels will be used
#'   in other functions
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#'   separated by a ':' (e.g. x1:x2)
#' @return a labelled q x 1 data.frame of betas and alphas

convert2 <- function(beta, gamma, main.effect.names, interaction.names,
                     intercept = NULL) {
  
    betas.and.gammas <- rbind2(beta,gamma)
    #rownames(betas.and.gammas) <- c(main.effect.names, interaction.names)
    # create output matrix
    betas.and.alphas <- matrix(nrow = nrow(betas.and.gammas)) %>% 
        magrittr::set_rownames(rownames(betas.and.gammas))
    
    for (k in interaction.names) {
        
        # get names of main effects corresponding to interaction
        main <- betas.and.gammas[k, , drop = F] %>% 
            rownames %>% 
            stringr::str_split(":") %>% 
            unlist
        
        # convert gamma to alpha
        betas.and.alphas[k,] <- betas.and.gammas[k,]*prod(betas.and.gammas[main,]) 
    }
    
    # add back the main effects which dont need to be transformed
    for (j in main.effect.names) {
        betas.and.alphas[j,] <- betas.and.gammas[j,]
    }
    
    # add back intercept if it is non-NULL
    if (!is.null(intercept)) betas.and.alphas["(Intercept)",] <- intercept
    
    return(betas.and.alphas)
    
}



#' Calculate working X's to update Gammas.
#' 
#' @description function used to calculate working X's (xtilde) in step 3 of
#'   algorithm
#' @param interaction.names character vector of interaction names. must be 
#'   separated by a ':' (e.g. x1:x2)
#' @param data.main.effects data frame or matrix containing the main effects 
#'   data
#' @param beta.main.effects data frame or matrix containing the coefficients of
#'   main effects
#' @param nlambda number of tuning parameters
#' @return matrix of working X's (xtilde) of dimension n x (p*(p-1)/2)
xtilde <- function(interaction.names, data.main.effects, beta.main.effects){
  
  # create output matrix
  xtildas <- matrix(ncol = length(interaction.names), 
                    nrow = nrow(data.main.effects))
  colnames(xtildas) <- interaction.names
  
  for (k in interaction.names) {
    
    # get names of main effects corresponding to interaction
    main <- unlist(stringr::str_split(k,":"))
    
    # step 3 to calculate x tilda
    xtildas[,k] <- prod(beta.main.effects[main,]) * 
      data.main.effects[,main[1],drop = F] * 
      data.main.effects[,main[2],drop = F]
  }
  
  return(xtildas)
}



#' Calculate working X's to update Betas
#' 
#' @description function used to calculate working X's (xtilde) in step 4 of
#'   algorithm
#' @param interaction.names character vector of interaction names. must be 
#'   separated by a ':' (e.g. x1:x2)
#' @param data.main.effects data frame or matrix containing the main effects 
#'   data
#' @param beta.main.effects data frame or matrix containing the coefficients of
#'   main effects
#' @param gamma.interaction.effects data frame or matrix containing the gamma
#'   parameters
#' @return matrix of working X's (xtilde) of dimension n x (p*(p-1)/2)
#' @note this function is a modified x_tilde for step 4 because we thought maybe
#'   there was a typo. Math and results suggests that there IS a typo. This is
#'   now being used
xtilde_mod <- function(interaction.names, data.main.effects, beta.main.effects, 
                       gamma.interaction.effects){
  
  # create output matrix. no pipe is faster
  xtildas <- matrix(ncol = length(interaction.names), 
                    nrow = nrow(data.main.effects)) 
  colnames(xtildas) <- interaction.names
  
  for (k in interaction.names) {
    
    # get names of main effects corresponding to interaction
    main <- unlist(stringr::str_split(k, ":"))
    
    # step 4 to calculate x tilda
    xtildas[,k] <- prod(beta.main.effects[main,]) * gamma.interaction.effects[k,] *  
      data.main.effects[,main[1],drop = F] * 
      data.main.effects[,main[2],drop = F]
  }
  
  return(xtildas)
}




#' Calculate Adaptive Weights based on Ridge Regression
#' 
#' @description uses ridge regression from glmnet package to calculate the
#'   weights used in the fitting algorithm.
#' @param x design matrix of dimension n x q, where n is the number of subjects
#'   and q is the total number of variables. This must include all main effects
#'   and interactions as well, with column names corresponding to the names of
#'   the variables (e.g. x1, x2, ...) and their interactions (e.g. x1:x2, x1:x3,
#'   ...). All columns should be scaled to have mean 0 and variance 1
#' @param y response (matrix form) of dimension n x 1
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#'   separated by a ':' (e.g. x1:x2)
#' @param include.intercept logical if intercept should be fitted. default is 
#'   FALSE. Should be set to TRUE if y is not centered
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


#' Likelihood function
#' 
#' @description calculates likelihood function. Used to assess convergence of
#'   fitting algorithm. This corresponds to the Q(theta) function in the paper
#' @param x design matrix of dimension n x q, where n is the number of subjects
#'   and q is the total number of variables. This must include all main effects
#'   and interactions as well, with column names corresponding to the names of
#'   the variables (e.g. x1, x2, ...) and their interactions (e.g. x1:x2, x1:x3,
#'   ...). All columns should be scaled to have mean 0 and variance 1
#' @param y response (matrix form) of dimension n x 1
#' @param beta p x 1 matrix of main effect estimates
#' @param gamma p*(p-1)/2 x 1 matrix of gamma estimates
#' @param weights adaptive weights calculated by \code{ridge_weights} function 
#'   with rownames corresponding to column names of x
#' @param lambda.beta a single tuning parameter for main effects
#' @param lambda.gamma a single tuning parameter for gammas
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#'   separated by a ':' (e.g. x1:x2)
#' @return value of likelihood function
#' @note you dont use the intercept in the calculation of the Q function
#' because its not being penalized

Q_theta <- function(x, y, beta, gamma, weights, 
                    lambda.beta, lambda.gamma, main.effect.names, 
                    interaction.names){
  
    # first convert gammas to alphas which will be used to calculate
    # the linear predictor
    betas.and.alphas <- convert2(beta = beta,
                                 gamma = gamma, 
                                 main.effect.names = main.effect.names, 
                                 interaction.names = interaction.names)
    
    crossprod(y - x %*% betas.and.alphas) +  
        lambda.beta * (crossprod(weights[main.effect.names,], abs(beta))) + 
        lambda.gamma * (crossprod(weights[interaction.names,], abs(gamma)))
}

#' Fit the Strong Heredity Interactions Model
#' 
#' @description This is the main workhorse function that fits the Strong
#'   Heredity Interactions Model (SHIM) of Choi et al 2009 (JASA) for a given
#'   pair of tuning paramters
#' @param x design matrix of dimension n x q, where n is the number of subjects
#'   and q is the total number of variables. This must include all main effects
#'   and interactions as well, with column names corresponding to the names of
#'   the variables (e.g. x1, x2, ...) and their interactions (e.g. x1:x2, x1:x3,
#'   ...). All columns should be scaled to have mean 0 and variance 1
#' @param y response (matrix form) of dimension n x 1 (should be centered)
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#'   separated by a ':' (e.g. x1:x2)
#' @param lambda.beta a single tuning parameter for main effects
#' @param lambda.gamma a single tuning parameter for gammas
#' @param threshold this corresponds to delta in step 5 of the algorithm. It is
#'   the threhold for the relative difference between likelihoods at successive 
#'   iterations (i.e. the algortihm will stop once the relative difference is
#'   less than threshold)
#' @param max.iter the maximum number of iterations. If algorithm hasn't 
#'   converged, try increasing this number.
#' @param initialization.type The procedure used to initialize betas and gammas.
#'   must be the character "ridge" or "univariate". This argument is passed to
#'   the \code{uni_fun} function
#' @return A list containing the following \enumerate{ \item beta p x niter
#'   matrix of beta coefficients at each iteration \item gamma p*(p-1)/2 x niter
#'   matrix of gamma coefficients at each iteration \item Q a matrix with the
#'   iteration number and corresponding value of the likelihood function \item m
#'   the number of iterations }

shim <- function(x, y, main.effect.names, interaction.names, 
                 lambda.beta, lambda.gamma, threshold, max.iter, 
                 initialization.type = "ridge") {
  
  adaptive.weights <- ridge_weights(x = x, y = y, 
                                    main.effect.names = main.effect.names, 
                                    interaction.names = interaction.names)
  
  # initialization
  betas_and_alphas <- uni_fun(variables = colnames(x), x = x, y = y, 
                              include.intercept = F,
                              type = initialization.type)
  
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
                                           standardize = F, intercept = F)
    
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
      
      #             y_tilde_2 <- y - x[,j_prime_not_in_j] %*% beta_hat_next[j_prime_not_in_j,] - 
      #                 (xtilde(interaction.names = interaction.names[-grep(j, interaction.names)],
      #                         data.main.effects = x[,j_prime_not_in_j],
      #                         beta.main.effects = beta_hat_next[j_prime_not_in_j,,drop=F]) %>% 
      #                      rowSums() %>% as.matrix(ncol = 1))
      
      y_tilde_2 <- y - x[,j_prime_not_in_j] %*% beta_hat_next[j_prime_not_in_j,] - 
        (xtilde_mod(interaction.names = interaction.names[-grep(j, interaction.names)],
                    data.main.effects = x[,j_prime_not_in_j],
                    beta.main.effects = beta_hat_next[j_prime_not_in_j,,drop=F],
                    gamma.interaction.effects = gamma_hat_next) %>% 
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
                                  nlambda = 1, intercept = F,
                                  lambda = lambda.beta,
                                  standardize = F,
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

#' Repeat each column of a matrix, n times
repcol <- function(x, n) {
  s = NCOL(x)
  matrix(x[, rep(1:s, each = n)], nrow = NROW(x), ncol = NCOL(x) * n)
}


#' Calculate Standard Deviation with factor N
#' 
#' @description this is to calculate standard deviation but with divisor of n
#'   and not n-1
#' @param i a vector of numerics

mysd <- function(i) sqrt(crossprod(i - mean(i))/length(i))


#' Calculate Sequence of Tuning Parameters
#' 
#' @description function to calculate the sequence of tuning parameters. This 
#'   formula is taken from section 2.5 of the glmnet paper in Journal of Stat. 
#'   Software
#'   
#' @param x matrix with rows as subjects and columns as variables
#' @param y a 1 column matrix
#' @param lambda.factor The factor for getting the minimal lambda in lambda 
#'   sequence, where min(lambda) = lambda.factor * max(lambda). max(lambda) is 
#'   the smallest value of lambda for which all coefficients are zero. The 
#'   default depends on the relationship between N (the number of rows in the 
#'   matrix of predictors) and p (the number of predictors). If N > p, the 
#'   default is 0.0001, close to zero. If N<p, the default is 0.01. A very small
#'   value of lambda.factor will lead to a saturated fit.
#' @param nlambda the number of lambda values - default is 100.
#' @param scale_x should the columns of x be scaled - default is FALSE
#' @param center_y should y be mean centered - default is FALSE
#' @note The maximum lambda is calculated using the following inequality: 
#'   \deqn{\frac{1}{N*w_j}\abs{\sum_{i=1}^{n} x_{ij}y_i  } \leq \lambda_{max}}
#' @note The minimum lambda is given by lambda.factor*lambda_max. The sequence
#'   of nlambda values are decreasing from lambda_max to lambda_min on the log
#'   scale
#' @note the penalty factors are internally rescaled to sum to the number of 
#'   predictor variables in glmnet. Therefore, to get the correct sequence of 
#'   lambdas when there are weights, this function first rescales the weights 
#'   and then calclated the sequence of lambdas
#' @example \dontrun{ lambda_sequence(X,Y, weights = adaptive.weights)}

lambda_sequence <- function(x, y, weights = NULL,
                            lambda.factor = ifelse(nobs < nvars, 0.01, 1e-06),
                            nlambda = 100, scale_x = F, center_y = F) {
  
    # x = X; y = Y ; weights = adaptive.weights
    # lambda.factor =  1e-07
    # nlambda = 100; scale_x = F; center_y = F
    # this shows that the lambdas are different is we are using weights!!
    # glmnet(x,y,standardize = F, intercept = F)$lambda
    # glmnet(x,y,standardize = F, intercept = F, penalty.factor = adaptive.weights)$lambda
  
  
  # when scaling, first you center then you standardize
  if (any(as.vector(weights) < 0)) stop("Weights must be positive")
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  
  if (length(as.vector(weights)) < nvars ) stop("You must provide weights for 
                                                every column of x")
  
  # scale the weights to sum to nvars
  w <- if (is.null(weights)) rep(1, nvars) else as.vector(weights)/sum(as.vector(weights))*nvars
  
  sx <- if (scale_x) apply(x,2, function(i) scale(i, center = TRUE, scale = mysd(i))) else x
  sy <- if (center_y) as.vector(scale(y, center = T, scale = F)) else as.vector(y)
  lambda.max <- max(abs(colSums(sy * sx)/w))/nrow(sx)
  
  rev(exp(seq(log(lambda.factor*lambda.max), log(lambda.max), length.out = nlambda)))
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
    #(beta <- lm.fit(x = cbind2(rep(1, length(y)),scale(x)), y = y) %>% coef %>% magrittr::extract(2))
    (beta <- lm.fit(x = x[, 1, drop = F], y = y) %>% coef %>% magrittr::extract(1))
    
    #lm.fit(x = cbind2(rep(1, length(y_tilde_2)),x_tilde_2), y = y_tilde_2) %>% coef %>% magrittr::extract(2)
    b_lasso <- sign(beta)* pmax(0, abs(beta) - lambda*weight)
    #return(list("beta" = b_lasso, "df" = nonzero(b_lasso)))
    # need to return a matrix, because this is used in the step to calculate y_tilde in the 
    # shim_multiple function
    return(matrix(b_lasso, ncol = 1))
  }
  
}


check_col_0 <- function(M) { M[, colSums(abs(M)) != 0, drop = F] }


#' Fit Strong Heredity model with one iteration just to get the 
#' sequence of lambda_gamma and lambda_beta 
shim_once <- function(x, y, main.effect.names, interaction.names, 
                      initialization.type = "ridge", 
                      nlambda.gamma = 20, 
                      nlambda.beta = 20,
                      lambda.factor = ifelse(nobs<nvars,0.01,1e-6)) {
  
    # x = X; y = Y; main.effect.names = main_effect_names;
    # interaction.names = interaction_names;
    # lambda.beta = NULL ; lambda.gamma = NULL
    # threshold = 1e-5 ; max.iter = 500 ; initialization.type = "ridge";
    # nlambda.gamma = 10; nlambda.beta = 20; cores = 2
  
  np = dim(x)
  nobs = np[1]
  nvars = np[2]
  
  # total number of tuning parameters
  nlambda = nlambda.gamma * nlambda.beta
  
  adaptive.weights <- ridge_weights(x = x, y = y, 
                                    main.effect.names = main.effect.names, 
                                    interaction.names = interaction.names)
  
  # initialization
  betas_and_alphas <- uni_fun(variables = colnames(x), x = x, y = y, 
                              include.intercept = F,
                              type = initialization.type)
  
  # this converts the alphas to gammas
  uni_start <- convert(betas_and_alphas, main.effect.names = main.effect.names, 
                       interaction.names = interaction.names)
  
  beta_hat_previous <- uni_start[main.effect.names, , drop = F]
  gamma_hat_previous <- uni_start[interaction.names, , drop = F]
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # get tuning parameters for gamma
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # this is a nsubjects x lambda matrix for each tuning parameter stored in a list
  # each element of the list corresponds to a tuning parameter      
  y_tilde <- y - x[,main.effect.names,drop = F] %*% beta_hat_previous
  
  # calculate x_tilde for each beta vector corresponding to a diffent tuning parameter
  x_tilde <- xtilde(interaction.names = interaction.names, 
                    data.main.effects = x[,main.effect.names, drop = F],
                    beta.main.effects = beta_hat_previous) 
  
  # get the sequence of lambda_gammas using the first iteration of
  # x_tilde and y_tilde
  # x_tilde only has the interaction columns, therefore, penalty.factor 
  # must also only include the weights for the interaction terms
  
#   lambda_gamma <- lambda_sequence(x = x_tilde, 
#                                   y = y_tilde,
#                                   weights = adaptive.weights[colnames(x_tilde), ],
#                                   nlambda = nlambda.gamma,
#                                   scale_x = F, center_y = F)
  
  # check that this lambda sequence matches that of glmnet
  fit_gamma_hat_glmnet <- glmnet::glmnet(x = x_tilde, 
                                         y = y_tilde, 
                                         lambda = NULL,
                                         nlambda = nlambda.gamma,
                                         penalty.factor = adaptive.weights[colnames(x_tilde), ],
                                         standardize = F, intercept = F,
                                         lambda.min.ratio = lambda.factor)
  #   # record sequence of lambda_gammas. this will be used in subsequent iterations
  #   # the lambdas generated by glmnet go from large to small
  #   plot(lambda_gamma, type = "l")  
  lambda_gamma_glmnet <- fit_gamma_hat_glmnet$lambda
  
  # get gamma coefficients and remove intercept
  # this results in a matrix of size p*(p-1)/2 x nlambda_gamma i.e. 
  # the number of interaction variables by the number of lambda_gammas
  gamma_hat_next <- as.matrix(coef(fit_gamma_hat_glmnet))[-1, , drop = F]
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # get tuning parameters for beta
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  beta_hat_next <- beta_hat_previous
  
  # for the lambda_beta sequence, calculate the sequences for each 
  # beta, and then take the sequence that contains the maximum 
  # this will ensure that all betas will be 0 under the maximum lambda
  
  # this is to store the lambda_beta sequences for each main.effect.names
  # then we will determine the maximum and minimum lambda_beta across
  # all main effects, for each of the nlambda.beta*nlambda.gamma combinations
  lambda_beta_temp <- vector("list", length(main.effect.names))
  names(lambda_beta_temp) <- main.effect.names
  
  
  # for each main effect, and for each nlambda_gamma,
  # we get a sequence of nlambda_beta tuning parameters
  # only store the max and min for each main effect for each lambda_gamma
  lambda_beta_seq_for_every_lambda_gamma <- replicate(nlambda.gamma,
                                                      matrix(nrow = 2, 
                                                             ncol = length(main.effect.names),
                                                             dimnames = list(paste0("lambda_",c("min","max")),
                                                                             main.effect.names)), 
                                                      simplify = "array")
  
  for (k in seq_len(ncol(gamma_hat_next))) {
    
    for (j in main.effect.names) {
      
      # k=1
      # j = "x10"
      #print(paste(j))
      # determine the main effects not in j
      j_prime_not_in_j <- dplyr::setdiff(main.effect.names,j)
      
      y_tilde_2 <- y - 
        x[,j_prime_not_in_j, drop = F] %*% beta_hat_next[j_prime_not_in_j, , drop = F] - 
        as.matrix(rowSums(xtilde_mod(beta.main.effects = beta_hat_next[j_prime_not_in_j, , drop = F],
                                     gamma.interaction.effects = gamma_hat_next[,k,drop = F],
                                     interaction.names = interaction.names[-grep(j, interaction.names)], 
                                     data.main.effects = x[,j_prime_not_in_j, drop = F])), ncol = 1)
      
      # index data.frame to figure out which j < j'
      index <- data.frame(main.effect.names, seq_along(main.effect.names), 
                          stringsAsFactors = F) 
      colnames(index) <- c("main.effect.names","index")
      
      # j' less than j
      j.prime.less <- index[which(index[,"index"] < index[which(index$main.effect.names == j),2]),
                            "main.effect.names"]
      
      # the if conditions in term1 and term2 are to check if there are 
      # any variables greater or less than j            
      # lapply is faster than mclapply here
      term_1 <- if (length(j.prime.less) != 0) { 
        x[,paste(j.prime.less,j,sep = ":")] %*% 
          (gamma_hat_next[paste(j.prime.less,j,sep = ":"),k, drop = F] * 
             beta_hat_next[j.prime.less,,drop = F])} else matrix(rep(0,nrow(x)), ncol = 1)
      
      # j' greater than j
      j.prime.greater <- index[which(index[,"index"] > index[which(index$main.effect.names == j),2]),
                               "main.effect.names"]
      
      term_2 <- if (length(j.prime.greater) != 0) { 
        x[,paste(j,j.prime.greater,sep = ":")] %*% 
          (gamma_hat_next[paste(j, j.prime.greater,sep = ":"),k, drop = F] * 
             beta_hat_next[j.prime.greater,,drop = F]) 
      } else matrix(rep(0,nrow(x)), ncol = 1)
      
      x_tilde_2 <- x[,j, drop = F] + term_1 + term_2
      
      lambda_beta_seq <- lambda_sequence(x_tilde_2,
                                         y_tilde_2, 
                                         weights = adaptive.weights[colnames(x_tilde_2),],
                                         nlambda = nlambda.beta)
      
      # the seqeunce of lambda_betas for variable j 
      lambda_beta_seq_for_every_lambda_gamma[ , j ,k] <- c(min(lambda_beta_seq), max(lambda_beta_seq))
    }
  }
  
  return(list(lambda_gamma = lambda_gamma_glmnet, 
              lambda_beta = lapply(seq_len(nlambda.gamma), function(i) 
                rev(exp(seq(log(min(lambda_beta_seq_for_every_lambda_gamma[ , ,i])), 
                            log(max(lambda_beta_seq_for_every_lambda_gamma[ , ,i])), 
                            length.out = nlambda.beta))))))
}

#' Fit Strong Heredity Model for Multiple Lambdas
#' @param nlambda.gamma number of tuning parameters for gamma
#' @param nlambda.beta number of tuning parameters for beta
#' @param cores number of cores to use. this is used in the step to calculate
#'   
#' @note let glmnet choose the lambda_betas and lambda_gammas
#' @note the index of the tuning parameters is as follows. If for example there
#'   are 10 lambda_gammas, and 20 lambda_betas, then the first lambda_gamma gets
#'   repeated 20 times. So the first twenty entries of tuning parameters 
#'   correspond to 1 lambda_gamma and the 20 lambda_betas
shim_multiple <- function(x, y, main.effect.names, interaction.names, 
                          lambda.beta = NULL, lambda.gamma = NULL, threshold, max.iter, 
                          initialization.type = "ridge", 
                          nlambda.gamma = 20, 
                          nlambda.beta = 20,
                          cores = 2) {
  

#   x = X; y = Y; main.effect.names = main_effect_names;
#   interaction.names = interaction_names;
#   lambda.beta = NULL ; lambda.gamma = NULL
#   threshold = 1e-5 ; max.iter = 500 ; initialization.type = "ridge";
#   nlambda.gamma = 5; nlambda.beta = 10; cores = 2

  if (is.null(lambda.gamma) & is.null(lambda.beta)) {
    
    (tuning_params <- shim_once(x = x, y = y, 
                               main.effect.names = main.effect.names,
                               interaction.names = interaction.names,
                               initialization.type = "ridge",
                               nlambda.gamma = nlambda.gamma, nlambda.beta = nlambda.beta))
    
    # convert to a list. each element corresponds to a value of lambda_gamma
    lambda_gamma_list <- rep(lapply(seq_len(length(tuning_params$lambda_gamma)), 
                                    function(i) tuning_params$lambda_gamma[i]),
                             each = nlambda.beta)
    
    lambda_beta_list <- lapply(seq_len(length(unlist(tuning_params$lambda_beta))), 
                               function(i) unlist(tuning_params$lambda_beta)[i])
  } else {
    # convert to a list. each element corresponds to a value of lambda_gamma
    lambda_gamma_list <- rep(lapply(seq_len(length(lambda.gamma)), 
                                    function(i) lambda.gamma[i]),
                             each = nlambda.beta)
    
    lambda_beta_list <- lapply(seq_len(length(unlist(lambda.beta))), 
                               function(i) unlist(lambda.beta)[i])
    
  }
  
  # total number of tuning parameters
  nlambda = nlambda.gamma * nlambda.beta
  
  adaptive.weights <- ridge_weights(x = x, y = y, 
                                    main.effect.names = main.effect.names, 
                                    interaction.names = interaction.names)
  
  # initialization
  betas_and_alphas <- uni_fun(variables = colnames(x), x = x, y = y, 
                              include.intercept = F,
                              type = initialization.type)
  
  # this converts the alphas to gammas
  uni_start <- convert(betas_and_alphas, main.effect.names = main.effect.names, 
                       interaction.names = interaction.names)
  
  # need to create a matrix here instead of a 1 column vector  
  # dim1: # of variables, 
  # dim2: # of lambdas 
  
  beta_hat_previous <- replicate(nlambda, uni_start[main.effect.names, , drop = F], 
                                 simplify = "matrix")
  rownames(beta_hat_previous) <- main_effect_names
  
  gamma_hat_previous <- replicate(nlambda, uni_start[interaction.names, , drop = F], 
                                  simplify = "matrix")
  rownames(gamma_hat_previous) <- interaction_names
  
  # convert gamma and beta previous to lists each element corresponds to the
  # coefficients for each combination of lambda_gamma and lambda_beta
  beta_hat_previous_list <- lapply(seq_len(ncol(beta_hat_previous)), 
                                   function(i) beta_hat_previous[,i, drop = F])
  gamma_hat_previous_list <- lapply(seq_len(ncol(gamma_hat_previous)), 
                                    function(i) gamma_hat_previous[,i, drop = F])
  
  m = 1 # iteration counter
  delta = 1 # threshold initialization
  #converged = rep(1, nlambda) # vector for convergence of each tuning parameter
  
  # store likelihood values at each iteration in a matrix Q
  # piping using magrittr::set_colnames is slower here
  # rows are the iterations, columns are the index of the sequence of 
  # lambda_gammas and lambda_betas
  # rows: iteration number
  # columns: tuning parameter
  Q <- matrix(nrow = max.iter+1, ncol = nlambda)
  
  # store the value of the likelihood at the 0th iteration
  Q[1,] <- parallel::mcmapply(Q_theta,
                              beta = beta_hat_previous_list, 
                              gamma = gamma_hat_previous_list, 
                              lambda.beta = lambda_beta_list, 
                              lambda.gamma = lambda_gamma_list,
                              MoreArgs = list(x = x, y = y, 
                                              weights = adaptive.weights,
                                              main.effect.names = main.effect.names, 
                                              interaction.names = interaction.names),
                              mc.cores = cores)
  
  while (threshold < delta && m < max.iter){
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # update gamma (interaction parameter)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # this is a nsubjects x lambda matrix for each tuning parameter stored in a list
    # each element of the list corresponds to a tuning parameter      
    y_tilde_list <- lapply(beta_hat_previous_list, 
                           function(i) y - x[,main.effect.names,drop = F] %*% i)
    
    # calculate x_tilde for each beta vector corresponding to a diffent tuning parameter
    x_tilde_list <- lapply(beta_hat_previous_list, 
                           function(i) xtilde(interaction.names = interaction.names, 
                                              data.main.effects = x[,main.effect.names, drop = F],
                                              beta.main.effects = i)) 

    # adaptive weight for each tuning parameter. currently this is the
    # same for iterations, but I am coding it here
    # for flexibility in case we want to change the weights at each iteration
    adaptive_weights_list <- lapply(x_tilde_list, function(i)
      adaptive.weights[colnames(i),,drop = F])
    
    # indices of the x_tilde matrices that have all 0 columns
    zero_x_tilde <- which(sapply(x_tilde_list, 
                                 function(i) is.null(colnames(check_col_0(i)))))
    
    # for all the x_tilde in zero_x_tilde, return the following matrix with 0 for each coefficient
    # this is like a place holder.
    coef_zero_gamma_matrix <- matrix(data = 0, 
                                     nrow = length(interaction.names), 
                                     ncol = 1, 
                                     dimnames = list(interaction.names))
    
    gamma_hat_next_list <- parallel::mclapply(seq_len(nlambda), 
                                                    function(i) {
                                                      if (i %in% zero_x_tilde) coef_zero_gamma_matrix else
                                                        as.matrix(coef(glmnet::glmnet(
                                                          x = x_tilde_list[[i]], 
                                                          y = y_tilde_list[[i]],
                                                          penalty.factor = adaptive_weights_list[[i]],
                                                          lambda = lambda_gamma_list[[i]],
                                                          standardize = F, intercept = F))[-1,,drop = F])},
                                                    mc.cores = cores)
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # update beta (main effect parameter) step 4 of algortihm in Choi et al
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    beta_hat_next_list <- beta_hat_previous_list
    
    for (j in main.effect.names) {
      
      #j = "x10"
      #print(paste(j))
      # determine the main effects not in j
      j_prime_not_in_j <- dplyr::setdiff(main.effect.names,j)
      
      
      y_tilde_2_list_temp <- lapply(beta_hat_next_list, function(i) y - 
                                      x[,j_prime_not_in_j, drop = F] %*% i[j_prime_not_in_j, , drop = F]) 
      
      #length(y_tilde_2_list_temp)
      # mclapply is faster than lapply even with just two cores
      term_2_temp_list <- parallel::mclapply(seq_len(length(beta_hat_next_list)), function(i) 
        as.matrix(
          rowSums(
            xtilde_mod(beta.main.effects = beta_hat_next_list[[i]][j_prime_not_in_j, , drop = F],
                       gamma.interaction.effects = gamma_hat_next_list[[i]],
                       interaction.names = interaction.names[-grep(j, interaction.names)], 
                       data.main.effects = x[,j_prime_not_in_j, drop = F])
          ), 
          ncol = 1), 
        mc.cores = cores)
      
      # this is of length nlambda.beta*nlambda.gamma i.e. one set of y's for each tuning parameter
      y_tilde_2_list <- mapply("-", y_tilde_2_list_temp, term_2_temp_list, SIMPLIFY = F)
      
      # index data.frame to figure out which j < j'
      index <- data.frame(main.effect.names, seq_along(main.effect.names), 
                          stringsAsFactors = F) 
      colnames(index) <- c("main.effect.names","index")
      
      # j' less than j
      j.prime.less <- index[which(index[,"index"] < index[which(index$main.effect.names == j),2]),
                            "main.effect.names"]
      
      # the if conditions in term1 and term2 are to check if there are 
      # any variables greater or less than j            
      # lapply is faster than mclapply here
      term_1_list <- if (length(j.prime.less) != 0) { 
        lapply(seq_len(length(beta_hat_next_list)), function(i) 
          x[,paste(j.prime.less,j,sep = ":")] %*% 
            (gamma_hat_next_list[[i]][paste(j.prime.less,j,sep = ":"),, drop = F] * 
               beta_hat_next_list[[i]][j.prime.less, , drop = F]))} else 
                 matrix(rep(0,length(beta_hat_next_list)), ncol = 1)
      
      # j' greater than j
      j.prime.greater <- index[which(index[,"index"] > index[which(index$main.effect.names == j),2]),
                               "main.effect.names"]
      
      term_2_list <- if (length(j.prime.greater) != 0) {
        lapply(seq_len(length(beta_hat_next_list)), function(i) 
          x[,paste(j,j.prime.greater,sep = ":")] %*% 
            (gamma_hat_next_list[[i]][paste(j, j.prime.greater,sep = ":"),, drop = F] * 
               beta_hat_next_list[[i]][j.prime.greater,,drop = F])) } else 
                 matrix(rep(0,length(beta_hat_next_list)), ncol = 1)
      
      length(term_1_list)
      
      # lapply is faster than mclapply
      x_tilde_2_list <- lapply(seq_len(length(term_1_list)), 
                               function(i) x[,j, drop = F] + 
                                 term_1_list[[i]] + term_2_list[[i]])
      
      # glmnet is giving weired results for this... and is slower than using my
      # soft function. use this. non-parallel version is faster
      beta_hat_next_list_j <- mapply(soft, 
                                     x = x_tilde_2_list, 
                                     y = y_tilde_2_list,
                                     lambda = lambda_beta_list,
                                     MoreArgs = list(
                                       weight = adaptive.weights[j,]), SIMPLIFY = F)
      
      # update beta_j for each tuning parameter
      for (i in seq_len(length(beta_hat_next_list_j))) { 
        beta_hat_next_list[[i]][j,] <- beta_hat_next_list_j[[i]]
      }
    }
    
    Q[m+1,] <- parallel::mcmapply(Q_theta,
                                  beta = beta_hat_next_list, 
                                  gamma = gamma_hat_next_list, 
                                  lambda.beta = lambda_beta_list, 
                                  lambda.gamma = lambda_gamma_list,
                                  MoreArgs = list(x = x, y = y, 
                                                  weights = adaptive.weights,
                                                  main.effect.names = main.effect.names, 
                                                  interaction.names = interaction.names),
                                  mc.cores = cores)
    
    #betas[,m+1,] <- beta_hat_next_list
    #gammas[,m+1,] <- gamma_hat_next_list
    
    # this is the delta for each pair of tuning parameters
    delta_vector <- abs(Q[m,] - Q[m+1,])/abs(Q[m,])
    
    #print(paste("Iteration:",m, ", Q(theta):",Q[m+1,2]))
    
    m = m + 1
    
    beta_hat_previous_list <- beta_hat_next_list
    
  }
  
  return(list(beta = beta_hat_next_list, gamma = gamma_hat_next_list, Q = Q, m = m))
  
}





#' Fit Strong Heredity Model for Multiple Lambdas trying to make a faster 
#' version so that it stops for a converged lambda pair
#' @param lambda.beta sequence of tuning parameters for beta. If NULL, this 
#'   function will automatically calculate a sequence which will be over a grid 
#'   of tuning parameters for gamma as well. If the user specifies a sequence 
#'   then this function will not automatically perform the serach over a grid.
#'   You will need to create the grid yourself e.g. repeat the lambda.gamma for
#'   each value of lambda.beta
#' @param nlambda.gamma number of tuning parameters for gamma. This needs to be 
#'   specified even for user defined inputs
#' @param nlambda.beta number of tuning parameters for beta. This needs to be 
#'   specified even for user defined inputs
#' @param nlambda total number of tuning parameters. This is important to 
#'   specify especially when a user defined sequence of tuning parameters is 
#'   set.
#' @param cores number of cores to use. this is used in the step to calculate
#'   
#' @note let glmnet choose the lambda_betas and lambda_gammas
#' @note the index of the tuning parameters is as follows. If for example there 
#'   are 10 lambda_gammas, and 20 lambda_betas, then the first lambda_gamma gets
#'   repeated 20 times. So the first twenty entries of tuning parameters 
#'   correspond to 1 lambda_gamma and the 20 lambda_betas
#' @note if the user specifies lambda.beta and lambda.gamma then they this will 
#'   not take all possible combinations of lambda.beta and lambda.gamma. It will
#'   be the first element of each as a pair, and so on. This is done on purpose 
#'   for use with the cv.shim function which uses the same lambda sequences for 
#'   each fold.
shim_multiple_faster <- function(x, y, main.effect.names, interaction.names, 
                          lambda.beta = NULL, lambda.gamma = NULL, threshold, max.iter, 
                          initialization.type = "ridge",
                          intercept=TRUE, normalize=TRUE,
                          nlambda.gamma = 20, 
                          nlambda.beta = 20,
                          nlambda = 400, 
                          cores = 2) {
        
      # x = X; y = Y; main.effect.names = main_effect_names;
      # interaction.names = interaction_names;
      # lambda.beta = NULL ; lambda.gamma = NULL
      # lambda.beta <- glmnet.object$lambda.beta ; lambda.gamma <- glmnet.object$lambda.gamma
      # threshold = 1e-5 ; max.iter = 500 ; initialization.type = "ridge";
      # nlambda.gamma = 5; nlambda.beta = 10; cores = 1;
      # intercept=TRUE; normalize=TRUE
    
    this.call = match.call()
    obj = standardize(x = x, y = y,intercept = intercept, normalize = normalize)
    x = obj$x
    y = obj$y
    bx = obj$bx
    by = obj$by
    sx = obj$sx

    if (is.null(lambda.gamma) & is.null(lambda.beta)) {
        
        tuning_params <- shim_once(x = x, y = y, 
                                   main.effect.names = main.effect.names,
                                   interaction.names = interaction.names,
                                   initialization.type = "ridge",
                                   nlambda.gamma = nlambda.gamma, 
                                   nlambda.beta = nlambda.beta)
        
        # tuning_params$lambda_gamma
        # tuning_params$lambda_beta
        
        # convert to a list. each element corresponds to a value of lambda_gamma
        lambda_gamma_list <- rep(lapply(seq_len(length(tuning_params$lambda_gamma)), 
                                        function(i) tuning_params$lambda_gamma[i]),
                                 each = nlambda.beta)
        
        lambda_beta_list <- lapply(seq_len(length(unlist(tuning_params$lambda_beta))), 
                                   function(i) unlist(tuning_params$lambda_beta)[i])
    } else {

        # convert to a list. each element corresponds to a value of lambda_gamma
        # these are already of the proper length i.e., if the user specifies
        # lambda.beta and lambda.gamma then they this will not take all possible 
        # combinations of lambda.beta and lambda.gamma. It will be the first element
        # of each as a pair, and so on. This is done on purpose for use with
        # the cv.shim function which uses the same lambda sequences for each fold...
        lambda_gamma_list <- lapply(seq_len(length(lambda.gamma)), 
                                        function(i) lambda.gamma[i])
                                 
        lambda_beta_list <- lapply(seq_len(length(unlist(lambda.beta))), 
                                   function(i) unlist(lambda.beta)[i])
    }
    
    adaptive.weights <- ridge_weights(x = x, y = y, 
                                      main.effect.names = main.effect.names, 
                                      interaction.names = interaction.names)
    adaptive.weights.mat <- replicate(nlambda,adaptive.weights, simplify = "matrix")
    rownames(adaptive.weights.mat) <- rownames(adaptive.weights)
    adaptive_weights_list <- lapply(seq_len(ncol(adaptive.weights.mat)), 
                                             function(i) adaptive.weights.mat[,i, drop = F])

    # initialization
    betas_and_alphas <- uni_fun(variables = colnames(x), x = x, y = y, 
                                include.intercept = F,
                                type = initialization.type)
    
    # this converts the alphas to gammas
    uni_start <- convert(betas_and_alphas, main.effect.names = main.effect.names, 
                         interaction.names = interaction.names)
    
    # need to create a matrix here instead of a 1 column vector  
    # dim1: # of variables, 
    # dim2: # of lambdas 
    
    beta_hat_previous <- replicate(nlambda, uni_start[main.effect.names, , drop = F], 
                                   simplify = "matrix")
    rownames(beta_hat_previous) <- main_effect_names
    
    gamma_hat_previous <- replicate(nlambda, uni_start[interaction.names, , drop = F], 
                                    simplify = "matrix")
    rownames(gamma_hat_previous) <- interaction_names
    
    # convert gamma and beta previous to lists each element corresponds to the
    # coefficients for each combination of lambda_gamma and lambda_beta
    beta_hat_previous_list <- lapply(seq_len(ncol(beta_hat_previous)), 
                                     function(i) beta_hat_previous[,i, drop = F])
    gamma_hat_previous_list <- lapply(seq_len(ncol(gamma_hat_previous)), 
                                      function(i) gamma_hat_previous[,i, drop = F])
    
    # store likelihood values at each iteration in a matrix Q
    # piping using magrittr::set_colnames is slower here
    # rows are the iterations, columns are the index of the sequence of 
    # lambda_gammas and lambda_betas
    # rows: iteration number
    # columns: tuning parameter
    Q <- matrix(nrow = max.iter+1, ncol = nlambda)
    
    # store the value of the likelihood at the 0th iteration
    Q[1,] <- parallel::mcmapply(Q_theta,
                                beta = beta_hat_previous_list, 
                                gamma = gamma_hat_previous_list, 
                                lambda.beta = lambda_beta_list, 
                                lambda.gamma = lambda_gamma_list,
                                weights = adaptive_weights_list,
                                MoreArgs = list(x = x, y = y, 
                                                main.effect.names = main.effect.names, 
                                                interaction.names = interaction.names),
                                mc.cores = cores)
    
    m = 1 # iteration counter
    delta = 1 # threshold initialization
    # to see which lambdas have converged: 0=not converged, 1=converged
    converged = rep(0, nlambda)
    y_tilde_list <- vector("list", nlambda)
    x_tilde_list <- vector("list", nlambda)
    gamma_hat_next_list <- vector("list", nlambda)
    y_tilde_2_list_temp <- vector("list", nlambda)
    term_2_temp_list <- vector("list", nlambda)
    term_1_list <- vector("list", nlambda)    
    term_2_list <- vector("list", nlambda)

    
    # for all the x_tilde in zero_x_tilde, return the following matrix with 0 for each coefficient
    # this is like a place holder.
    coef_zero_gamma_matrix <- matrix(data = 0, 
                                     nrow = length(interaction.names), 
                                     ncol = 1, 
                                     dimnames = list(interaction.names))
    
    # index data.frame to figure out which j < j'
    index <- data.frame(main.effect.names, seq_along(main.effect.names), 
                        stringsAsFactors = F) 
    colnames(index) <- c("main.effect.names","index")
    
    while (any(converged==0) && m < max.iter){
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # update gamma (interaction parameter)
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        # this is a nsubjects x lambda matrix for each tuning parameter stored in a list
        # each element of the list corresponds to a tuning parameter      
        # need to keep y_tilde_list and x_tilde_list of length nlambda
        
        not_converged <- which(converged==0)
        for (j in not_converged) {
            y_tilde_list[[j]] <- y - x[,main.effect.names,drop = F] %*% beta_hat_previous_list[[j]]
        }
        
        for (j in not_converged) {
            
            x_tilde_list[[j]] <- xtilde(interaction.names = interaction.names, 
                                        data.main.effects = x[,main.effect.names, drop = F],
                                        beta.main.effects = beta_hat_previous_list[[j]])
        }
        

        # indices of the x_tilde matrices that have all 0 columns
        zero_x_tilde <- which(sapply(x_tilde_list, 
                                     function(i) is.null(colnames(check_col_0(i)))))
        
        # this will store the results but will be shorter than nlambda
        gamma_hat_next_list_not_converged <- parallel::mclapply(seq_len(nlambda)[not_converged], 
            function(i) {
                if (i %in% zero_x_tilde) coef_zero_gamma_matrix else
                    as.matrix(coef(glmnet::glmnet(
                        x = x_tilde_list[[i]], 
                        y = y_tilde_list[[i]],
                        penalty.factor = adaptive_weights_list[[i]][interaction.names,,drop=F],
                        lambda = lambda_gamma_list[[i]],
                        standardize = F, intercept = F))[-1,,drop = F])},
            mc.cores = cores)

        gamma_hat_next_list <- replace(gamma_hat_next_list, not_converged, 
                                       gamma_hat_next_list_not_converged)
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # update beta (main effect parameter) step 4 of algortihm in Choi et al
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        beta_hat_next_list <- beta_hat_previous_list
        
        for (j in main.effect.names) {
            
            #j = "x5"
            #print(paste(j))
            # determine the main effects not in j
            j_prime_not_in_j <- dplyr::setdiff(main.effect.names,j)
            
            for (notconverged in not_converged) {
                y_tilde_2_list_temp[[notconverged]] <- y - 
                    x[,j_prime_not_in_j, drop = F] %*% 
                    beta_hat_next_list[[notconverged]][j_prime_not_in_j, , drop = F]
            }
            
            # length(y_tilde_2_list_temp)
            # mclapply is faster than lapply even with just two cores
            term_2_temp_list_not_converged <- parallel::mclapply(seq_len(nlambda)[not_converged], 
            function(i) 
                as.matrix(
                    rowSums(
                        xtilde_mod(beta.main.effects = beta_hat_next_list[[i]][j_prime_not_in_j, , drop = F],
                                   gamma.interaction.effects = gamma_hat_next_list[[i]],
                                   interaction.names = interaction.names[-grep(j, interaction.names)],
                                   data.main.effects = x[,j_prime_not_in_j, drop = F])
                    ), 
                    ncol = 1), 
            mc.cores = cores)
            
            term_2_temp_list <- replace(term_2_temp_list, not_converged, term_2_temp_list_not_converged)
            
            # this is of length nlambda.beta*nlambda.gamma i.e. one set of y's for each tuning parameter
            y_tilde_2_list <- mapply("-", y_tilde_2_list_temp, term_2_temp_list, SIMPLIFY = F)
            
            
            
            # j' less than j
            j.prime.less <- index[which(index[,"index"] < index[which(index$main.effect.names == j),2]),
                                  "main.effect.names"]
            
            # the if conditions in term1 and term2 are to check if there are 
            # any variables greater or less than j            
            # lapply is faster than mclapply here
            term_1_list_not_converged <- if (length(j.prime.less) != 0) { 
                lapply(seq_len(nlambda)[not_converged], function(i) 
                    x[,paste(j.prime.less,j,sep = ":")] %*% 
                        (gamma_hat_next_list[[i]][paste(j.prime.less,j,sep = ":"),, drop = F] * 
                             beta_hat_next_list[[i]][j.prime.less, , drop = F]))} else 
                                 matrix(rep(0,length(beta_hat_next_list[not_converged])), ncol = 1)
            
            term_1_list <- replace(term_1_list, not_converged, term_1_list_not_converged)
            
            # j' greater than j
            j.prime.greater <- index[which(index[,"index"] > 
                                               index[which(index$main.effect.names == j),2]),
                                     "main.effect.names"]
            
            term_2_list_not_converged <- if (length(j.prime.greater) != 0) {
                lapply(seq_len(nlambda)[not_converged], function(i) 
                    x[,paste(j,j.prime.greater,sep = ":")] %*% 
                        (gamma_hat_next_list[[i]][paste(j, j.prime.greater,sep = ":"),, drop = F] * 
                             beta_hat_next_list[[i]][j.prime.greater,,drop = F])) } else 
                                 matrix(rep(0,length(beta_hat_next_list[not_converged])), ncol = 1)
            
            
            term_2_list <- replace(term_2_list, not_converged, term_2_list_not_converged)
            
            
            # lapply is faster than mclapply
            x_tilde_2_list <- lapply(seq_len(length(term_1_list)), 
                                     function(i) x[,j, drop = F] + 
                                         term_1_list[[i]] + term_2_list[[i]])
            
            # glmnet is giving weired results for this... and is slower than using my
            # soft function. use this. non-parallel version is faster
            # the result of this should give 1 beta for each tuningn parameter
            # This calculates for all tuning parameters
            beta_hat_next_list_j <- lapply(seq_len(length(x_tilde_2_list)), function(i) 
                                            soft(x = x_tilde_2_list[[i]], 
                                           y = y_tilde_2_list[[i]],
                                           weight = adaptive_weights_list[[i]][j,,drop=F],
                                           lambda = lambda_beta_list[[i]]))
            
            # update beta_j for each tuning parameter but only those that
            # have not converged
            for (i in seq_len(nlambda)[not_converged]) { 
                beta_hat_next_list[[i]][j,] <- beta_hat_next_list_j[[i]]
            }
        }
        
        Q[m+1,not_converged] <- parallel::mcmapply(Q_theta,
                                      beta = beta_hat_next_list[not_converged], 
                                      gamma = gamma_hat_next_list[not_converged], 
                                      lambda.beta = lambda_beta_list[not_converged], 
                                      lambda.gamma = lambda_gamma_list[not_converged],
                                      weights = adaptive_weights_list[not_converged],
                                      MoreArgs = list(x = x, y = y, 
                                                      main.effect.names = main.effect.names, 
                                                      interaction.names = interaction.names),
                                      mc.cores = cores)
        

        delta <- abs(Q[m,] - Q[m+1,])/abs(Q[m,])
        # if delta is NA, this means Q wasnt calculated for the previous iteration
        # because the algorithm converged, therefore replace with threshold
        # so that it stays as converged
        delta[is.na(delta)] <- threshold 
        converged <- as.numeric(delta<=threshold)
        print(paste("Iteration:",m, ", Q(theta):",Q[m+1,2]))
        print(converged)
        #print(paste(converged))
        
        m = m + 1
        
        beta_hat_previous_list <- beta_hat_next_list
        
        # adaptive weight for each tuning parameter. currently this is the
        # same for iterations, but I am coding it here
        # for flexibility in case we want to change the weights at each iteration
        # adaptive_weights_list <- lapply(x_tilde_list, function(i)
        #     adaptive.weights[colnames(i),,drop = F])
        # 
        # colnames(x_tilde_list[[1]])
        
        adaptive_weights_list_not_converged <- lapply(seq_len(nlambda)[not_converged], 
                    function(i) 
                        update_weights(betas = beta_hat_previous_list[[i]],
                                       gammas = gamma_hat_previous_list[[i]],
                                       main.effect.names = main.effect.names, 
                                       interaction.names = interaction.names))
        
        adaptive_weights_list <- replace(adaptive_weights_list, not_converged, 
                                         adaptive_weights_list_not_converged)
        
        
    }
    
    # convert to original scale
    betas_original_scale_list <- lapply(beta_hat_next_list, function(i) i/sx[main.effect.names])
    gammas_original_scale_list <- lapply(gamma_hat_next_list, function(i) i/sx[interaction.names])
    
    
    # convert gammas to alphas
    betas_alphas_original_scale <- mapply(convert2, 
           beta = betas_original_scale_list, 
           gamma = gammas_original_scale_list, 
           MoreArgs = list(main.effect.names = main.effect.names, 
                           interaction.names = interaction.names))

    dimnames(betas_alphas_original_scale) <- list(c(main.effect.names, interaction.names),
                                   paste0("s",1:nlambda))
    
    betas_original_scale <- betas_alphas_original_scale[main.effect.names,]
    alphas_original_scale <- betas_alphas_original_scale[interaction.names,]
    
    b0 <- vector(length = nlambda)
    for (lam in seq_len(nlambda)) {
        b0[lam] <- by - sum(betas_original_scale[,lam,drop = F] * bx[main.effect.names]) - 
            sum(alphas_original_scale[,lam,drop=F]*bx[interaction.names]) 
    }
    names(b0) <- paste0("s",1:nlambda)
    
    gamma_final <- as(matrix(unlist(gammas_original_scale_list, use.names = F),
                    ncol = nlambda,
                    byrow = T,
                    dimnames = list(interaction.names, paste0("s",1:nlambda))),
                "dgCMatrix")
    beta_final = as(betas_original_scale,"dgCMatrix")
    alpha_final = as(alphas_original_scale,"dgCMatrix")

    lambda.beta = unlist(lambda_beta_list)
    names(lambda.beta) = paste0("s",1:nlambda,".beta")
    lambda.gamma = unlist(lambda_gamma_list)
    names(lambda.gamma) = paste0("s",1:nlambda, ".gamma")
    
    tuning.parameters <- matrix(nrow = 2, ncol = nlambda, 
                                dimnames = list(c("lambda.beta", "lambda.gamma"), 
                                                paste0("s",1:nlambda)))
    
    tuning.parameters["lambda.beta",] <- lambda.beta
    tuning.parameters["lambda.gamma",] <- lambda.gamma
    
    out = list(b0 = b0,
               beta = beta_final, 
               alpha = alpha_final, 
               gamma = gamma_final,
               lambda.beta = lambda.beta,
               lambda.gamma = lambda.gamma,
               tuning.parameters = tuning.parameters,
               Q = Q, m = m,
               dfbeta = nonzero(beta_final),
               dfalpha = nonzero(alpha_final),
               converged = converged, x=x,y=y,bx=bx,by=by,sx=sx,
               intercept = intercept, normalize = normalize,call=this.call,
               nlambda.gamma = nlambda.gamma,
               nlambda.beta = nlambda.beta,
               nlambda = nlambda, 
               interaction.names = interaction.names,
               main.effect.names = main.effect.names) 
    class(out) = "shim"
    return(out)
    

}



#' @description this function only works for tuning parameter values defined by 
#'   the shim_multiple_faster function. The interpolation feature is not working
#'   yet
#' @param s index of tuning parameter. Must be a character and an element of
#'   "s1","s2",...."s100", where "s100" is the index of the last pair of tuning
#'   parameters
#'   @param object is of class shim


predict.shim <- function(object, newx, 
                          #s.beta = NULL, 
                          #s.gamma = NULL,
                          s = NULL,
                          type = c("link", "response", "coefficients", 
                                   "nonzero", "class")) {
  
  #type = "coefficients"
  #object = res2
  #nlambda = length(object$lambda.beta)
  #lambda.beta = object$lambda.beta
  #lambda.gamma = object$lambda.gamma
  #names(lambda.gamma) <- paste0("s",1:nlambda)
  #names(lambda.beta) <- paste0("s",1:nlambda)
    
  type = match.arg(type)

  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE))
      stop("You need to supply a value for 'newx'")
  }
  
  # length.s.beta <- length(s.beta)
  # length.s.gamma <- length(s.gamma)
  # 
  # if (length.s.beta != length.s.gamma) stop("Length of s.beta must be equal to length of s.gamma")
  
#   if (!is.null(s)) {
#     if (any(s %ni% paste0("s",1:nlambda))) stop(paste("s must be one of ",paste0("s",1:nlambda, collapse = ",") ))
#   }
  
  # object <- c()
  # object$beta <- cvfit$glmnet.fit$beta
  # object$alpha <- cvfit$glmnet.fit$alpha
  # object$gamma <- cvfit$glmnet.fit$gamma
  # object$b0 <- cvfit$glmnet.fit$b0
  # object$lambda.beta <- cvfit$lambda.beta
  # object$lambda.gamma <- cvfit$lambda.gamma
  # object$interaction.names <- interaction_names
  # object$main.effect.names <- main_effect_names
  # object$nlambda <- cvfit$glmnet.fit$nlambda
  # s.gamma = cvfit$lambda.min.gamma ; s.beta = cvfit$lambda.min.beta
  # object <- cvfit$glmnet.fit

  a0 = t(as.matrix(object$b0))
  rownames(a0) = "(Intercept)"
  # this includes tuning parameters pairs that didnt converge
  nbeta = rbind(a0, object$beta, object$alpha)
  nbeta@Dimnames <- list(X = c("(Intercept)",object$main.effect.names, object$interaction.names),
                         Y = paste0("s",seq_len(object$nlambda)))  
  
  # this is the default returned by coef.shim i.e. any object of class shim
  # it will return all tuning parameters (including those that didnt converge)
  if (type == "coefficients" && is.null(s)) {
    return(nbeta)
  }
  
  if (type == "coefficients" && !is.null(s)) {
    return(nbeta[ , s, drop = F])
  }
  
  if (type == "nonzero") {
    nbeta = rbind(a0, object$beta, object$alpha)
    return(list(main = nonzero(nbeta[object$main.effect.names, , drop = FALSE], bystep = TRUE),
                interaction = nonzero(nbeta[object$interaction.names, , drop = FALSE], bystep = TRUE)))
  }
  
  if (inherits(newx, "sparseMatrix")) {   
    newx = as(newx, "dgCMatrix")
  }
  
  # this is used by the cv_lspath function to calculate predicted values
  # which will subsequently be used for calculating MSE for each fold
  if (type == "link") {
   
    nfit = as.matrix(cbind2(1, newx) %*% nbeta) 
    
    return(nfit)
  }
  
}

coef.shim <- function (object, s = NULL) {
    predict(object, s = s, type = "coefficients")
}

#' @param object object of class cv.shim from cv.shim function 
coef.cv.shim <- function(object, 
                         #s.beta = c("lambda.1se.beta", "lambda.min.beta"), 
                         #s.gamma = c("lambda.1se.gamma", "lambda.min.gamma")
                         s = c("lambda.1se", "lambda.min"), ...) {
  
  if (is.numeric(s) || s %ni% c("lambda.1se", "lambda.min")) stop("s must be in lambda.1se or lambda.min") 
  
  s = match.arg(s)
  
  lambda = names(which(object$glmnet.fit$tuning.parameters["lambda.beta",]==cvfit[[paste0(s,".beta")]] & 
                         object$glmnet.fit$tuning.parameters["lambda.gamma",]==cvfit[[paste0(s,".gamma")]]))
  
  coef(object$glmnet.fit, s = lambda, ...)
}

# coef(res2, s=c("s10","s11"))
# predict(res2, type = "nonzero")
# predict(res2, type = "coefficients", s = "s10")
# predict(res2, type = "link", newx = X)


plot.shim <- function(x, xvar = c("norm", "lambda", "dev"), label = T, 
          ...) {
    xvar = match.arg(xvar)
    plotShim(x$beta, 
             # lambda = x$lambda, 
             df = x$dfbeta, 
             # dev = x$dev.ratio, 
             label = label, 
             xvar = xvar, ...)
}


plotShim <- function(beta, norm, 
                      #lambda, 
                      df, 
                      #dev, 
                      label = T, 
                       xvar = c("norm", 
                                "lambda", "dev"), xlab = iname, ylab = "Coefficients", ...) {
    #beta <- res2$beta
    which = nonzero(beta)
    nwhich = length(which)
    switch(nwhich + 1, `0` = {
        warning("No plot produced since all coefficients zero")
        return()
    }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
    beta = as.matrix(beta[which, , drop = FALSE])
    xvar = match.arg(xvar)
    switch(xvar, norm = {
        index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
        iname = "L1 Norm"
        approx.f = 1
    }, lambda = {
        index = log(lambda)
        iname = "Log Lambda"
        approx.f = 0
    }, dev = {
        index = dev
        iname = "Fraction Deviance Explained"
        approx.f = 1
    })
    dotlist = list(...)
    type = dotlist$type
    if (is.null(type)) 
        matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, 
                type = "l", ...)
    else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, 
                 ...)
    atdf = pretty(index)
    # prettydf = approx(x = index, y = df, xout = atdf, rule = 2, 
    #                   method = "constant", f = approx.f)$y
    # axis(3, at = atdf, labels = prettydf, tcl = NA)
    if (label) {
        nnz = length(which)
        xpos = max(index)
        pos = 4
        if (xvar == "lambda") {
            xpos = min(index)
            pos = 2
        }
        xpos = rep(xpos, nnz)
        ypos = beta[, ncol(beta)]
        text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
    }
}

nonzero <- function(beta, bystep = FALSE) {
    beta <- as.matrix(beta)
    nr = nrow(beta)
    if (nr == 1) {
        if (bystep) 
            apply(beta, 2, function(x) if (abs(x) > 0) 
                1
                else NULL)
        else {
            if (any(abs(beta) > 0)) 
                1
            else NULL
        }
    }
    else {
        beta = abs(beta) > 0
        which = seq(nr)
        ones = rep(1, ncol(beta))
        nz = as.vector((beta %*% ones) > 0)
        which = which[nz]
        if (bystep) {
            if (length(which) > 0) {
                beta = as.matrix(beta[which, , drop = FALSE])
                nzel = function(x, which) if (any(x)) 
                    which[x]
                else NULL
                which = apply(beta, 2, nzel, which)
                if (!is.list(which)) 
                    which = data.frame(which)
                which
            }
            else {
                dn = dimnames(beta)[[2]]
                which = vector("list", length(dn))
                names(which) = dn
                which
            }
        }
        else which
    }
}




cv.shim <- function(x, y, main.effect.names, interaction.names, 
                     lambda.beta = NULL, lambda.gamma = NULL, threshold, max.iter, 
                     initialization.type = "ridge",
                     intercept=TRUE, normalize=TRUE,
                     nlambda.gamma = 5, 
                     nlambda.beta = 20,
                     nlambda = 100,
                     cores = 1,
                     #x, y, offset = NULL, 
                     #lambda = NULL,
                     weights, 
                     type.measure = c("mse", "deviance", "class", "auc", "mae"), 
                     nfolds = 10, 
                     foldid, 
                     grouped = TRUE, keep = FALSE, parallel = TRUE, ...) {
    
    # x = X; y = Y; main.effect.names = main_effect_names;
    # interaction.names = interaction_names;
    # lambda.beta = NULL ; lambda.gamma = NULL
    # threshold = 1e-4 ; max.iter = 500 ; initialization.type = "ridge";
    # nlambda.gamma = 5; nlambda.beta = 20; cores = 1;
    # intercept=TRUE; normalize=TRUE
    # 
    # nfolds = 5
    # grouped = TRUE; keep = FALSE; parallel = TRUE
    
    if (missing(type.measure)) 
        type.measure = "default"
    else type.measure = match.arg(type.measure)
    if (!is.null(lambda.beta) && length(lambda.beta) < 2) 
        stop("Need more than one value of lambda.beta for cv.shim")
    if (!is.null(lambda.gamma) && length(lambda.gamma) < 2) 
        stop("Need more than one value of lambda.gamma for cv.shim")
    N = nrow(x)
    if (missing(weights)) 
        weights = rep(1, N)
    else weights = as.double(weights)
    y = drop(y)
    glmnet.call = match.call(expand.dots = TRUE)
    which = match(c("type.measure", "nfolds", "foldid", "grouped", 
                    "keep"), names(glmnet.call), F)
    if (any(which)) 
        glmnet.call = glmnet.call[-which]
    glmnet.call[[1]] = as.name("glmnet")
    glmnet.object = shim_multiple_faster(x = x, y = y, 
                                         main.effect.names = main.effect.names,
                                         interaction.names = interaction.names,
                                         lambda.beta = lambda.beta, lambda.gamma = lambda.gamma,
                                         nlambda = nlambda,
                                         threshold = threshold, max.iter = max.iter, 
                                         initialization.type = initialization.type,
                                         nlambda.gamma = nlambda.gamma, 
                                         nlambda.beta = nlambda.beta, cores = 1)
    
    glmnet.object$call = glmnet.call
    #glmnet.object$nlambda.beta
    nz.main = sapply(predict(glmnet.object, type = "nonzero")[["main"]], length)
    nz.interaction = sapply(predict(glmnet.object, type = "nonzero")[["interaction"]], length)
    
    if (missing(foldid)) 
        foldid = sample(rep(seq(nfolds), length = N)) else nfolds = max(foldid)
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist = as.list(seq(nfolds))
    if (parallel) {
        outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar% 
        {
            which = foldid == i
            if (is.matrix(y)) y_sub = y[!which, ] else y_sub = y[!which]
            # glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda, 
            #        offset = offset_sub, weights = weights[!which], 
            #        ...)
            print(paste("Foldid = ",i))
            shim_multiple_faster(x = x[!which, , drop = FALSE], 
                                 y = y_sub, 
                                 #x = x , y = y,
                                 main.effect.names = main.effect.names,
                                 interaction.names = interaction.names,
                                 lambda.beta = glmnet.object$lambda.beta, 
                                 lambda.gamma = glmnet.object$lambda.gamma,
                                 nlambda = glmnet.object$nlambda,
                                 threshold = threshold, 
                                 max.iter = max.iter , 
                                 initialization.type = initialization.type,
                                 nlambda.gamma = nlambda.gamma, 
                                 nlambda.beta = nlambda.beta, cores = 1)
        }
    }
    else {
        for (i in seq(nfolds)) {
            which = foldid == i
            if (is.matrix(y)) 
                y_sub = y[!which, ]
            else y_sub = y[!which]
            if (is.offset) 
                offset_sub = as.matrix(offset)[!which, ]
            else offset_sub = NULL
            outlist[[i]] = shim_multiple_faster(x[!which, , drop = FALSE], y_sub, 
                                                main.effect.names,
                                                interaction.names,
                                                lambda.beta = glmnet.object$lambda.beta, 
                                                lambda.gamma = glmnet.object$lambda.gamma,
                                                nlambda = glmnet.object$nlambda,
                                                threshold, max.iter , initialization.type,
                                                nlambda.gamma, nlambda.beta, cores = 1)
        }
    }
    
    #outlist[[1]]
    
    #fun = paste("cv", class(glmnet.object)[[1]], sep = ".")
    lambda.beta = glmnet.object$lambda.beta
    lambda.gamma = glmnet.object$lambda.gamma
    
    cvstuff = do.call(cv_lspath, list(outlist = outlist, 
                                      x = x, y = y, foldid = foldid, 
                                      nlambda = glmnet.object$nlambda,
                                      nlambda.beta = glmnet.object$nlambda.beta,
                                      nlambda.gamma = glmnet.object$nlambda.gamma))
    
    cvm = cvstuff$cvm
    cvsd = cvstuff$cvsd
    # the following checks is any of the tunining parameter pairs
    # have a cvsd==NA or did not converge. cvstuff$converged should be
    # equal to the number of folds because it is the sum of booleans 
    # that converged over all folds.
    nas = (is.na(cvsd) + (cvstuff$converged != nfolds)) != 0
    if (any(nas)) {
        lambda.beta = lambda.beta[!nas]
        lambda.gamma = lambda.gamma[!nas]
        cvm = cvm[!nas]
        cvsd = cvsd[!nas]
        # this is the total number of non-zero parameters (both betas and alphas)
        nz.main = nz.main[!nas]
        nz.interaction = nz.interaction[!nas]
    }
    cvname = cvstuff$name
    
    df <- as.data.frame(cbind(lambda.beta = lambda.beta, 
          lambda.gamma = lambda.gamma,
          mse = cvm,
          upper = cvm+cvsd,
          lower = cvm-cvsd,
          nz.main = nz.main,
          nz.interaction = nz.interaction,
          lg = round(log(lambda.gamma),2)))
    
    out = list(lambda.beta = lambda.beta, lambda.gamma = lambda.gamma, 
               cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd, 
               cvlo = cvm - cvsd, nz.main = nz.main, name = cvname,
               nz.interaction = nz.interaction,
               glmnet.fit = glmnet.object, converged = cvstuff$converged, cvm.mat.all = cvstuff$cvm.mat.all,
               df = df)
    # if (keep) 
    #     out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
    lamin.beta = if (cvname == "AUC") 
        getmin(lambda.beta, -cvm, cvsd, type = "beta") else getmin(lambda.beta, cvm, cvsd, type = "beta")
    lamin.gamma = if (cvname == "AUC") 
         getmin(lambda.gamma, -cvm, cvsd, type = "gamma") else getmin(lambda.gamma, cvm, cvsd, type = "gamma")
    obj = c(out, as.list(lamin.beta), as.list(lamin.gamma)) #, as.list(glmnet.object$lambda.beta),
            #as.list(glmnet.object$lambda.gamma))
    class(obj) = "cv.shim"
    obj
}

#' @param nlambda total number of tuning parameter pairs. This includes those 
#'   pairs of tuning parameters that didn't converge in the CV folds. The output
#'   of this function only returns values for those tuning paramters that DID 
#'   converge
#' @param outlist list containing results from cv.shim function. each element of
#'   the list is a run of shim_multiple_faster for each CV fold

cv_lspath <- function(outlist, x, y, foldid,
                      nlambda, nlambda.beta, nlambda.gamma) {
    #   typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    #   if (pred.loss == "default") 
    #     pred.loss <- "loss"
    #   if (!match(pred.loss, c("loss"), FALSE)) {
    #     warning("Only 'loss' available for least squares regression; 'loss' used")
    #     pred.loss <- "loss"
    #   }
    
    # lambda.beta = glmnet.object$lambda.beta
    # lambda.gamma = glmnet.object$lambda.gamma
    # nlambda = glmnet.object$nlambda
    
    y <- as.double(y)
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), nlambda)
    nlams <- double(nfolds)
    converged <- matrix(nrow = nfolds, ncol = nlambda)
    for (i in seq(nfolds)) {
        #i=1
        which <- foldid == i
        # this gives be the fitted object for each CV fold
        fitobj <- outlist[[i]]
        
        # this gives the predicted responses for the subjects in the held-out fold for each lambda
        # so if each fold has 20 subjects, and there are 100 lambdas, then this will return a 
        # 20 x 100 matrix
        #preds <- x[which, , drop = FALSE]  %*% rbind2(fitobj$beta, fitobj$alpha)
        preds <- predict(fitobj, newx = x[which, ,drop = F], type = "link")
        #preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
        nlami <- fitobj$nlambda
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
        converged[i,] <- fitobj$converged
    }
    
    conv <- colSums(converged)
    cvraw <- (y - predmat)^2
    cvob <- cvcompute(cvraw, foldid, nlams)
    cvraw <- cvob$cvraw
    N <- cvob$N
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvm_mat_all <- matrix(cvm, ncol = nlambda.gamma, nrow = nlambda.beta, byrow = T)

    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
    list(cvm = cvm, cvsd = cvsd, name = "MSE", converged = conv, cvm.mat.all = cvm_mat_all)
}

cvcompute <- function(mat, foldid, nlams) {
    nfolds <- max(foldid)
    outmat <- matrix(NA, nfolds, ncol(mat))
    good <- matrix(0, nfolds, ncol(mat))
    mat[is.infinite(mat)] <- NA
    for (i in seq(nfolds)) {
        mati <- mat[foldid == i, ]
        outmat[i, ] <- apply(mati, 2, mean, na.rm = TRUE)
        good[i, seq(nlams[i])] <- 1
    }
    N <- apply(good, 2, sum)
    list(cvraw = outmat, N = N)
}


getmin <- function (lambda, cvm, cvsd, type) {
    cvmin = min(cvm, na.rm = TRUE)
    idmin = cvm <= cvmin
    lambda.min = max(lambda[idmin], na.rm = TRUE)
    idmin = match(lambda.min, lambda)
    semin = (cvm + cvsd)[idmin]
    idmin = cvm <= semin
    lambda.1se = max(lambda[idmin], na.rm = TRUE)
    # this is to get the index of the tuning parameter pair which is labelled by "s"
    # e.g. "s25" corresponds to the 25th pair of tuning parameters
    lambda.min.name = gsub(".beta","",names(lambda.min))
    lambda.1se.name = gsub(".beta","",names(lambda.1se))
    res <- list(lambda.min = lambda.min, lambda.min.name, lambda.1se = lambda.1se, lambda.1se.name )
    names(res) <- c(paste0("lambda.min.",type),"lambda.min.name", paste0("lambda.1se.",type),"lambda.1se.name")
    res
}

lambda.interp <- function (lambda, s) {
  if (length(lambda) == 1) {
    nums = length(s)
    left = rep(1, nums)
    right = left
    sfrac = rep(1, nums)
  }
  else {
    s[s > max(lambda)] = max(lambda)
    s[s < min(lambda)] = min(lambda)
    k = length(lambda)
    sfrac <- (lambda[1] - s)/(lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac = (sfrac - lambda[right])/(lambda[left] - lambda[right])
    sfrac[left == right] = 1
  }
  list(left = left, right = right, frac = sfrac)
}

plot.cv.shim <- function(x, sign.lambda = 1, ...) {
  
  pckg = try(require(ggplot2))
  if(!pckg) {
    cat("Installing 'ggplot2' from CRAN\n")
    getPckg("ggplot2")
    require(ggplot2)
  }
  
  pckg = try(require(latex2exp))
  if(!pckg) {
    cat("Installing 'latex2exp' from CRAN\n")
    getPckg("latex2exp")
    require(latex2exp)
  }
  
  pckg = try(require(data.table))
  if(!pckg) {
    cat("Installing 'data.table' from CRAN\n")
    getPckg("data.table")
    require(data.table)
  }

  cvobj = x
  
  d <- cvobj$df %>% 
    as.data.table %>% 
    mutate(lambda.min.beta = cvobj$lambda.min.beta,
           lambda.min.gamma = cvobj$lambda.min.gamma, 
           lambda.1se.beta = cvobj$lambda.1se.beta,
           lambda.1se.gamma = cvobj$lambda.1se.gamma)  
  
  # needed to get colored lines
  d2 <- d[(lambda.beta==lambda.min.beta & lambda.gamma == lambda.min.gamma) | 
            (lambda.beta==lambda.1se.beta & lambda.gamma == lambda.1se.gamma)] %>% 
    melt(measure.vars = c("lambda.min.beta","lambda.1se.beta"))
  
  d2[,variable:=gsub(".beta", "",variable)]
  
  appender <- function(string) TeX(paste("$\\log(\\lambda_{\\gamma}) = $",string))  
  
  p <- ggplot(d, 
         aes(log(lambda.beta), 
             ymin = lower, 
             ymax = upper)) 
  
  l <- ggplot_build(p)
  p + 
    geom_errorbar(color = "grey", width = 0.5) + 
    geom_point(aes(x = log(lambda.beta), y = mse), colour = "red") +
    theme_bw() + 
    ylim(c(min(d$lower) - 10 , max(d$upper) + 500)) + 
    facet_wrap(~lg, scales = "fixed",
               #switch = "x",
               labeller = as_labeller(appender, default = label_parsed)) + 
    theme(strip.background = element_blank(), 
          strip.text.x = element_text(size = rel(1.3)),
          legend.position = "bottom") + 
    xlab(TeX("$\\log(\\lambda_{\\beta})$")) + 
    geom_vline(data = d2[lambda.gamma==lambda.1se.gamma & variable == "lambda.1se"], 
               aes(xintercept = log(value), colour = variable), size = 0.7, linetype = 1) +
    geom_vline(data = d2[lambda.gamma==lambda.min.gamma & variable == "lambda.min"], 
               aes(xintercept = log(value), colour = variable),size = 0.7, linetype = 1) + 
    scale_color_discrete(name="") + 
    geom_text(aes(label = nz.main, x = log(lambda.beta), y = Inf, vjust = 1)) +
    geom_text(aes(label = nz.interaction, x = log(lambda.beta), y = Inf,
                  vjust = 2)) +
    ylab(c("10 fold CV MSE")) + 
    coord_cartesian(ylim = c(l$panel$ranges[[1]]$y.range[1], l$panel$ranges[[1]]$y.range[2]*1.1))

  # cvobj = x
  # xlab = "log(Lambda)"
  # if (sign.lambda < 0) 
  #   xlab = paste("-", xlab, sep = "")
  # plot.args = list(x = sign.lambda * log(cvobj$lambda.beta), y = cvobj$cvm, 
  #                  ylim = range(cvobj$cvup, cvobj$cvlo), 
  #                  #ylim = quantile(c(cvobj$cvlo, cvobj$cvup), probs = c(0,0.85)),
  #                  xlab = xlab, ylab = cvobj$name, 
  #                  type = "n")
  # new.args = list(...)
  # if (length(new.args)) 
  #   plot.args[names(new.args)] = new.args
  # do.call("plot", plot.args)
  # error.bars(sign.lambda * log(cvobj$lambda.beta), cvobj$cvup, 
  #            cvobj$cvlo, width = 0.01, col = "darkgrey")
  # points(sign.lambda * log(cvobj$lambda.beta), cvobj$cvm, pch = 20, 
  #        col = "red")
  # axis(side = 3, at = sign.lambda * log(cvobj$lambda.beta), labels = paste("main:",cvobj$nz.main), 
  #      tick = FALSE, line = 1)
  # axis(side = 3, at = sign.lambda * log(cvobj$lambda.beta), labels = paste("interaction:",cvobj$nz.interaction), 
  #      tick = FALSE, line = 0)
  # abline(v = sign.lambda * log(cvobj$lambda.min.beta), lty = 3)
  # abline(v = sign.lambda * log(cvobj$lambda.1se.beta), lty = 3)
  # invisible()
}

error.bars <- function (x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}






#' Fit the Strong Heredity Interactions Model with Fixed Betas
#' 
#' @description This is a test function to see if the algorithm is converging if
#'   it only needs to estimate one set of parameters. In this function we are 
#'   fixing the betas at the true betas, and trying to estimate the gammas.
#' @param fixed.beta p x 1 matrix of betas, with rows labelled accordingly
#' @note refer to \code{shim} function for argument definitions and return value

shim_fix_betas <- function(x, y, main.effect.names, interaction.names, 
                           lambda.beta, lambda.gamma, threshold, max.iter, 
                           initialization.type = "ridge", fixed.betas) {
  
  adaptive.weights <- ridge_weights(x = x, y = y, 
                                    main.effect.names = main.effect.names, 
                                    interaction.names = interaction.names)
  
  # initialization
  betas_and_alphas <- uni_fun(variables = colnames(x), x = x, y = y, 
                              include.intercept = F,
                              type = initialization.type)
  
  # this converts the alphas to gammas
  uni_start <- betas_and_alphas %>% 
    convert(., main.effect.names = main.effect.names, 
            interaction.names = interaction.names)
  
  # initialize beta_hat_next also because we are updating each beta individually
  #beta_hat_previous <- beta_hat_next <- uni_start[main.effect.names, , drop = F]
  
  beta_hat_previous <- fixed.betas
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
    
    beta_hat_next <- fixed.betas
    
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
    
  }
  
  return(list(beta = betas, gamma = gammas, Q = Q, m = m))
  
}


#' Fit the Strong Heredity Interactions Model with Fixed Gammas
#' 
#' @description This is a test function to see if the algorithm is converging if
#'   it only needs to estimate one set of parameters. In this function we are 
#'   fixing the gammas at the true gammas, and trying to estimate the betas.
#' @param fixed.gamma p(p-1)/2 x 1 matrix of betas, with rows labelled
#'   accordingly
#' @note refer to \code{shim} function for argument definitions and return value
shim_fix_gamma <- function(x, y, main.effect.names, interaction.names, 
                           lambda.beta, lambda.gamma, threshold, max.iter, 
                           initialization.type = "ridge", fixed.gamma) {
  
  adaptive.weights <- ridge_weights(x = x, y = y, 
                                    main.effect.names = main.effect.names, 
                                    interaction.names = interaction.names)
  
  # initialization
  betas_and_alphas <- uni_fun(variables = colnames(x), x = x, y = y, 
                              include.intercept = F,
                              type = initialization.type)
  
  # this converts the alphas to gammas
  uni_start <- betas_and_alphas %>% 
    convert(., main.effect.names = main.effect.names, 
            interaction.names = interaction.names)
  
  # initialize beta_hat_next also because we are updating each beta individually
  #beta_hat_previous <- beta_hat_next <- uni_start[main.effect.names, , drop = F]
  
  beta_hat_previous <- uni_start[main.effect.names, , drop = F]
  gamma_hat_previous <- fixed.gamma
  
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
    
    gamma_hat_next <- fixed.gamma
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # update beta (main effect parameter) step 4 of algortihm in Choi et al
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    beta_hat_next <- beta_hat_previous
    
    for (j in main.effect.names) {
      
      # determine the main effects not in j
      j_prime_not_in_j <- dplyr::setdiff(main.effect.names,j)
      
      #       y_tilde_2 <- y - x[,j_prime_not_in_j] %*% beta_hat_next[j_prime_not_in_j,] - 
      #         (xtilde(interaction.names = interaction.names[-grep(j, interaction.names)],
      #                 data.main.effects = x[,j_prime_not_in_j],
      #                 beta.main.effects = beta_hat_next[j_prime_not_in_j,,drop=F]) %>% 
      #            rowSums() %>% as.matrix(ncol = 1))
      
      y_tilde_2 <- y - x[,j_prime_not_in_j] %*% beta_hat_next[j_prime_not_in_j,] - 
        (xtilde_mod(interaction.names = interaction.names[-grep(j, interaction.names)],
                    data.main.effects = x[,j_prime_not_in_j],
                    beta.main.effects = beta_hat_next[j_prime_not_in_j,,drop=F],
                    gamma.interaction.effects = gamma_hat_next) %>% 
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
                                  nlambda = 1, intercept = F,
                                  lambda = lambda.beta,
                                  standardize = F,
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








#' Update Weights based on betas and gammas. Currently not being used.
#' 
#' @description uses betas and gammas to update weights. this is used to update 
#'   the weights at each iteration of the fitting algorithm in the \code{shim} 
#'   function
#' @param betas.and.gammas q x 1 data.frame or matrix of betas and gamma
#'   estimates. The rownames must be appropriately labelled because these labels
#'   are be used in this function and must match those in the arguments 
#'   \code{main.effect.names} and \code{interaction.names}
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be 
#'   separated by a ':' (e.g. x1:x2)
#' @param include.intercept logical if intercept should be fitted. default is 
#'   FALSE. Should be set to TRUE if y is not centered
#' @note Currently this is not being used in the shim function i.e. we are not 
#'   updating the weights at each iteration
#' @return q x 1 matrix of weights

update_weights <- function(betas,
                           gammas, 
                           main.effect.names, 
                           interaction.names,
                           epsilon = 1e-5) {
  
    betas.and.gammas <- rbind(betas, gammas)
    
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






#' Check Arguments Functions
#' 
#' @description function to check inputs of shim function
#' @note adapted from
#'   \link{https://github.com/selective-inference/R-software/blob/master/selectiveInference/R/funs.common.R}
#' 

checkargs.xy <- function(x, y) {
  if (missing(x)) stop("x is missing")
  if (is.null(x) || !is.matrix(x)) stop("x must be a matrix")
  if (missing(y)) stop("y is missing")
  if (is.null(y) || !is.numeric(y)) stop("y must be numeric")
  if (ncol(x) == 0) stop("There must be at least one predictor [must have ncol(x) > 0]")
  if (checkcols(x)) stop("x cannot have duplicate columns")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 0]")
  if (length(y)!=nrow(x)) stop("Dimensions don't match [length(y) != nrow(x)]")
}



checkargs.misc <- function(sigma=NULL, alpha=NULL, k=NULL,
                           gridrange=NULL, gridpts=NULL, griddepth=NULL,
                           mult=NULL, ntimes=NULL,
                           beta=NULL, lambda=NULL, tol.beta=NULL, tol.kkt=NULL,
                           bh.q=NULL) {
  
  if (!is.null(sigma) && sigma <= 0) stop("sigma must be > 0")
  if (!is.null(lambda) && lambda < 0) stop("lambda must be >= 0")
  if (!is.null(alpha) && (alpha <= 0 || alpha >= 1)) stop("alpha must be between 0 and 1")
  if (!is.null(k) && length(k) != 1) stop("k must be a single number")
  if (!is.null(k) && (k < 1 || k != floor(k))) stop("k must be an integer >= 1")
  if (!is.null(gridrange) && (length(gridrange) != 2 || gridrange[1] > gridrange[2]))
    stop("gridrange must be an interval of the form c(a,b) with a <= b")
  if (!is.null(gridpts) && (gridpts < 20 || gridpts != round(gridpts)))
    stop("gridpts must be an integer >= 20")
  if (!is.null(griddepth) && (griddepth > 10 || griddepth != round(griddepth)))
    stop("griddepth must be an integer <= 10")
  if (!is.null(mult) && mult < 0) stop("mult must be >= 0")
  if (!is.null(ntimes) && (ntimes <= 0 || ntimes != round(ntimes)))
    stop("ntimes must be an integer > 0")
  if (!is.null(beta) && sum(beta!=0)==0) stop("Value of lambda too large, beta is zero")
  if (!is.null(lambda) && length(lambda) != 1) stop("lambda must be a single number")
  if (!is.null(lambda) && lambda < 0) stop("lambda must be >=0")
  if (!is.null(tol.beta) && tol.beta <= 0) stop("tol.beta must be > 0")
  if (!is.null(tol.kkt) && tol.kkt <= 0) stop("tol.kkt must be > 0")
}

# Make sure that no two columms of A are the same
# (this works with probability one).

checkcols <- function(A) {
  b = rnorm(nrow(A))
  a = sort(t(A)%*%b)
  if (any(diff(a)==0)) return(TRUE)
  return(FALSE)
}

standardize <- function(x, y, intercept, normalize) {
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  if (intercept) {
    bx = colMeans(x)
    by = mean(y)
    x = scale(x,bx,FALSE)
    y = y-mean(y)
  } else {
    bx = rep(0,p)
    by = 0
  }
  if (normalize) {
    sx = sqrt(colSums(x^2)/n)
    x = scale(x,FALSE,sx)
  } else {
    sx = rep(1,p)
  }
  
  return(list(x=x,y=y,bx=bx,by=by,sx=sx))
}

"%ni%" <- Negate("%in%")



grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}


isEmpty <- function(x) {
  return(length(x)==0)
}
