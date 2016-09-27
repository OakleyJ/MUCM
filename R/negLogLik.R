#' @title Negative Log Likelihood Functions

#' @description These functions are minimised to find the optimal value of theta (or phi and sigma in LMC case). 
#'              \code{negLogLikNugget} is used in the univariate case for uncertain input parameters.
#'              \code{negLogLikLMC} is used when applying the LMC Emulator for Multivariate models
#' @param theta Initial values for the parameters to be optimized over. (includes values for all the parameters that need to be estimated.)
#'        For \code{negLogLik}, theta represents just phi, for \code{negLogLikNugget}, theta represents just phi and sigma combined.
#' @param inputs A data frame, matrix or vector containing the input values of the training data.
#' @param H A matrix of prior mean regressors for the training data.
#' @param outputs A data frame, matrix or vector containing the output values of the training data. In \code{negLogLikLMCOptim}, the outputs should be stacked (either a vector or a matrix with 1 column). 
#' @param cor.function Specifies a correlation function used as part of the prior information for the emulator
#' @param ... additional arguments to be passed on to correlation functions (see \code{\link{corGaussian}})
#' @param nugget  For noisy data, a vector giving the observation variance for each training data point. 
#' @return The function returns the negetive log-likelihood of \code{theta}
#' @seealso \code{\link{corGaussian}}
#' @author Sajni Malde, Jeremy Oakley


# \prod_{i=1}^{n.theta} \exp[-\{(x_i-x_i')/\delta_i}^2] where delta = correlation length parameters -ADD IN DESC
#' @export
negLogLik <- function(theta, nugget = NULL, inputs, H, outputs, cor.function, ...) {
    
    n <- nrow(inputs)
    n.regressors <- ncol(H)
    
    A <- cor.function(inputs, phi = theta, ...)
    
    # ensures A=c(inputs,inputs) positive definite
    L <- try(chol(A), silent = TRUE)
    if (class(L) == "try-error")
        return(negloglik.lim)
    
    # Calculate Q=H^T A^{-1}H via A=LL^T
    # tL <- t(L)
    w <- try(backsolve(L, H, transpose = TRUE))     # = solve(tL, H)
    if (class(w) == "try-error")
        return(negloglik.lim)
    Q <- crossprod(w)                               # = t(w) %*% w
    
    # Calculate A^{-1}outputs via A=LL^T
    mat1 <- backsolve(L, outputs, transpose = TRUE) 
    mat2 <- backsolve(L, mat1)                      # = solve(L, mat1)
    
    # Check that Q is positive definite and Calculate H^TA^{-1}
    K <- try(chol(Q)) 
    if (class(K) == "try-error") 
        return(negloglik.lim)
    mat3 <- backsolve(K, t(H), transpose = TRUE) 
    mat4 <- mat3 %*% mat2
    betahat <- backsolve(K, mat4) 
    Hbeta <- H %*% betahat
    out.min.Hbeta <- outputs - Hbeta
    w2 <- backsolve(L, out.min.Hbeta, transpose = TRUE) 
    det.sigmasq.hat.prop <- det(crossprod(w2))  # constant ((n - n.regressors - 2) ^ (-1)) in sigmahat ignored
    
    n.outputs <- ncol(outputs)  # where n.outputs is number of outputs
    negloglik <- n.outputs * sum(log(diag(L))) + n.outputs * sum(log(diag(K))) + 0.5 * (n - n.regressors) * log(det.sigmasq.hat.prop)
    
    if (missing(negloglik) || (negloglik > negloglik.lim) || is.infinite(negloglik) || is.na(negloglik)) 
        return(negloglik.lim) 
    
    negloglik
}



#########---------- NEg Log Likelihood - Nugget version! ---------- ###########

#' @rdname negLogLik
#' @export
negLogLikNugget <- function(theta, nugget, inputs, H, outputs, cor.function, ...) {
    # extract number of inputs and regressors
    n <- nrow(inputs)
    n.regressors <- ncol(H)
    
    # extract phi and sigmasq from theta
    phi <- theta[1:(length(theta) - 1)]
    log.sigmasq <- tail(theta, 1)
    sigmasq <- exp(log.sigmasq)
    
    # calculate correlation matrix
    A <- cor.function(inputs, phi = phi, ...) + diag(nugget) / sigmasq
    
    # return if any infinite values
    if (any(is.infinite(A)))
        return(negloglik.lim)
    
    # ensures A=c(inputs,inputs) positive definite
    L <- try(chol(A), silent = TRUE)
    if (class(L) == "try-error") 
        return(negloglik.lim)
    
    # Calculate Q=H^T A^{-1}H via A=LL^T
    w <- try(backsolve(L, H, transpose = TRUE))  # = solve(tL, H)
    if (class(w) == "try-error") 
        return(negloglik.lim)
    Q <- crossprod(w)                            # = t(w) %*% w
    
    # Calculate A^{-1}outputs via A=LL^T
    mat1 <- backsolve(L, outputs, transpose = TRUE)  
    mat2 <- backsolve(L, mat1)
    
    # Check that Q is positive definite and Calculate H^TA^{-1}
    K <- try(chol(Q))  
    if (class(K) == "try-error") 
        return(negloglik.lim)
    mat3 <- backsolve(K, t(H), transpose = TRUE)  
    mat4 <- mat3 %*% mat2
    betahat <- backsolve(K, mat4)  
    Hbeta <- H %*% betahat
    out.min.Hbeta <- outputs - Hbeta
    w2 <- backsolve(L, out.min.Hbeta, transpose = TRUE)
    
    # calculate log likelihood
    negloglik <- (n + 2 - n.regressors) * log.sigmasq + sum(log(diag(L))) + sum(log(diag(K))) + crossprod(w2) / sigmasq
    
    if (missing(negloglik) || (negloglik > negloglik.lim) || is.infinite(negloglik) || is.na(negloglik)) 
        return(negloglik.lim)
    
    negloglik
}
