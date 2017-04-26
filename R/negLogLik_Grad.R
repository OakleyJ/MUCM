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

#' @export
negLogLik_Grad <- function(theta, inputs, H, outputs, cor.function, nugget = NULL) {
    
    # Returns the derivative of -log likelihood, with respect to phi = log (delta),
    # assuming Gaussian correlation function
    # c(x, x') = exp[-{(x-x') / delta}^2]
    # Note: original code differentiatied w.r.t. 2*phi, so end result is multiplied by 2
    
    ncol.inputs <- ncol(inputs)
    n <- nrow(inputs)
    n.regressors <- ncol(H)
    # negative log likelihood of 2*log(roughness parameter), the
    # transformed roughness parameters

    #part of makeAdat
    #Phi <- diag(1/exp(theta[1:ncol.inputs]/2)) # if theta = 2 * log delta
    Phi <- diag(1/exp(theta[1:ncol.inputs])) # if theta = log delta
    nug <- 0 # 1/(1 + exp(-theta[ncol.inputs + 1]))
    inputs.phi <- inputs %*% Phi
    Dk <- lapply(1:ncol.inputs, function(k) doA(inputs.phi[, k, drop = FALSE]))
    D <- matrix(rowSums(sapply(Dk, function(x) x)), n)
    A <- (1 - nug) * exp(-D)
    diag(A) <- 1

    # ensures A=c(inputs,inputs) positive definite
    L <- try(chol(A), silent = TRUE)
    if (class(L) == "try-error")
        return(negloglik.lim)

    # # Calculate Q=H^T A^{-1}H via A=LL^T
    w <- try(backsolve(L, H, transpose = TRUE))     # = solve(tL, H)
    if (class(w) == "try-error")
        return(negloglik.lim)

    Q <- crossprod(w)                               # = t(w) %*% w
    # 
    # calculating the gradient
    iA <- chol2inv(L)
    P <- iA - iA %*% H %*% solve(Q, crossprod(H, iA))
    P.outputs <- P %*% outputs
    R <- P - P.outputs %*% solve(crossprod(outputs, P.outputs), crossprod(outputs, P))
    # R2 <- P - tcrossprod(P.outputs)/(n - n.regressors - 2)/sigmahat.prop[1, 1]
    # R - R2
    
    gradA <- lapply(1:ncol.inputs, function(i) Phi[i, i] * Dk[[i]] * A)
    
    # # the extra bit for the nugget
    # gradA[[ncol.inputs + 1]] <- -0.5 * exp(-D)
    # diag(gradA[[ncol.inputs + 1]]) <- 0
    
    gout <- sapply(gradA, function(x) (1 - n + n.regressors) * sum(diag(P %*% x)) + (n - n.regressors) * sum(diag(R %*% x)))
    # gout <- c(gout, sum(diag(iA %*% gnug)))
    #J <- c(0.5 * exp(theta[1:ncol.inputs]/2))# , exp(-theta[ncol.inputs + 1]) * nug^2) # if theta = 2 * log delta
    J <- c(0.5 * exp(theta[1:ncol.inputs]))#  if theta = log delta
    
    # attr(negloglik, "gradient") <- J * gout
    
    return(2 * J * gout)
}
