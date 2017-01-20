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
doA <- function(inputs){(rdist(inputs) ^ 2)}

# negloglik.lim <- 1e+06

#' @export
negLogLikGrad <- function(theta, inputs, H, outputs) {
    
    ncol.inputs <- ncol(inputs)
    n <- nrow(inputs)
    n.regressors <- ncol(H)
    
    # negative log likelihood of 2*log(roughness parameter), the
    # transformed roughness parameters
    
    #part of makeAdat
    Phi <- diag(1/exp(theta[1:ncol.inputs]/2))
    nug <- 1/(1 + exp(-theta[ncol.inputs + 1]))
    inputs.phi <- inputs %*% Phi 
    Dk <- lapply(1:ncol.inputs, function(k) doA(inputs.phi[, k, drop = FALSE]))  
    D <- matrix(rowSums(sapply(Dk, function(x) x)), n)
    A <- (1 - nug) * exp(-D)
    diag(A) <- 1
    
    # ensures A=c(inputs,inputs) positive definite
    L <- try(chol(A), silent = TRUE)
    if (class(L) == "try-error")
        return(negloglik.lim)
    
    # Calculate Q=H^T A^{-1}H via A=LL^T
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
    # multiple output ---> negloglik <- n.outputs * sum(log(diag(L))) + n.outputs * sum(log(diag(K))) + 0.5 * (n - n.regressors) * log(det.sigmasq.hat.prop)
    negloglik <- sum(log(diag(L))) + sum(log(diag(K))) + 0.5 * (n - n.regressors) * log(det.sigmasq.hat.prop)
    
    if (missing(negloglik) || (negloglik > negloglik.lim) || is.infinite(negloglik) || is.na(negloglik)) 
        return(negloglik.lim) 
    
    # calculating the gradient
    iA <- chol2inv(L)
    P <- iA - iA %*% H %*% solve(Q, crossprod(H, iA))
    P.outputs <- P %*% outputs
    R <- P - P.outputs %*% solve(crossprod(outputs, P.outputs), crossprod(outputs, P))
    # R2 <- P - tcrossprod(P.outputs)/(n - n.regressors - 2)/sigmahat.prop[1, 1]
    # R - R2
    
    gradA <- lapply(1:ncol.inputs, function(i) Phi[i, i] * Dk[[i]] * A)
    gradA[[ncol.inputs + 1]] <- -0.5 * exp(-D)
    diag(gradA[[ncol.inputs + 1]]) <- 0
    gout <- sapply(gradA, function(x) (1 - n + n.regressors) * sum(diag(P %*% x)) + (n - n.regressors) * sum(diag(R %*% x)))
    # gout <- c(gout, sum(diag(iA %*% gnug)))
    J <- c(0.5 * exp(theta[1:ncol.inputs]/2), exp(-theta[ncol.inputs + 1]) * nug^2)
    
    attr(negloglik, "gradient") <- J * gout
    
    return(negloglik)
}
