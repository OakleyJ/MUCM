negloglik.lim <- 10 ^ 6

#' @rdname negLogLik
#' @param phi Phi is the roughness parameter for each input and output
#' @param sigma Sigma is matrix of the covariance between different outputs.

#' @export
negLogLikLMC <- function(phi, sigma, inputs, H, outputs, cor.function, ...) {

    # phi and sigma need to be matrices
    V <- buildCovMat(phi, sigma, inputs = inputs, cor.function = cor.function, ...) 
    
    if (is.null(V)) # indirectly checking if V if Sigma is positive def
        return(negloglik.lim)
    
    # ensures V = c(inputs,inputs) positive definite
    L <- try(chol(V), silent = TRUE)
    if (class(L) == "try-error") 
        return(negloglik.lim)
    
    # Calculate Q=H^T V^{-1}H via V=LL^T
    w <- try(backsolve(L, H, transpose = TRUE))  # = solve(tL, H, tol = 1e-200)
    if (class(w) == "try-error") 
        return(negloglik.lim)
    Q <- crossprod(w)                            # = t(w) %*% w
    
    # Calculate V^{-1}y via V=LL^T
    mat1 <- backsolve(L, outputs, transpose = TRUE) 
    mat2 <- backsolve(L, mat1)                   # = solve(L, mat1)
    
    # Check that Q is positive definite and Calculate Betahat
    K <- try(chol(Q)) 
    if (class(K) == "try-error") 
        return(negloglik.lim)
    mat3 <- backsolve(K, t(H), transpose = TRUE) 
    mat4 <- mat3 %*% mat2
    betahat <- backsolve(K, mat4)
    Hbeta <- H %*% betahat
    out.min.Hbeta <- outputs - Hbeta
    
    # calculate (y-hbeta)^T (V)^{1/2} (y -hbeta)
    w2 <- backsolve(L, out.min.Hbeta, transpose = TRUE) 
    sigmasq.hat.prop <- crossprod(w2) 
    negloglik <- sum(log(diag(L))) + sum(log(diag(K))) + 0.5 * sigmasq.hat.prop
    
    if (missing(negloglik) || (negloglik > negloglik.lim) || is.infinite(negloglik) || is.na(negloglik)) 
        return(negloglik.lim)
    
    negloglik
    
}


