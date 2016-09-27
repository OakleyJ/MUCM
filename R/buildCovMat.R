#' @title Building a Covariance Matrix for the LMC method

#' @description This function is used to build a covariance matrix for multivariate outputs. 

#' @param sigma Estimate of the between outputs variance covariance matrix.
#' @param phi Estimate of the roughness parameter for each input and output.
#' @param inputs A data frame, matrix or vector containing the input values of the data.
#' @param inputs2 A data frame, matrix or vector containing additional input values. (Used when predicting the model at new points)
#' @param cor.function Specifies a correlation function used as part of the prior information for the emulator.
#' @param ... Additional arguments to be passed on to \code{cor.function}.

#' @return The function returns a covariance matrix, eqn{V}.

#' @author Sajni Malde
#' @export

buildCovMat <- function(phi, sigma, inputs, inputs2, cor.function, ...) { 

    n.outputs <- ncol(sigma)
    
    #eigendecomposition of Sigma ==> R
    EG <- eigen(sigma, symmetric = TRUE)
    
    if (any(EG$values < 0))
        return(NULL)
    
    R <- crossprod(t(EG$vectors) * sqrt(sqrt(EG$values))) # EG$vectors %*% sqrt(diag(EG$values)) %*% t(EG$vectors) 
    
    # initialise V - variance matrix
    V <- 0
    
    # looping through outputs
    for (i in seq_len(n.outputs)) {
        A <- cor.function(inputs,inputs2, phi = phi[ , i], ...)
        big.sigma <- tcrossprod(R[, i, drop = FALSE]) # R[, i, drop = FALSE] %*% t(R[, i, drop = FALSE])
        V <- V + kronecker(big.sigma, A, make.dimnames = TRUE)
    }
    V
}    
    