#' @title  Make Emulator Constants 
#' @description  Function for collating and calculating constants needed for fitting the emulator
#' @param phi Optimum Maximum likelihood estimate of the correlation length parameter
#' @param inputs  A data frame, matrix or vector containing the input values of the training data.
#' @param outputs A data frame, matrix or vector containing the output values of the training data.In \code{makeEmulatorConstantLMC} the outputs should be stacked (either a vector or a matrix with 1 column). 
#' @param H  A matrix of prior mean regressors for the training data
#' @param cor.function Specifies a correlation function used as part of the prior information for the emulator. 
#'        This package has options of: \code{\link{corGaussian}}, \code{\link{corMatern2.5}}, \code{\link{corGaussianPeriodic}},
#'         \code{\link{corGaussianPeriodic2}} and \code{\link{corCombined}}. One can also specify a user defined function. 
#' @param ... additional arguments to be passed on to correlation functions (see \code{\link{corGaussian}})
#' @param sigmasq.hat Optimum Maximum likelihood estimate of the unknown scale parameter.
#' @param nugget For noisy data, a vector giving the observation variance for each training data point. 
#' @return The function returns a list containting some components as the \code{emulatorFit} class. 
#'          See \code{\link{fitEmulator}} for details 
#' @author Originally written by Jeremy Oakley. Modified by Sajni Malde
#' @references J. Oakley 1999. 'Bayesian uncertainty analysis for complex computer codes', PhD thesis, University of Sheffield.
makeEmulatorConstants <- function(phi, inputs, H, outputs, cor.function, nugget = NULL, sigmasq.hat = NULL, ...) {
    
    n <- nrow(inputs)
    n.regressors <- ncol(H)
    n.outputs <- ncol(outputs)
    A <- cor.function(inputs, phi = phi, ...)
    
    colnames(A) <- rownames(inputs)
    rownames(A) <- rownames(inputs)
    
    if (!is.null(sigmasq.hat))
        A <- A + diag(nugget)/sigmasq.hat
    
    L <- chol(A)
    w <- backsolve(L, H, transpose = TRUE) # solve(tL, H, tol = 1e-100)
    Q <- crossprod(w) #t(w) %*% w
    mat1 <- backsolve(L, outputs, transpose = TRUE) # solve(tL, outputs)
    mat2 <- backsolve(L, mat1) # solve(L, mat1)
    K <- chol(Q) # K <- t(chol(Q))
    mat3 <- backsolve(K, t(H), transpose = TRUE) # solve(K, t(H))                                          
    mat4 <- mat3 %*% mat2
    betahat <- backsolve(K, mat4) # solve(t(K), mat4)
    Hbeta <- H %*% betahat
    out.min.Hbeta <- outputs - Hbeta
    w2 <- backsolve(L, out.min.Hbeta, transpose = TRUE) # solve(tL, out.min.Hbeta)
    
    if (is.null(sigmasq.hat)) {
        sigmasq.hat <- ((n - n.regressors - n.outputs - 1) ^ (-1)) * crossprod(w2) #t(w2) %*% w2
    }
    Ainv.H <- solve(A, H, tol = 1e-100)
    Ainv.e <- solve(A, out.min.Hbeta, tol = 1e-100)
    
    dimnames(betahat) <- list(colnames(H), colnames(outputs))
    dimnames(K) <- list(colnames(H), colnames(H))
    if (!is.null(sigmasq.hat) && !is.null(dim(sigmasq.hat)))
        dimnames(sigmasq.hat) <- list(colnames(outputs), colnames(outputs))
    
    fit <- list(betahat = betahat, sigmasq.hat = sigmasq.hat, out.min.Hbeta = out.min.Hbeta, Ainv.H = Ainv.H, n.outputs = n.outputs,
        training.outputs = outputs, H.training = H, L = L, K = K, Ainv.e = Ainv.e,
        cor.function = cor.function, cor.function.args = list(...), 
        n.regressors = n.regressors, nugget = nugget, A = A)
    
    fit
} 


############### --------- make Emulator Constants LMC ---------- #############


#' @rdname makeEmulatorConstants
#' @param sigma Optimum Maximum likelihood estimate of the between outputs variance covariance matrix.
makeEmulatorConstantsLMC <- function(phi, sigma, inputs, H, outputs, cor.function, ...) {
    
    V <- buildCovMat(phi = phi, sigma = sigma, inputs = inputs, cor.function = cor.function, ...)
    dimnames(V) <- list(rownames(outputs), rownames(outputs))
    
    L <- chol(V) # tL <- t(L)
    w <- backsolve(L, H, transpose = TRUE) # solve(tL, H, tol = 1e-100)
    Q <- crossprod(w) # t(w) %*% w
    mat1 <- backsolve(L, outputs, transpose = TRUE) # solve(tL, outputs) 
    mat2 <- backsolve(L, mat1) # solve(L, mat1) 
    K <- chol(Q) # K <- t(chol(Q))
    mat3 <- backsolve(K, t(H), transpose = TRUE) # solve(K, t(H))                                          
    mat4 <- mat3 %*% mat2
    betahat <- backsolve(K, mat4) # solve(t(K), mat4)
    Hbeta <- H %*% betahat
    out.min.Hbeta <- outputs - Hbeta
    Vinv.H <- solve(V, H, tol = 1e-100)
    Vinv.e <- solve(V, out.min.Hbeta, tol = 1e-100)
    
    dimnames(betahat) <- list(colnames(H), 1)
    dimnames(K) <- list(colnames(H), colnames(H))
    
    fit <- list(betahat = betahat, sigmasq.hat = sigma, out.min.Hbeta = out.min.Hbeta, Vinv.H = Vinv.H,
                phi.hat = phi, training.inputs = inputs, training.outputs = outputs, H.training = H, L = L, K = K, Vinv.e = Vinv.e, 
                cor.function = cor.function, cor.function.args = list(...), 
                V = V)
    
    fit
} 
