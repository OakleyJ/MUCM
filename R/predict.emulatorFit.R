#' @title Prediction using Emulators

#' @description Predicts value and confidence interval at new inputs using Gaussian Process Emulation. 
#' This function should be preceded by the \code{\link{fitEmulator}} function.

#' @param object  A fit object of class inheriting from \code{'emulatorFit'}.
#' @param newdata  A data matrix of input(s) at which emulation is desired (new inputs). 
#'        Must contain at least all parameters given in \code{object$training.inputs}.
#'        If missing, the fitted inputs \code{object$training.inputs} are used.
#' @param var.cov  Optionally calculates posterior variance covariance matrix. Default is set to FALSE. For large numbers of training and prediction data, this is quite time consuming.
#' @param sd  Optionally calculates only the posterior standard deviation. Default is set to \code{TRUE}.
#' @param tol The tolerance for capping negative small values of posterior standard deviation to zero.
#'        The default is \eqn{-10^{-11}}{-10^-11}.
#' @param ...  Further arguments not used and an error is thrown if provided.
#' @return The function returns a list containting the following components: 
#' \tabular{ll}{
#'     \code{posterior.mean}  \tab    Approximation of the outputs for the given inputs in \code{newdata}   \cr
#'     \code{posterior.variance}  \tab Variance covariance matrix around this approximation \cr
#'     \code{standard.deviation} \tab  Standard Deviation of the approximation. It equals the square-root
#'                                    of the diagonal of the \code{posterior.variance}
#' }
#' When the number of outputs to emulate is more than 1, \code{method = 'separable'}, and object is of class \code{"emulatorFit"} 
#' two extra values are returned from this function. These are 
#' \tabular{ll}{
#'     \code{correlation.Matrix}  \tab A spatial correlation matrix.  \cr
#'     \code{sigmahat}  \tab    A between outputs covariance matrix.  \cr
#' }
#' @details Note that when using the LMC method, calculating the posterior variance is quite time-consuming. 
#' @references Oakley, J. (1999). Bayesian uncertainty analysis for complex computer codes, Ph.D. thesis, University of Sheffield.
#' @author Originally written by Jeremy Oakley. Modified by Sajni Malde.
#' @export
predict.emulatorFit <- function(object, newdata, var.cov = FALSE, sd = TRUE, tol = -1e-11, ...) {
    # check no further arguments are provided
    l <- list(...)
    if (length(l) > 0) 
        stop("unused argument (", names(l[1]), " = ", l[[1]], ")")
    
    # ensure new and old fits are backward compatible
    object <- updateFit(object)
    
    # set newdata to fitted data if missing
    if (missing(newdata) || is.null(newdata))
        newdata <- object$training.inputs
    
    # checking if data is a dataframe
    newdata <- as.data.frame(newdata)
    
    n.train <- object$n.train
    n.pred <- nrow(newdata)
    n.all <- n.train + n.pred
    
    # extract inputs from newdata
    input.names <- colnames(object$training.inputs)
    # check if all inputs are given
    if (!all(input.names %in% colnames(newdata))) 
        stop("one or more input variables are missing from the newdata matrix")
    
    new.inputs <- as.matrix(newdata[, input.names, drop = FALSE])
    
    H.pred <- model.matrix(object$formula, data = as.data.frame(new.inputs), lhs = 0, drop = FALSE)
    
    cov.training.prediction <- do.call(object$cor.function, c(list(new.inputs, object$training.inputs, object$phi), object$cor.function.args))
    
    post.mean <- H.pred %*% object$betahat + cov.training.prediction %*% object$Ainv.e
    predict <- list(posterior.mean = post.mean)
    
    if (var.cov || sd) {
        temp.constant <- H.pred - cov.training.prediction %*% object$Ainv.H
        w3 <- backsolve(object$L, t(cov.training.prediction), transpose = TRUE) # solve(object$tL, t(cov.training.prediction))
        w4 <- backsolve(object$K, t(temp.constant), transpose = TRUE) # solve(object$K, t(temp.constant))
        
        if (var.cov) {
            cov.pred <- do.call(object$cor.function, c(list(new.inputs, phi = object$phi), object$cor.function.args))
            colnames(cov.pred) <- rownames(new.inputs)
            rownames(cov.pred) <- rownames(new.inputs)
            correlation.matrix <- (cov.pred - crossprod(w3) + crossprod(w4)) #(t(w3) %*% w3) + (t(w4) %*% w4))
            predict$posterior.variance <- kronecker(object$sigmasq, correlation.matrix, make.dimnames = TRUE)
            
            if (object$n.outputs > 1) {
                predict$correlation.matrix <- correlation.matrix
                predict$sigmasq <- object$sigmasq
            }
            if ((min(diag(predict$posterior.variance), na.rm = TRUE) < 0) && !sd)  {
                warning(paste("Negative minimum posterior variance:", min(diag(predict$posterior.variance), na.rm = TRUE)))
            }
        }
        
        if (sd) {
            if (var.cov == TRUE) {
                var <-  diag(predict$posterior.variance)
            } else {          
                diag.corr.mat <- 1 - colSums(w3 ^ 2) + colSums(w4 ^ 2)
                if (is.null(object$nugget))
                    rep.sigmasq <- rep(diag(object$sigmasq), each = n.pred)
                else
                    rep.sigmasq <- rep(object$sigmasq, each = n.pred)
                
                rep.diag.corr.mat <- rep(diag.corr.mat, times = object$n.outputs)
                var <- rep.sigmasq * rep.diag.corr.mat
            }
            
            var[(var < 0) & (var > tol)] <- 0
            predict$standard.deviation <- sqrt(var)
            
            #unstack standard.deviation and add names
            predict$standard.deviation <- matrix(predict$standard.deviation, nrow = n.pred, ncol = object$n.outputs)
            colnames(predict$standard.deviation) <- colnames(object$sigmasq.hat)
            rownames(predict$standard.deviation) <- rownames(new.inputs)
            
            if ((!var.cov) && (min(var, na.rm = TRUE) < 0))
                warning("missing posterior standard deviation values as some variances are negetive")
            else if (var.cov && (min(var, na.rm = TRUE) < 0))
                warning(paste("Negative minimum posterior variance:", min(diag(predict$posterior.variance), na.rm = TRUE),
                              "and missing some posterior standard deviation values"))
            else if (var.cov && (min(diag(predict$posterior.variance), na.rm = TRUE) < 0))
                warning(paste("Negative minimum posterior variance:", min(diag(predict$posterior.variance), na.rm = TRUE)))
        }
    }
    
    # return
    predict
}
