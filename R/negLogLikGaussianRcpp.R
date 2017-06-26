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
#' 
#' 
#' @export
negLogLikGaussianRcpp <- function(theta, inputs, H, outputs, nugget = NULL, ...) {
    GetNllDelta(delta = theta, matX = inputs, matH = H, matY = outputs)
}