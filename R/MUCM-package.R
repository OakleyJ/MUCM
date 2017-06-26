## MUCM-package.R : roxygen documentation for MUCM package

#' @title MUCM
#'
#' @description This package allows the user to estimate the output of a simulator at any point,
#' without actually running it, given a few hundred simulator runs (referred to as the training runs). The emulator essentially 
#' determines using Bayesian Analysis the posterior mean and variance of a Gaussian Process for a given input data set,
#' conditioned on the training runs and a user-specified prior mean function.  The emphasis is on complex codes that take 
#' weeks or months to run, and that have a large number of input parameters; many metrological prediction models fall into 
#' this class. A working example is given here for the main functions of this package, which should be the first point of reference.
#'
#' @name MUCM-package
#' @aliases MUCM
#' @docType package
#' 
#' @references 
#' J. Oakley 1999. 'Bayesian uncertainty analysis for complex computer codes', PhD thesis, University of Sheffield.
#' \cr Bastos, L. S. and O'Hagan, A. (2009). Diagnostics for gaussian process emulators, Technometrics, 51 (4): 425-438.
#' \cr Fricker, T. E., Oakley, J. E., & Urban, N. M. (2013). Multivariate Gaussian process emulators with nonseparable 
#'      covariance structures. Technometrics, 55(1), 47-56.
#' 
#' @author
#' Sajni Malde, Jeremy Oakley (\email{j.oakley@@sheffield.ac.uk}) and David Wyncoll.
#' 
#' @examples
#' # Plot the training data to look for trends
#' plot(surfebm[1:25, 1:3], lower.panel = NULL)
#' 
#' # Fit the emulator using a linear prior mean (default) and
#' # a Gaussian correlation function
#' fit <- fitEmulator(inputs = surfebm[1:25, 1:2], 
#'                    outputs = surfebm[1:25, 3, drop = FALSE],
#'                    cor.function = corGaussian)
#' 
#' # Use fitted emulator to predict posterior means and variances at the new points
#' predictions <- predict(fit, surfebm[26:35, 1:2], sd = FALSE, var.cov = TRUE)
#' 
#' # Compare predictions with true values for the new inputs
#' # Can also compare accuracy of prediction based on posterior variance
#' validateEmulator(fit, surfebm[26:35, 3], predictions, plot = TRUE)
#' @importFrom Rcpp evalCpp
#' @useDynLib MUCM
{
} 


