#' @title Fit Emulator

#' @description A function that fits an emulator to the data given resulting in an object of class \code{emulatorFit}.

#' @param inputs  A data frame, matrix or vector containing the input values of the training data.
#'        For matrices and data frames, use one column per input variable.  
#' @param outputs  A data frame, matrix or vector containing the output values of the training data. 
#'        For multivariate output, use one column per output variable.
#' @param prior.mean  Either an object of class \code{'formula'} or \code{'Formula'} that specifies the regressors in the
#'        prior mean function, or one of the two strings \code{'linear'} or \code{'constant'} that specify
#'        a linear prior mean in each input and a constant prior mean respectively.
#'        See details if \code{method = 'LMC'} and you want to specify a formula.
#' @param cor.function Specifies a correlation function used as part of the prior information for the emulator. 
#'        This package has options of: \code{\link{corGaussian}}, \code{\link{corMatern2.5}}, \code{\link{corGaussianPeriodic}},
#'         \code{\link{corGaussianPeriodic2}} and \code{\link{corCombined}}. One can also specify a user defined function, however this must follow the signature of, for example \code{\link{corGaussian}}.
#' @param phi.opt  An optional vector that indicates the optimal phi values used to fit the emulator. It avoids running the optimiser to generate optimal phi values (in this case, the phi.init values are ignored).
#'        If not given phi.init values are optimised.
#' @param sigmasq.opt  An optional scalar that indicates the optimal sigmasq value used to fit the emulator. It avoids running the optimiser to generate optimal sigmasq value (in this case, the sigmasq.init values are ignored).
#'        If not given sigmasq.init values are optimised.
#' @param MCMC.iterations An optional value that indicates how many iterations to run the MCMC for. 
#'        A non-positive value indicates not to apply the MCMC to estimate the starting values for the optimisation. 
#'        In this case it will use the values given by phi.init as the starting values for the optimisation.
#'        Note: MCMC is run in order to get better starting values for the optim function that evaluates an optimum value for the hyperparameters,
#' @param phi.init  An optional vector used as starting values for either the optimisation or the MCMC output depending on \code{MCMC.iterations} given. 
#'        If provided, it should be either a vector of length 1, where it will be assumed that all inputs have the same phi.init values, or a vector of length tequal to the number of inputs. The default value 
#'        is a vector of zeros, however when \code{method = 'LMC'}, the default value is derived by fitting individual optimized estimates from the univariate emulator. 
#'        Phi refers to the roughness parameter used for the calcuations.
#' @param sigmasq.init A optional scalar used as a starting value for either the optimisation or the MCMC output depending on \code{MCMC.iterations} given. 
#'        If provided, it should be a vector of length the same as the number of inputs. The default value is a vector of zeros, however when \code{method = 'LMC'}, 
#'        the default value is derived by fitting a multivariate separate emulator and obtaining the estimates of the optimised \code{sigmasq.hat} values. For \code{method = 'separable'} sigmsq refers to an unknown 
#'        scale parameter to calculate the covariance matirx. For \code{method = 'LMC'}, sigmasq refers to the between outputs covaraince matrix.
#' @param MC.plot  If \code{TRUE}, produces a trace plot of the MCMC output of log likelihood against the number of iterations. (default: \code{FALSE})
#' @param nugget For noisy data, a scalar or a vector giving the observation variance for each training data point. 
#'        If scalar, nugget for each training data point will be set to the same scalar value as specified.
#'        Currently only available with univariate output.
#' @param method Can take string \code{'seperable'} (default) or \code{'LMC'}. The seperable method is used for univariate output. For multivariate output, this distinction is more important. 
#'        The separable emulator restricts the model to have just one set of hyper-parameters that are shared across all outputs, however the LMC (Linear model of coregionalisation) emulator
#'        allows for flexibilty of individual hyper-parameters for each output. See the references for further information.
#' @param optim.method  The optimisation method to use for estimation of the unknown hyper-parameters. See \code{\link{optim}}.
#' @param optimise.interval  For \code{optim.method = "Brent"} with only one parameter to be estimated,
#'        a vector passed on to \code{\link{optimise}} containing the interval to search within.
#' @param ... Additional arguments to be used for \code{fitEmulator} and passed on to correlation functions.

#' @return \code{fitEmulator} returns an object of class 'emulatorFit'. 
#'  An object of class 'emulatorFit' is a list containing at least the following components:
#' \tabular{ll}{
#'  \code{betahat} \tab Posterior mean of \eqn{\beta}.  \cr
#'  \code{sigmasq.hat} \tab Posterior mean of \eqn{\sigma^2}{sigmasq}.  \cr
#'  \code{out.min.Hbeta} \tab \eqn{(y - H \beta)}. \cr
#'  \code{phi.hat} \tab Optimum maximum likelihood estimate of the correlation length parameter  \cr 
#'  \code{training.inputs} \tab Inputs used to train the emulator.  \cr
#'  \code{training.outputs} \tab Outputs used to train the emulator.  \cr
#'  \code{H.training} \tab Part of the \code{prior.mean}.  \cr
#'  \code{tL} \tab The cholesky decomposition of matrix \eqn{A}.  \cr
#'  \code{K} \tab The cholesky decomposition of the inverse of \code{tL} multiplied with \code{H}.  \cr
#'  \code{n.train} \tab Number of training runs (rows of \code{D}).  \cr
#'  \code{n.regressors} \tab Number of columns of \code{H}.  \cr
#'  \code{n.outputs} \tab Number of outputs.  \cr
#'  \code{cor.function} \tab The correlation function chosen for the prior variance. \cr 
#'  \code{cor.functions.args} \tab Extra arguments used for \code{cor.functions} and passed on for further evaluation.  \cr
#'  \code{formula} \tab The regressor terms in the prior mean function \code{data}.  \cr
#'  \code{log.lik} \tab  It gives the log-likelihood value.  \cr
#'  \code{opt.convergence} \tab Only given if the optimisation was run. It is An integer code where 0 indicates
#'       successful completion. See \code{\link{optim}} for more information.  \cr

#'  \code{opt.message} \tab Only given if the optimisation was run.
#'       A character string giving any additional information returned by the optimizer, or \code{NULL}. 
#' }
#' 
#' 
#' If \code{method = 'separable'} the following components are included.
#' \tabular{ll}{
#'  \code{Ainv.H} \tab \eqn{A^{-1} H}{Ainv * H}. \cr
#'  \code{Ainv.e} \tab \eqn{A^{-1} (y - H \beta)}{Ainv * out.min.Hbeta}.  \cr
#'  \code{nugget} \tab A vector giving the observation variance for each training data point.
#' }
#' 
#' 
#' If \code{method = 'LMC'}, \code{fitEmulator} returns an object of class \code{'emulatorFit'} and \code{'emulatorFitLMC'}. The following components are also included.
#' \tabular{ll}{
#'  \code{Vinv.H} \tab \eqn{V^{-1} H}{Vinv * H}. \cr
#'  \code{Vinv.e} \tab \eqn{V^{-1} (y - H \beta)}{Vinv * out.min.Hbeta}.  \cr
#'  \code{sigmasq.param} \tab The type of parametrization chosen to appopriately parametrize the covariance matrix (sigmasq)
#' }


#' @seealso
#' If all data is in single matrix, see \code{\link{fitEmulatorData}} for a more conveniant method. 
#' To predict at new points, see \code{\link{predict.emulatorFit}}
#'  
#' @details This function implements emulator fitting based on the procedures described in
#' \url{http://mucm.aston.ac.uk/MUCM/MUCMToolkit}.
#' 
#' For the \code{prior.mean} arguement, if \code{method = 'separable'} then the left hand side of the equation is not required, 
#' but if specified, should include variables that are colnames of \code{outputs}.
#' However, if \code{method = 'LMC'} and a formula is specified, it should have a similar format to the following.
#' Given 3 output variables \code{(y1, y2} and \code{y3)} and 3 input variables \code{(x1, x2} and \code{x3):  y1 | y2 | y3 ~ x1 + x2 | x1 + x3 | x2 + x3}. 
#' (The actual equations can vary). Note: the vertical lines separate prior mean formula for 
#' each output, for instance in the example above the prior mean formula for \code{y1} is given as (\code{y1 = x1 + x2}), and for \code{y3} is given as (\code{y3 = x2 + x3}).
#' This format allows the user to specify different prior mean formulas for each output. 

#' @references Fricker, T. E., Oakley, J. E., & Urban, N. M. (2013). Multivariate Gaussian process emulators with nonseparable 
#'             covariance structures. Technometrics, 55(1), 47-56.
#' 
#' @examples
#' fit1 <- fitEmulator(inputs = surfebm[1:25,1:2], outputs = surfebm[1:25,3]) 
#' 
#' fit2 <- fitEmulator(inputs = surfebm[1:25,1:2], outputs = surfebm[1:25,3],
#'                     prior.mean = 'constant', cor.function = corMatern2.5,
#'                     MCMC.iterations = 1000)

#' @export
fitEmulator <- function(inputs, outputs, prior.mean = "linear", 
                        cor.function = corGaussian, method = c("separable", "LMC"), ...) {
    
    method <- match.arg(method)
    
    if (method == "LMC")
        fit <- fitEmulatorLMC(inputs, outputs, prior.mean, cor.function, ...)
    else 
        fit <- fitEmulatorSEP(inputs, outputs, prior.mean, cor.function, ...)
    
    
    fit$call <- match.call()
    
    fit
}

