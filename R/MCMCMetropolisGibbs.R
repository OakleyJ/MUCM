#' @title MCMC using Metroplis Hastings within Gibbs

#' @description Generate a Markov Chain of the parameters in the correlation function Using Metropolis-Hastings within Gibbs 

#' @param inputs A data frame, matrix or vector containing the input values of the training data.
#' @param outputs A data frame, matrix or vector containing the output values of the training data. 
#' @param fn A function used to maximise the negetive log likelihood
#' @param H A matrix of prior mean regressors from the training data
#' @param MCMC.iterations The number of iterations that MCMC should be run for
#' @param starting.values the starting values for which the MCMC can start running
#' @param proposal.sd is the standard deviation of the random walk proposal (default \code{0.1})
#' @param cor.function Specifies a correlation function used as part of the prior information for the emulator. 
#'        This package has options of: \code{\link{corGaussian}}, \code{\link{corMatern2.5}}, \code{\link{corGaussianPeriodic}},
#'         \code{\link{corGaussianPeriodic2}} and \code{\link{corCombined}}. One can also specify a user defined function. 
#' @param MC.plot If \code{TRUE}, produces a trace plot of the MCMC output of log likelihood against the number of iterations. (default=\code{TRUE})
#' @param ... additional arguments to be passed on to correlation functions (see \code{\link{corGaussian}})

#' @return The function returns a list containting the following components: 
#' \tabular{ll}{
#' \code{density.sample} \tab  The negetive log likelihood of the MCMC output at the starting value and at each iteration \cr
#' \code{theta.sample} \tab  A matrix of the theta sample at the starting value and at each iteration   \cr }
#' @note Note that this function first calculates the negetive log likelihood of the starting values so returns \code{MCMC.iterations} + 1 values. 
#' @author Originally written by Jeremy Oakley. Modified by Sajni Malde
## CONCERN IF IT IS TITLED APPROPRIATELY ##
#' @export
MCMCMetropolisGibbs <- function(inputs, outputs, fn, H, MCMC.iterations, starting.values, proposal.sd = 0.1, cor.function, MC.plot = TRUE, ...) {
    
    
    # H<-make.prior.mean.regressors(inputs) m <- ncol(H)
    theta.sample <- matrix(0, MCMC.iterations + 1, length(starting.values))
    n.theta <- ncol(theta.sample)
    density.sample <- rep(0, MCMC.iterations + 1)
    theta.sample[1, ] <- starting.values
    density.sample[1] <- -fn(starting.values, inputs = inputs, H = H, outputs = outputs, cor.function = cor.function, ...)
    current.density <- density.sample[1]
    
    if (MCMC.iterations >= 1) {
        for (i in 1:MCMC.iterations) {
            current.theta <- theta.sample[i, ]
            for (j in 1:n.theta) {
                proposal <- current.theta
                proposal[j] <- rnorm(1, current.theta[j], proposal.sd)
                new.density <- -fn(proposal, inputs = inputs, H = H, outputs = outputs, cor.function = cor.function, ...)
                
                u <- runif(1)
                acceptance.prob <- exp(new.density - current.density)
                if (u < acceptance.prob) {
                    current.theta <- proposal
                    current.density <- new.density
                }
            }
            theta.sample[i + 1, ] <- current.theta
            density.sample[i + 1] <- current.density
        }
        
        if (MC.plot == TRUE) 
            try(plot(0:MCMC.iterations, density.sample, type = "l", xlab = "iteration", ylab = "log likelihood"))
    }
    list(density.sample = density.sample, theta.sample = theta.sample)
} 
