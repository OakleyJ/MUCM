#' @title Correlation functions
#' 
#' @description Functions to calculate the correlation between two points. 
#' The correlation functions give the correlation between the simulator output at any given vector of input values and 
#' the output at any other given input vector. Within this package, when these two vectors are the same then the correlation
#' function \code{C} at input vector \code{x} is \code{C(x,x) = 1}
#' In principal, the correlation function can take a very wide variety of forms. Only a few are included in this package, however, 
#' a user is able to define his own function. A user defined function must follow the same signature as those included in this package. 

#' @param inputs A data frame, matrix or vector containing the input values of the data where each row gives the observations. 
#' @param inputs2 An additional data frame, matrix or vector containing the input values of the data where each row gives the observations.
#' @param phi Estimate of the roughness parameter (a vector). 

#' @return Returns the correlation matrix \eqn{A} calculated by the specific function.

#' @details \code{corGaussian} uses the Gaussian correlation function defined as 
#' \deqn{\prod _{i=1} ^ {n} \exp [- \frac{(x_i - x_i')}{\delta _i} ^2]}{\prod _{i=1}^{n} exp(-(x_i - x_i') / \delta _i^2)}
#' @seealso See \code{\link{optim}} or  \code{\link{MCMCMetropolisGibbs}} where the function is used.  

#' @note These functions use \code{\link{rdist}} from package \pkg{fields} for speed.

#' @importFrom fields rdist

#' @author Sajni Malde and Jeremy Oakley
#' @export
corGaussian <- function(inputs, inputs2, phi) {
    
    if (missing(inputs2) || is.null(inputs2))
        return(corGaussianSquare(inputs, phi))
    
    delta <- exp(phi)
    exp(-(rdist(inputs / rep(delta, each = nrow(inputs)), inputs2 / rep(delta, each = nrow(inputs2))) ^ 2))
}

corGaussianSquare <- function(inputs, phi) {
    delta <- exp(phi)
    exp(-(rdist(inputs / rep(delta, each = nrow(inputs))) ^ 2))
}




############ ----------------   MATERN   ------------------- ################

#' @rdname corGaussian
#' @details \code{corMatern2.5} uses the Matern correlation function with \eqn{\nu = 2.5} defined as 
#' \deqn{(1 + 5^{0.5}r + \frac{5}{3}r^2)\exp(-5^{0.5}r)}{(1 + (5^0.5)*r + (5*r^2)/3) * exp(-(5^0.5)*r)} where \eqn{r} is distance between inputs \eqn{x_i} and \eqn{x_j}, scaled by
#' delta: \deqn{[\frac{(x_{i1}-x_{j1})^2}{\delta _1^2} +...+ \frac{(x_{id}-x_{jd})^2}{\delta_d^2}]^{0.5}}{[(x_{i1}-x_{j1})^2 / \delta_1^2 +...+ (x_{id}-x_{jd})^2 / \delta_d^2 ]^{0.5}}
#' @export
corMatern2.5 <- function(inputs, inputs2, phi) {
    
    if (missing(inputs2) || is.null(inputs2))
        return(corMaternSquare2.5(inputs, phi))
    
    delta <- exp(phi)
    
    # calculate distance matrix
    distances <- rdist(inputs / rep(delta, each = nrow(inputs)), inputs2/rep(delta, each = nrow(inputs2)))
    squared.distances <- distances ^ 2
    
    # return correlation matrix
    (1 + sqrt(5) * distances + (5 / 3) * squared.distances) * exp(-sqrt(5) * distances)
}

corMaternSquare2.5 <- function(inputs, phi) {
    n <- nrow(inputs)
    delta <- exp(phi)
    
    # calculate distance matrix
    distances <- rdist(inputs / rep(delta, each = n))
    squared.distances <- distances ^ 2
    
    # return correlation matrix
    (1 + sqrt(5) * distances + (5 / 3) * squared.distances) * exp(-sqrt(5) * distances)
}

