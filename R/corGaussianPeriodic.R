# correlation function is exp(-sin (pi/period * abs (x_i - x_j))^2 / delta^2) 

#' @rdname corGaussian
#' @param period A scalar or a vector indicating the period for each input parameter. If scalar, period for all input 
#' parameters will be set to the same scalar value as specified. Set period such that \eqn{cor(x_i, x_i + period) = 1}
#' @details \code{corGaussianPeriodic} uses the Guassian periodic correlation function defined as 
#' \deqn{\exp(- \delta ^2 sin(\frac{\pi}{period} |x_i - x_j|)^2)}{exp( - \delta ^2 * sin (\pi/period * |(x_i - x_j)|)^2)}
#' \code{corGaussianPeriodic2} uses a different way of accounting for 
#' the periodicity in the data. Both \code{corGaussianPeriodic} and \code{corGaussianPeriodic2} give the same results. 
#' @export
corGaussianPeriodic <- function(inputs, inputs2,  phi, period) {
    
    if (missing(inputs2) || is.null(inputs2))
        return(corGaussianPeriodicSquare(inputs, phi, period))
    
    n.inputs <- ncol(inputs)
    delta <- exp(phi)
    Ap <- 1
    
    if (length(period) == 1) {
        period <- c(rep(period, n.inputs))
    } else if (length(period) != n.inputs) {
        stop("length of period should equal the number of columns in inputs")
    }
    
    for (i in 1:n.inputs) {
        dist.mat <- as.matrix(rdist(inputs[, i], inputs2[, i]))
        Ap <- Ap * exp(-delta[i] ^ (-2) * sin(pi/period[i] * dist.mat) ^ 2)
    }
    Ap
}

corGaussianPeriodicSquare <- function(inputs, phi, period) {
    n.inputs <- ncol(inputs)
    delta <- exp(phi)
    Ap <- 1
    
    if (length(period) == 1) {
        period <- c(rep(period, n.inputs))
    } else if (length(period) != n.inputs) {
        stop("length of period should equal the number of columns in inputs")
    }
    
    for (i in 1:n.inputs) {
        dist.mat <- rdist(inputs[, i])
        Ap <- Ap * exp(-delta[i] ^ (-2) * sin(pi/period[i] * dist.mat) ^ 2)
        }
    Ap
}


#' @rdname corGaussian
#' @export
corGaussianPeriodic2 <- function(inputs, inputs2, phi, period) {
    
    if (missing(inputs2) || is.null(inputs2))
        return(corGaussianPeriodicSquare2(inputs, phi, period))
    
    n.inputs <- ncol(inputs)
    delta <- exp(phi)
    Ap2 <- 1
    
    for (i in 1:n.inputs) {
        dist.mat <- rdist(inputs[, i], inputs2[, i])
        stopifnot(dist.mat >= 0, dist.mat <= period)
        min.dist <- pmin(period - dist.mat, dist.mat)
        Ap2 <- Ap2 * exp(-(min.dist/delta[i]) ^ 2)
    }
    Ap2
}

corGaussianPeriodicSquare2 <- function(inputs, phi, period) {
    
    n.inputs <- ncol(inputs)
    delta <- exp(phi)
    Ap2 <- 1
    
    for (i in 1:n.inputs) {
        dist.mat <- rdist(inputs[, i])
        stopifnot(dist.mat >= 0, dist.mat <= period)
        min.dist <- pmin(period - dist.mat, dist.mat)
        Ap2 <- Ap2 * exp(-(min.dist/delta[i]) ^ 2)
    }
    Ap2
}

#' @rdname corGaussian
#' @param cor.funcs An ordered vector of length equal to number of inputs containing characters \code{'g'}, \code{'m2.5'}, \code{'gp'}, or \code{'gp2'}
#' as an indicator for whether a particular input should be calculated using 
#' \code{corGaussian}, \code{corMatern2.5}, \code{corGaussianPeriodic}, or \code{corGaussianPeriodic2} respectively. 
#' @param ... Additional arguments to be passed on to correlation functions.
#' @details \code{corCombined} is a function that allows the user to specify a different correlation function for each output 
#' (from the list of predefined functions). 
#' @export
corCombined <- function(inputs, inputs2, phi, cor.funcs, ...) {
    
    if (missing(inputs2) || is.null(inputs2))
        return(corCombinedSquare(inputs, phi, cor.funcs, ...))
    
    if (length(cor.funcs) != ncol(inputs)) 
        stop("length of cor.funcs should equal the number of columns in inputs")
    
    indexGaussian <- which(cor.funcs == "g")
    if (length(indexGaussian) > 0) {
        Ag <- corGaussian(inputs[, indexGaussian, drop = FALSE], inputs2[, indexGaussian, drop = FALSE], phi[indexGaussian])
    } else {
        Ag <- 1
    }
    
    indexMatern2.5 <- which(cor.funcs == "m2.5")
    if (length(indexMatern2.5) > 0) {
        Am2.5 <- corMatern2.5(inputs[, indexMatern2.5, drop = FALSE], inputs2[, indexMatern2.5, drop = FALSE], phi[indexMatern2.5])
    } else {
        Am2.5 <- 1
    }
    
    indexPeriodic <- which(cor.funcs == "gp")
    if (length(indexPeriodic) > 0) {
        Ap <- corGaussianPeriodic(inputs[, indexPeriodic, drop = FALSE], inputs2[, indexPeriodic, drop = FALSE], phi[indexPeriodic], ...)
    } else {
        Ap <- 1
    }
    
    indexPeriodic2 <- which(cor.funcs == "gp2")
    if (length(indexPeriodic2) > 0) {
        Ap2 <- corGaussianPeriodic2(inputs[, indexPeriodic2, drop = FALSE], inputs2[, indexPeriodic2, drop = FALSE], phi[indexPeriodic2], ...)
    } else {
        Ap2 <- 1
    }
    
    Ag * Am2.5 * Ap * Ap2
} 

corCombinedSquare <- function(inputs, phi, cor.funcs, ...) {
    
    if (length(cor.funcs) != ncol(inputs)) 
        stop("length of cor.funcs should equal the number of columns in inputs")
    
    indexGaussian <- which(cor.funcs == "g")
    if (length(indexGaussian) > 0) {
        Ag <- corGaussianSquare(inputs[, indexGaussian, drop = FALSE], phi[indexGaussian])
    } else {
        Ag <- 1
    }
    
    indexMatern2.5 <- which(cor.funcs == "m2.5")
    if (length(indexMatern2.5) > 0) {
        Am2.5 <- corMaternSquare2.5(inputs[, indexMatern2.5, drop = FALSE], phi[indexMatern2.5])
    } else {
        Am2.5 <- 1
    }
    
    indexPeriodic <- which(cor.funcs == "gp")
    if (length(indexPeriodic) > 0) {
        Ap <- corGaussianPeriodicSquare(inputs[, indexPeriodic, drop = FALSE], phi[indexPeriodic], ...)
    } else {
        Ap <- 1
    }
    
    indexPeriodic2 <- which(cor.funcs == "gp2")
    if (length(indexPeriodic2) > 0) {
        Ap2 <- corGaussianPeriodicSquare2(inputs[, indexPeriodic2, drop = FALSE], phi[indexPeriodic2], ...)
    } else {
        Ap2 <- 1
    }
    
    Ag * Am2.5 * Ap * Ap2
} 
