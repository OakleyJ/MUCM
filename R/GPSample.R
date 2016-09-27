#' @title GP Sample

#' @description Function to sample from posterior Gaussian process

#' @param n the number of samples required
#' @param mean a vector giving the means of the variables
#' @param var a positive-definite symmetric matrix specifying the covariance matrix of the variables

#' @author Jeremy Oakley.
#' 
#' @importFrom mvtnorm rmvnorm
#' @export
GPSample <- function(n, mean, var) {
    
  # Use pivoted Cholesky to decide which outputs should be sampled
  # Use the mean of the remaining outputs, conditional on the sampled values
    U <- chol(var, pivot = T)
    index.in <- attr(U, "pivot")[1:attr(U, "rank")]  # index of which points will be sampled
    index.out <- setdiff(1:length(mean), index.in)  # index of which points will not be sampled
    
    V1 <- var[index.in, index.in]
    m1 <- mean[index.in, ]
    m2 <- mean[index.out, ]
    
    new.y1 <- t(rmvnorm(n, m1, V1))
    
    # Calculate conditional mean of unsampled points, given the sampled points
    new.y2 <- matrix(t(m2), length(m2), n) + var[index.out, index.in] %*% solve(V1, new.y1 - matrix(t(m1), length(m1), n))
    
    f.sample <- matrix(0, length(mean), n)
    f.sample[index.in, ] <- new.y1
    f.sample[index.out, ] <- new.y2
    f.sample
} 
