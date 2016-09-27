thetaSigmaParam <- function(theta, param, n.outputs) {
    
    if (param == "original") {
        ###original Parametrization inverse
        Sigma <- matrix(nrow = n.outputs, ncol = n.outputs)
        Sigma[upper.tri(Sigma, diag = TRUE)] <- theta
        Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
        
    } else if (param == "cholesky") {
        ###cholesky Parametrization inverse
        L <- matrix(0, nrow = n.outputs, ncol = n.outputs) 
        L[upper.tri(L, diag = TRUE)] <- theta 
        Sigma <- crossprod(L) # t(L) %*% L
        
    } else if (param == "log-cholesky") {
        ###log - cholesky Parametrization -inverse
        L <- matrix(0, nrow = n.outputs, ncol = n.outputs) 
        L[upper.tri(L, diag = TRUE)] <- theta 
        diag(L) <- exp(diag(L))
        Sigma <- crossprod(L) # t(L) %*% L
        
    } else if (param == "spherical") {
        ###spherical Parametrization -inverse
        
        l <- matrix(0,nrow = n.outputs, ncol = n.outputs)
        
        for (i in 1:n.outputs) {
            l[1,i] <- exp(theta[i])
        }
        for (i in 2:n.outputs) {
            for (j in 2:i) {
                exp1 <- exp(theta[n.outputs + (i - 2)*(i - 1)/2 + (j - 1)])
                l[j,i] <- (pi*exp1)/(1 + exp1)
            }
        }
        
        L <- matrix(0, nrow = n.outputs, ncol = n.outputs)
        for (i in 1:n.outputs) {
            for (j in 1:i) {
                L[j,i] <- l[1,i] 
                if (j < i) {
                    L[j,i] <- L[j,i]  * cos(l[j + 1,i])
                }
                if (j > 1) {
                    for (k in 2:j)
                        L[j,i] <- L[j,i] * sin(l[k,i])
                }
            }
        }
        Sigma <- crossprod(L) # t(L) %*% L
        
    } else if (param == "matrix-log") {
        ###Matrix logarithm Parametrization - inverse 
        log.SSigma2 <- matrix(0, nrow = n.outputs, ncol = n.outputs)
        log.SSigma2[upper.tri(log.SSigma2, diag = TRUE)] <- theta
        log.SSigma2[lower.tri(log.SSigma2)] <- t(log.SSigma2)[lower.tri(log.SSigma2)] 
        EG2 <- eigen(log.SSigma2, symmetric = TRUE)
        log.Lambda <- EG2$values
        Sigma <- EG2$vectors %*% diag(exp(log.Lambda)) %*% t(EG2$vectors)
    } else 
        stop("Select appropriate param!")
    Sigma
}






sigmaThetaParam <- function(Sigma, param) {
    if (param == "original") {
        ###original Parametrization inverse
        theta <- as.vector(Sigma[upper.tri(Sigma, diag = TRUE)])
        
    } else if (param == "cholesky") {
        ###cholesky Parametrization
        L <- chol(Sigma)
        theta <- L[upper.tri(L, diag = TRUE)]
        
        
    } else if (param == "log-cholesky") {
        ###log - cholesky Parametrization
        L <- chol(Sigma)
        diag(L) <- log(diag(L))
        theta <- L[upper.tri(L, diag = TRUE)]
        
        
    } else if (param == "spherical") {
        
        L <- chol(Sigma)
        n.outputs <- ncol(Sigma)
        l <- matrix(0, ncol = n.outputs, nrow = n.outputs)
        for (i in 1:n.outputs) {
            l[1,i] <- sqrt(sum((L[, i] ^ 2)))
            if (i > 1) {
                for (j in 2:i) {
                    l[j,i] <- acos(L[j - 1,i] / sqrt(sum(L[(j - 1):i,i] ^ 2))) 
                }
            }
        }
        
        theta <- vector(length = n.outputs*(n.outputs + 1)/2)
        for (i in 1:n.outputs) {
            theta[i] <- log(l[1,i])
        }
        for (i in 2:n.outputs) {
                for (j in 2:i)  {
                    theta[n.outputs + (i - 2)*(i - 1)/2 + (j - 1)] = log(l[j,i]/(pi - l[j,i]))
            }
        }
        
    } else if (param == "matrix-log") {
        ###Matrix logarithm Parametrization
        EG <- eigen(Sigma, symmetric = TRUE)
        Lambda <- EG$values
        log.SSigma <- EG$vectors %*% diag(log(Lambda)) %*% t(EG$vectors)
        theta <- log.SSigma[upper.tri(log.SSigma, diag = TRUE)]

    } else 
        stop("Select appropriate param!")
    
    theta
}
