#' @title Transform output from Prcomp / Princomp

#' @description  These functions reverse the outcome of the \code{predict} function 
#' when the class of the first argument is either \code{prcomp} or \code{princomp}. 
#' These functions are used to transform data where the variables are principal components of a dataset 
#' of original variables, given the \code{object}. 

#' @param object    Object of class inheriting from \code{princomp} or \code{prcomp}
#' @param newdata   A matrix of principal components which is to be converted back to the original output variables. 
#' An important feature here is that \code{newdata} can have a reduced number of dimensions with which to predict, however it has to be a matrix.
#' @return \code{backPredict} returns a matrix of data in terms of the original output variables.

#' @details  \code{backPredict} defines a function to reverse the outcome of the \code{predict.prcomp} or \code{predict.princomp} function.  
#' It is used to transform data (typically values predicted using regressions) where the variables are principal components of a dataset of original variables, given in the \code{object}. 

#' @seealso  \code{\link{princomp}}, \code{\link{prcomp}}, \code{\link{predict.prcomp}}

#' @author Sajni Malde

#' @examples 
#' library(mvtnorm)
#' 
#' #defining data
#' Sigma <- matrix(c(10.3, 3.6,3.6,2.4),2,2) 
#' out <- rmvnorm(n=100000, c(1.2,2.4), Sigma); colnames(out) <- c('Y1', 'Y2')
#' cov(out)
#' new.out <- rmvnorm(n=100000, c(1.2,2.4), Sigma); colnames(new.out) <- c('Y1', 'Y2')
#' 
#' # using the 'princomp' function to convert the data
#' out.pca <- princomp(out) # converting the data to principal components
#' # predicting principal components for a new data set
#' new.out.pca<- predict(out.pca, new.out) 
#' 
#' # applying backPredict and comparing values
#' backpredict <- backPredict(out.pca, new.out.pca)
#' head(backpredict) # should be quite similar to head(new.out)
#' # reducing the dimensions # only 1 components used to 'backPredict'
#' backpredict2 <- backPredict(out.pca, new.out.pca[, 1, drop = FALSE])
#' 
#' # applying backPredictVar and comparing values
#' # check to see if assuming off diagonals are equal to zero is a reasonable assumption
#' cov(new.out.pca) 
#' (backpredict.var <- backPredictVar(out.pca, new.vardata = matrix(diag(cov(new.out.pca)), 
#'                                    ncol = 2), only.var = FALSE)[1,,])
#' # backpredict.var should be very close to cov(out)
#' @export 
backPredict <- function(object, newdata) {
    
    # define n to be the number of components chosen illustrated by ncol(newdata)
    if (!is.matrix(newdata)) 
        stop(" \"newdata\" has to be a matrix")
    
    n <- ncol(newdata)
    
    # reverse PCA renaming object$loadings and object$rotations
    if (class(object) == "princomp") {
        loadings <- object$loadings[, 1:n, drop = FALSE]
    } else if (class(object) == "prcomp") {
        loadings <- object$rotation[, 1:n, drop = FALSE]
    } else {
        stop("object must be either \"princomp\" or \"prcomp\"")
    }
    
    # scaled reverse PCA
    ret <- loadings %*% t(newdata)
    
    # reverse the scaling appropriately if it was used
    if (class(object) == "princomp") {
        ret <- object$center + ret * object$scale
	} else if (class(object) == "prcomp") {
        if (!identical(object$scale, FALSE)) 
            ret <- ret * object$scale
        if (!identical(object$center, FALSE)) 
            ret <- object$center + ret
    }
    
    # return the transposed dataset
    ret <- t(ret)
    ret
}

#' @rdname backPredict
#' @param new.vardata A matrix of variances where each column represents the variances of one principal component's prediction for many observations, 
#' and each row represents the variances of multiple principal components of one observation. 
#' @param only.var This parameter lets the user choose how the output of the function \code{backPredictVar} should appear. 
#' The default option is set to \code{TRUE}. In this case the output of the function only contains the variances of each observation with the output variable. 
#' If \code{FALSE} this function returns an array of 3 dimensions representing the cross covariances of output variables and observations. 
#' The 1st dimension is the number of rows in \code{new.vardata} and the 2nd and 3rd dimensions represent a variance covariance matrix 
#' between output \eqn{j} and output \eqn{k} for row \eqn{i} in \code{new.vardata}. 

#' @return \code{backPredictVar} returns a matrix of variances (and covariances if \code{only.var = FALSE}) of each observation with the original output variable

#' @details  \code{backPredictVar} defines a similar function to reverse the outcome of the \code{predict.prcomp} or \code{predict.princomp} function.  
#' It is used to transform data that represents variances of principal components of an original dataset given in the \code{object}. 
#' It uses the linear variance transformation property. 
#' @export
backPredictVar <- function(object, new.vardata, only.var = TRUE) {
    
    # new.vardata is a n x p matrix. however we need to transform the data one datapoint at a time
    
    # define n to be the number of components chosen illustrated by ncol(new.vardata)
    if (!is.matrix(new.vardata)) 
        stop(" \"new.vardata\" has to be a matrix")
    
    n <- ncol(new.vardata)
    
    # renaming object$loadings and object$rotations as 'loadings' to make it for 'prcomp' and 'princomp'
    if (class(object) == "princomp") {
        loadings <- object$loadings[, 1:n, drop = F]
    } else if (class(object) == "prcomp") {
        loadings <- object$rotation[, 1:n, drop = F]
    } else {
        stop("object must be either \"princomp\" or \"prcomp\"")
    }
    
    # Defining the output dimensions depending on only.var
    if (!only.var) {
        PCA.var.data <- array(dim = c(nrow(new.vardata), nrow(loadings), nrow(loadings)))
        dimnames(PCA.var.data) <- list(rownames(new.vardata), rownames(loadings), rownames(loadings))
    } else {
        PCA.var.data <- matrix(nrow = nrow(new.vardata), ncol = nrow(loadings))
        colnames(PCA.var.data) <- rownames(loadings)
        rownames(PCA.var.data) <- rownames(new.vardata)
    }
    
    # for loop to loop through all the data points
    for (i in 1:nrow(new.vardata)) {
        
        # creating the variance matrix for ith data point
        diag.var <- new.vardata[i, ]
        fin.data <- diag(diag.var, n, n)
        
        # output before scaling
        var.transform <- loadings %*% fin.data %*% t(loadings)
        
        # if scaling is used unscale PCA transformation before returning final result
        if ((((class(object) == "prcomp") && (!identical(object$scale, FALSE)))) || ((class(object) == "princomp") && any(object$scale != 
            1))) {
            scale <- outer(object$scale, object$scale)
            var.transform <- var.transform * scale
        }
        
        # If only.var is TRUE return output with covariances, else return only variances.
        if (!only.var) {
            PCA.var.data[i, , ] <- var.transform
        } else {
            PCA.var.data[i, ] <- diag(var.transform)
        }
    }
    
    PCA.var.data
} 
