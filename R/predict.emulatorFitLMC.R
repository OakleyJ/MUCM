#' @export
predict.emulatorFitLMC <- function(object, newdata, var.cov = FALSE, sd = TRUE, tol = -1e-11, ...) {
    l <- list(...)
    if (length(l) > 0) 
        stop("unused Argument (", names(l[1]), " = ", l[[1]], ")")
    
    # ensure new and old fits are backward compatible
    object <- updateFit(object)
    
    # set newdata to fitted data if missing
    if (missing(newdata) || is.null(newdata))
        newdata <- object$training.inputs
    
    # checking if data is a dataframe
    newdata <- as.data.frame(newdata)
    n.pred <- nrow(newdata)
    n.all <- (object$n.train + n.pred) * object$n.outputs
    
    # extract inputs from newdata
    input.names <- colnames(object$training.inputs)
    # check if all inputs are given
    if (!all(input.names %in% colnames(newdata))) 
        stop("one or more input variables are missing from the newdata matrix")
    
    new.inputs <- as.matrix(newdata[, input.names, drop = FALSE])
    
    # initialise MV.H.pred
    MV.H.pred <- matrix(0, nrow = (n.pred * object$n.outputs), ncol = ncol(object$H.training))
    rownames(MV.H.pred) <- rep(rownames(newdata), object$n.outputs)
    colnames(MV.H.pred) <- colnames(object$H.training)
    
    # overwriting MV.H with values from small.h
    for (i in 1:object$n.outputs) {
        MV.H.pred[(((i - 1)*n.pred) + 1):(i*n.pred) , (sum(object$n.regressors[0:(i - 1)]) + 1):(sum(object$n.regressors[1:i]))] <- 
            model.matrix(object$formula, data = newdata, rhs = i)
    }
    
    cov.training.prediction <- t(do.call(buildCovMat, c(list(sigma = object$sigmasq.hat, phi = object$phi.hat, inputs = object$training.inputs, 
                                                             inputs2 = new.inputs, cor.function = object$cor.function),
                                                        object$cor.function.args)))
    
    post.mean <- MV.H.pred %*% object$betahat + cov.training.prediction %*% object$Vinv.e
    post.mean <- matrix(post.mean, nrow = n.pred, ncol = object$n.outputs)
    rownames(post.mean) <- rownames(new.inputs)
    colnames(post.mean) <- colnames(object$sigmasq.hat)
    predict <- list(posterior.mean = post.mean)
    
    if (var.cov || sd) {
        temp.constant <- MV.H.pred - cov.training.prediction %*% object$Vinv.H
        w3 <- backsolve(object$L, t(cov.training.prediction), transpose = TRUE) # solve(object$tL, t(cov.training.prediction))
        w4 <- backsolve(object$K, t(temp.constant), transpose = TRUE) # solve(object$K, t(temp.constant))
        
        if (var.cov) {
            V.pred <- do.call(buildCovMat, c(list(sigma = object$sigmasq.hat, phi = object$phi.hat, inputs = new.inputs, 
                                                  cor.function = object$cor.function),
                                             object$cor.function.args))
            
            predict$posterior.variance <- (V.pred - crossprod(w3) + crossprod(w4)) # (t(w3) %*% w3) + (t(w4) %*% w4))
            
            names <- paste(rep(colnames(object$sigmasq.hat), each = nrow(new.inputs)), rownames(new.inputs), sep = ":")
            dimnames(predict$posterior.variance) <- list(names, names)
            
            if ((min(diag(predict$posterior.variance), na.rm = TRUE) < 0) && !sd)
                warning(paste("Negative minimum posterior variance:", min(diag(predict$posterior.variance), na.rm = TRUE)))
        }
        
        if (sd) {
            if (var.cov == TRUE) {
                var <- diag(predict$posterior.variance)
            } else {    
                var <- rep(diag(object$sigmasq.hat), each = n.pred) - colSums(w3 ^ 2) + colSums(w4 ^ 2)
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
    predict
}

