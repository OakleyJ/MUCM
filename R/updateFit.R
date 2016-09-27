#' @title Update fit
#' 
#' @description Updates old \code{"emulatorFit"} objects to work with latest 
#' code to allow for improvements to be made within the code. 
#' 
#' @param object  An \code{"emulatorFit"} object to update.
#' 
#' @return Returns the updated \code{"emulatorFit"} object.
#' 
#' @note  A warning is issued if the object is updated.
updateFit <- function(object) {
    # return if OK
    if (!is.null(object$L))
        return(object)
    
    # old fit had tL - new fit has L
    object$L <- t(object$tL)
    object$tL <- NULL
    
    # old fit had K = t(chol(Q)) - new fit has K = chol(Q)
    object$K <- t(object$K)
    
    # add old cor.function to args to be used by new wrapper function
    object$cor.function.args$old.cor.function <- object$cor.function
    
    # create wrapper for old cor.function to provide 2 input args
    cor.func <- function(inputs, inputs2 = NULL, phi, ...) {
        # extract additional arguments into list
        args <- list(...)
        
        # remove old cor function from arguments
        old.cor.func <- args$old.cor.function
        old.args <- args[names(args) != "old.cor.function"]
        
        # use old cor.function if second argument missing
        if (missing(inputs2) || is.null(inputs2)) {
            return(do.call(old.cor.func, c(list(inputs, phi = phi), old.args)))
        }
        
        # call old cor.function to calculate entire cov matrix of all inputs
        combined.inputs <- rbind(inputs2, inputs)
        combined.A <- do.call(old.cor.func, c(list(combined.inputs, phi = phi), old.args))
        
        # return requested subset
        n.all <- nrow(combined.inputs)
        n.2 <- nrow(inputs2)
        combined.A[(n.2 + 1):n.all, 1:n.2, drop = FALSE]
    }
    object$cor.function <- cor.func
    
    # issue warning
    warning("You are using an old version of the fit object.  Consider refitting to benefit from faster predictions.")
    
    # return modified fit
    object
}
