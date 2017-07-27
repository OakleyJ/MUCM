#' @rdname fitEmulator
#' @param x A fit object of class inheriting from \code{'emulatorFit'}.
#' @export
print.emulatorFit <- function(x, ...) {
    # print the call funciton
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
    if (all(dim(x$sigmasq) > 0)) {
        if (all(dim(x$sigmasq) == c(1, 1))) {
            cat("Sigmasq: ", x$sigmasq.hat, "\n")
        } else {
            cat("Sigmasq: \n")
            print(format(x$sigmasq.hat), print.gap = 2L, quote = FALSE)
        }
    } else {
        cat("No Sigmasq.Hat\n")}
    
    if (any(class(x) == "emulatorFitLMC")) {
        #unstack beta
        beta.mat <- matrix(x$betahat, nrow = length(x$betahat)/x$n.outputs, ncol = x$n.outputs)
        dimnames(beta.mat) <- list(c(rownames(x$betahat[1:(length(x$betahat)/x$n.outputs), ,drop = FALSE])), rownames(x$sigmasq.hat))
        cat("\nPosterior mean of beta:\n")
        print(format(beta.mat), print.gap = 2L, quote = FALSE)
    } else if (all(dim(x$betahat) > 0)) {
        cat("\nPosterior mean of beta:\n")
        print(format(x$betahat), print.gap = 2L, quote = FALSE)
    } else cat("No Beta calculated\n")
    
    cat("\nNumber of outputs: ", x$n.outputs, "\n")
    
    cat("\nMLE of correlation parameters:\n")
    print(exp(x$phi.hat), print.gap = 2L, quote = FALSE)
    
    cat("\nMaximum log likelihood:\n")
    print(x$log.lik, print.gap = 2L, quote = FALSE)
    
    if (!is.null(x$opt.convergence)) {
        cat("\nOptimisation information:\n")
        if (x$opt.convergence == 0) {
            conv <- "Successful"
        } else if (x$opt.convergence == 1) {
            conv <- "Iteration limit reached"
        } else conv <- cat(x$opt.convergence, "- optimisation may not have succeeded")
        cat("Convergence: ", conv, "\n")
        if (!is.null(x$opt.message)) {
            cat("Message:", x$opt.message, "\n")
        }
    }
    
    cat("\n")
    invisible(x)
} 
