## fitted.emulatorFit : fitted values method for emulatorFit objects
##


#' @rdname fitEmulator
#' @param object  A fit object of class inheriting from \code{'emulatorFit'}.
#' @return  \code{fitted} returns the fitted values (i.e. the posterior mean estimate).
#'          This should equal \code{outputs} unless \code{nugget} is used.
#' @export
fitted.emulatorFit <- function(object, ...) {
    predict(object, var.cov = FALSE, sd = FALSE)$posterior.mean
}

