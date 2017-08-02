#' Generate a report for the GPE analysis
#' @param emulator A fit object of class inheriting from \code{'emulatorFit'}.
#' @param new.inputs A data matrix of input(s) at which emulation is desired (new inputs). 
#'        Must contain at least all parameters given in \code{emulator$training.inputs}.
#'        (Ensure they are the same used to generate the \code{emulator.predictions}.)
#' @param new.outputs A matrix outputs used to validate the emulator predictions. This Matrix must have columns equal to the number of outputs. 
#' @param emulator.predictions Output generated from the \code{\link{predict.emulatorFit}} function. It is an optional argument. 
#'        If this is not provided it will be calculated using arguments provided. Note Posterior Variance will be calculated, unless you provide this argument. 
#' @param ... Additional arguments to be passed to the \code{\link{render}} function (from rmarkdown).
#' @export
generateReport <- function(fit, ...){
    
    if (!require(rmarkdown)) 
        stop("this function requires R package rmarkdown to be installed")
    if (!require(pander)) 
        stop("this function requires R package pander to be installed")
    
    
    # TODO: to loop over all outputs for the validation bit. 
    
    path <- render(input="GPE_FitSummary.Rmd", output_file ="summaryReport.html", ...)
    message("File saved to ", path)
    
}