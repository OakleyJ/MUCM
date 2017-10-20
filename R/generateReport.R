#' Generate a report for the GPE analysis
#' @param emulator A fit object of class inheriting from \code{'emulatorFit'}.
#' @param output.filename A filename (including the full path and file extension) to where the output should be saved. 
#' @export
generateReport <- function(fit, output.filename){
    
    # load packages
    if (!require(rmarkdown)) 
        stop("this function requires R package rmarkdown to be installed")
    if (!require(pander)) 
        stop("this function requires R package pander to be installed")
    
    # save fit to a temporary directory
    temporary.dir <-  tempdir()
    filename <- file.path(temporary.dir, "TempFit.RData")
    save(fit, file = filename)
 
    # Copy the rmd template to temp.rmd.fil
    filepath <- file.path(temporary.dir, "GPE_FitSummary.Rmd")
    file.copy(from = system.file("rmd/GPE_FitSummary.Rmd", package = "MUCM"), to = filepath)
    
    # Render the moved RMD file
    output.file.name <- basename(output.filename)
    output.file.dir <- dirname(output.filename)
    path <- render(input = filepath, output_file = output.file.name, output_dir = output.file.dir)
    # message("File saved to ", path)
    
    # delete temporary folder
    unlink(temporary.dir)
    
    # TODO: to loop over all outputs for the validation bit. 
}