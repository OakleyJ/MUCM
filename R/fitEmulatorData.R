#' @title Fit Emulator (Data)

#' @description A function equivalent to \code{fitEmulator}. In \code{fitEmulatorData} the user can provide
#' the inputs and outputs in one matrix/dataframe. To provide them seperately, use \code{\link{fitEmulator}}

#' @param formula Formula to include the RHS and LHS. The LHS indicates the outputs and the RHS defines the prior mean function defined by the column names of \code{data}. 
#'        A multivariate formula can also be given (see \code{\link{as.Formula}}) . 
#' @param data A data frame or matrix containing at least the inputs and outputs of the model. 
#' @param input.names An optional vector of character containing the names of variables in data to be considered as input variables. 
#'        By default, all variables in the data matrix, except those on the left hand side of the formula, are treated as inputs. 
#' @param ... Additional arguments to be passed to the more specific \code{fitEmulator} function.
#'        Arguments such as \code{MCMC.iterations}, \code{phi.opt}, \code{MC.plot}, etc.   

#' @return \code{fitEmulatorData} returns an object of class 'emulatorFit'. 
#'          See \code{\link{fitEmulator}} for details.


#' @seealso
#' If input and output data is in different objects, see \code{\link{fitEmulator}} as a more conveniant method. 
#
#' @examples 
#' # specify the prior mean function as a object of class 'formula'
#' # ensure colnames of datasets include at least these variables.
#' emulator.prior.mean.formula <- Y ~ X.1 + X.2  
#' 
#' # fit the emulator using a prior mean specified above and a Gaussian correlation function.
#' fit <- fitEmulatorData(formula = emulator.prior.mean.formula, data = surfebm[1:25, ], 
#'                        cor.function = corGaussian)
#' 
#' # Use fitted emulator to predict posterior means and variances at the new inputs
#' predictions <- predict(fit, surfebm[26:35, 1:2], sd = FALSE, var.cov = TRUE)
#' 
#' # view results
#' print(predictions$posterior.mean)

#' @importFrom Formula as.Formula
#' @export
fitEmulatorData <- function(formula, data, input.names = NULL, ...) {
    
    # convert formula to Formula (MV)
    formula <- as.Formula(formula)
    
    # extract inputs and output from data. if LHS = 0 throw an error to say Y must be provided
    if (length(formula)[1] == 0) {
        stop("the output variable must be provided in formula")
    }
    
    # extract output names from left hand side of formula (without transformations e.g. log)
    output.names <- all.vars(formula[[2]])
    
    # set input names to all remaining variables, if missing
    if (is.null(input.names)) {
        input.names <- setdiff(colnames(data), output.names)
    }
    
    # check output names
    if (!all(output.names %in% colnames(data))) {
        stop("formula states output variables that are not given in dataset")
    }
    
    # check input names
    if (!all(input.names %in% colnames(data))) {
        stop("provided input variables are not given in dataset")
    }
    
    ## NOTE: should check all variables stated in RHS of formula are in input.names
    
    # extract outputs (potentially transformed)
    outputs <- model.frame(formula, rhs = 0, data = as.data.frame(data))
    
    # extract inputs
    inputs <- data[, input.names, drop = FALSE]
    
    # fit emulator
    fit <- fitEmulator(inputs = inputs, outputs = outputs, prior.mean = formula,   ...)
    
    # update call
    fit$call <- match.call()
    
    # return
    fit
} 
