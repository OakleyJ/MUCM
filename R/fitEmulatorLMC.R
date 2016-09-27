#' @rdname fitEmulator
#' @param param The type of parametrization chosen to appopriately parametrize the covariance matrix. 
#' param can take one of the following values: \code{'original'}, \code{'cholesky'}, \code{'log-cholesky'}, \code{'spherical'}, \code{'matrix-log'} (defult), which refer to the name of the method used to parametrise according to the paper (1996) listed in the reference.
#' @references Pinheiro, J. C., & Bates, D. M. (1996). Unconstrained parametrizations for variance-covariance matrices. Statistics and Computing, 6(3), 289-296.
#' @export
fitEmulatorLMC <- function(inputs, outputs, prior.mean = "linear",
                           cor.function = corGaussian, phi.opt, sigmasq.opt, 
                           MCMC.iterations = 50,   
                           phi.init, sigmasq.init,
                           MC.plot = FALSE, param = "matrix-log", ...) {
    
    inputs <- data.frame(inputs)
    outputs <- data.frame(outputs)
    n.outputs <- ncol(outputs)
    n.inputs <-  ncol(inputs)
    n.train <- nrow(inputs)
    n.sigma <- n.outputs * (n.outputs + 1)/2
    
    # Checking if outputs and inputs have the same no of observations
    if (n.train != nrow(outputs))
        stop("Number of rows in inputs does not match the number of rows in outputs.")
    
    # Checking if outputs have colnames otherwise throw error 
    #( TODO: rethink! - could be equal to lhs of formula)
    if (is.null(colnames(outputs))) 
        stop("colnames of outputs must be provided")
    
    # matching "prior.mean" to strings (constant/linear) or allowing user to provide a CORRECT "formula"
    if (is.character(prior.mean)) {
        match <- pmatch(prior.mean, c("constant", "linear"))
        if (is.na(match)) 
            stop("prior.mean should be a string of either 'linear' or 'constant', or an object with class 'formula'")
        if (match == 1) 
            formula <- as.Formula(paste(paste(colnames(outputs), collapse =  " | "),  "~", paste(rep("1", n.outputs), collapse = " | ")))
        else 
            formula <- as.Formula(paste(paste(colnames(outputs), collapse =  " | "), "~", 
                                        paste(rep(paste(colnames(inputs), collapse = " + "), n.outputs), collapse = " | ")))
        # if formula is provided... 
    } else if (any(class(prior.mean) == "formula")) {
        prior.mean <- as.Formula(prior.mean)
        # if no LHS provided, add it
        if (length(prior.mean)[1] == 0) {
            formula <- as.Formula(paste(paste(colnames(outputs), collapse = " | "), paste(prior.mean, collapse = " ")))
        } else if (length(prior.mean)[1] == n.outputs) { # checking the correct number of sub LHS's provided
            output.names <- all.vars(prior.mean[[2]]) # checking the names match colnames of output
            if (!all(output.names %in% colnames(outputs)))
                stop("some output names provided do not match those in colnames(outputs)")
            else 
                formula <- prior.mean    
        } else 
            stop("Length of LHS of formula should match ncol(outputs)")
        # checking RHS
        if (length(formula)[2] == n.outputs) {
            input.names <- all.vars(formula[[3]]) #checking the names match colnames of output
            if (!all(input.names %in% colnames(inputs))) 
                stop("some input names provided do not match those in colnames(inputs)")
        }
    } else 
        stop("prior.mean should be a string of either 'linear' or 'constant', or an object with class 'formula'")
    
    # create a MV H matrix
    # get the dimensions of MV.H
    small.h <- list()
    names <- NULL
    for (i in 1:n.outputs) {
        small.h[[i]] <- model.matrix(formula, data = inputs, rhs = i)
        names <- c(names, colnames(small.h[[i]]))
    }    
    n.regressors <- sapply(small.h, ncol)
    n.beta <- sum(n.regressors)
    
    # initialise MV.H
    MV.H <- matrix(0, nrow = (n.train * n.outputs), ncol = n.beta)
    rownames(MV.H) <- rep(rownames(inputs), n.outputs)
    colnames(MV.H) <- names
    # overwriting MV.H with values from small.h
    for (i in 1:n.outputs) {
        MV.H[(((i - 1)*n.train) + 1):(i*n.train) , (sum(n.regressors[0:(i - 1)]) + 1):(sum(n.regressors[1:i]))] <- small.h[[i]]
    }
    
    outputs <- as.matrix(outputs)
    inputs <- as.matrix(inputs)
    stacked.outputs <- matrix(as.matrix(outputs), ncol = 1)
    rownames(stacked.outputs) <- paste(rep(colnames(outputs), each = nrow(inputs)), rownames(inputs), sep = ":")
    
    phi.init.miss <- missing(phi.init) || is.null(phi.init)
    sigmasq.init.miss <- missing(sigmasq.init) || is.null(sigmasq.init)
    phi.opt.miss <- missing(phi.opt) || is.null(phi.opt)
    sigmasq.opt.miss <- missing(sigmasq.opt) || is.null(sigmasq.opt)
    
    # if phi.opt and phi.init are provided, issue a warning 
    if (!phi.opt.miss && !phi.init.miss)
        warning("phi.init value ignored as phi.opt has been provided.")
    
    # if sigmasq.opt and sigmasq.init are provided, issue a warning 
    if (!sigmasq.opt.miss && !sigmasq.init.miss )
        warning("sigmasq.init value ignored as sigmasq.opt has been provided.")
    
    # checking sigmasq.opt/init is of correct length
    if (!sigmasq.opt.miss && all(dim(sigmasq.opt) != c(n.outputs, n.outputs)))
        stop("sigmasq.opt has incorrect dimensions")
    else if (!sigmasq.opt.miss && !isSymmetric(sigmasq.opt))
        stop("sigmasq.opt must be symmetric")
    if (!sigmasq.init.miss && all(dim(sigmasq.init) != c(n.outputs, n.outputs)))
        stop("sigmasq.init has incorrect dimensions")
    else if (!sigmasq.init.miss && !isSymmetric(sigmasq.init))
        stop("sigmasq.init must be symmetric")
    
    if (!phi.opt.miss && ncol(phi.opt) != n.outputs)
        stop("phi.opt has incorrect number of columns")
    if (!phi.init.miss && ncol(phi.init) != n.outputs) 
        stop("phi.init has incorrect number of columns")
    
    # Checking whether phi or sigmasq needs to be estimated
    if (phi.opt.miss || sigmasq.opt.miss) {
        
        if (!phi.opt.miss) {
            phi.init <- NULL
            n.phi = nrow(phi.opt)
        } else if (phi.init.miss) {
            phi.init <- matrix(0, ncol = n.outputs, nrow = n.inputs)
            for (i in 1:ncol(outputs)) {
                fit.phi <- fitEmulatorSEP(inputs, outputs[,i, drop = FALSE], prior.mean = "linear", cor.function = cor.function, 
                                          MCMC.iterations = MCMC.iterations, MC.plot = FALSE, ...)            
                phi.init[,i] <- fit.phi$phi.hat
                cat("finding default values for phi... ", colnames(outputs)[i], "\n")
            }
        } 
        
        if (!is.null(phi.init))
            n.phi = nrow(phi.init)
        
        scaled.inputs <- inputs
        scaled.index <- NULL        
        
        if (phi.opt.miss) {
            # Scale inputs to [0, 1] for estimating parameters in the correlation function
            if (identical(cor.function, corGaussian) || identical(cor.function, corMatern2.5)) {
                min.x <- apply(inputs, 2, min)
                max.x <- apply(inputs, 2, max)
                range.x <- matrix(max.x - min.x, n.train, n.inputs, byrow = TRUE)
                min.x.matrix <- matrix(min.x, n.train, n.inputs, byrow = TRUE)
                scaled.inputs <- as.matrix((inputs - min.x.matrix)/range.x)
                scaled.index <- 1:n.inputs
                phi.init[scaled.index, ] <- phi.init[scaled.index, ] - log(max.x - min.x)[scaled.index]
            } else if (identical(cor.function, corCombined)) {
                l <- list(...)
                
                # checking cor.funcs is provided and if so checking its of correct length
                if (is.null(l$cor.funcs))
                    stop("cor.funcs argument missing")
                else if (length(l$cor.funcs) != n.inputs)
                    stop("cor.funcs has incorrect length. Length of cor.funcs should equal the number of inputs")
                
                min.x <- apply(inputs, 2, min)
                max.x <- apply(inputs, 2, max)
                range.x <- matrix(max.x - min.x, n.train, n.inputs, byrow = TRUE)
                min.x.matrix <- matrix(min.x, n.train, n.inputs, byrow = TRUE)
                scaled.inputs <- as.matrix((inputs - min.x.matrix)/range.x)
                scaled.index <- l$cor.funcs == "gp" | l$cor.funcs == "gp2"
                scaled.inputs[, scaled.index] <- inputs[, scaled.index]
                phi.init[scaled.index, ] <- phi.init[scaled.index, ] - log(max.x - min.x)[scaled.index]
            }
        } 
        
        # if sigmasq.opt is given, then we dont need sigmasq.init, otherwise we use a defult sigmasq.init value, or user provided one
        if (!sigmasq.opt.miss) {
            sigmasq.init <- NULL
        } else if (sigmasq.init.miss) {
            cat("finding default values for sigma", "\n")
            fit.sigma <- fitEmulatorSEP(inputs, outputs, cor.function = cor.function, 
                                        MCMC.iterations = MCMC.iterations, MC.plot = FALSE, ...)
            sigmasq.init <- fit.sigma$sigmasq.hat
        } 
        
        # Set likelihood function according to what needs to be estimated - phi, and/or sigma, or none
        if (sigmasq.opt.miss && phi.opt.miss) {
            # theta is phi + uppertriangular sigma (inc diagonal)
            fn <- function(theta.vect, inputs, H, outputs, cor.function, n.outputs, param, ...){
                # extracting phi matrix
                n.phi.total = (length(theta.vect) - (n.outputs * (n.outputs + 1)/2))
                phi.vect <- theta.vect[1:(n.phi.total)]
                phi <- matrix(phi.vect, ncol = n.outputs, nrow = (n.phi.total/n.outputs))
                # creating and extracting sigma matrix
                sigma <- thetaSigmaParam(theta = theta.vect[(n.phi.total + 1):length(theta.vect)], param = param, n.outputs = n.outputs)
                # Note: Going to check sigma is posive definite in negLogLikLMC
                if (any(is.nan(sigma)) || !all(is.finite(sigma)))
                    return(negloglik.lim)
                negLogLikLMC(phi, sigma, inputs, H, outputs, cor.function, ...) 
            }

            theta.vect.init <- c(as.vector(phi.init), sigmaThetaParam(sigmasq.init, param))
        } else if (!sigmasq.opt.miss && phi.opt.miss) {
            # theta = phi
            fn <- function(theta.vect, inputs, H, outputs, cor.function, n.outputs, param, ...){
                phi <- matrix(theta.vect, ncol = n.outputs, nrow = n.phi)
                negLogLikLMC(phi, sigma = sigmasq.opt, inputs, H, outputs, cor.function, ...) 
            } 
            theta.vect.init <- as.vector(phi.init)
        } else if (sigmasq.opt.miss && !phi.opt.miss) {
            # theta = sigma
            fn <- function(theta.vect, inputs, H, outputs, cor.function, n.outputs, param, ...) {
                sigma <- thetaSigmaParam(theta = theta.vect, param = param, n.outputs = n.outputs)
                # Note: Going to check sigma is posive definite in negLogLikLMC
                if (any(is.nan(sigma)) || !all(is.finite(sigma)))
                    return(negloglik.lim)
                negLogLikLMC(phi = phi.opt, sigma, inputs, H, outputs, cor.function, ...) 
            }
            theta.vect.init <- sigmaThetaParam(sigmasq.init, param)
        }
        
        # option to run mcmc if nec
        if (MCMC.iterations > 0) {
            # Estimate parameters in the correlation function -MCMC to get the starting values
            cat("MCMCing...", "\n")
            MCMCoutput <- MCMCMetropolisGibbs(inputs = scaled.inputs,
                                              outputs = stacked.outputs, fn = fn, param = param,
                                              H = MV.H, MCMC.iterations = MCMC.iterations,
                                              starting.values = theta.vect.init, 
                                              cor.function = cor.function, 
                                              MC.plot = MC.plot, n.outputs = n.outputs, ...)
            index <- which(MCMCoutput$density.sample == max(MCMCoutput$density.sample))[1] 

            if (all(theta.vect.init == MCMCoutput$theta.sample[index, ]))
                warning("MCMC output equals default value of theta.init. You may want to increase MCMC.iterations or provide a different theta.init.")
            theta.vect.init <- MCMCoutput$theta.sample[index, ]  # initial values of optimisation
        }
        
        #  estimate correlation parameters - optimisation bit
        cat("Optimising...", "\n")
        opt <- optim(theta.vect.init, fn, method = "BFGS", inputs = scaled.inputs, H = MV.H, param = param,
                     outputs = stacked.outputs, cor.function = cor.function, n.outputs = n.outputs, ...)
        
        theta.vect.opt <- opt$par
        log.lik <- -opt$value
        
        if (opt$value == negloglik.lim)
            warning(paste("Likelihood is equal zero - Consider using different initial values"))
        else if (opt$convergence != 0) 
            warning("Optimisation for correlation parameters not Converged - Consider using different initial values. ", "Log likelihood value: ", log.lik)
        
        # extracting phi matrix
        if (phi.opt.miss && sigmasq.opt.miss ) {
            phi.opt <- matrix(theta.vect.opt[1:(n.outputs*n.phi)], ncol = n.outputs, nrow = n.phi)
            if (!is.null(scaled.index))
                phi.opt[scaled.index, ] <- phi.opt[scaled.index, ] + log(max.x - min.x)[scaled.index]
            # creating and extracting sigmasq.opt matrix
            sigmasq.opt <- thetaSigmaParam(theta = theta.vect.opt[((n.outputs*n.phi) + 1):length(theta.vect.opt)], param = param, n.outputs = n.outputs)
        } else if (phi.opt.miss && !sigmasq.opt.miss ) {
            phi.opt <- matrix(theta.vect.opt, ncol = n.outputs, nrow = n.phi)
            if (!is.null(scaled.index))
                phi.opt[scaled.index, ] <- phi.opt[scaled.index, ] + log(max.x - min.x)[scaled.index]
        } else if (!phi.opt.miss && sigmasq.opt.miss ) { 
            # creating and extracting sigmasq.opt matrix
            sigmasq.opt <- thetaSigmaParam(theta = theta.vect.opt, param = param, n.outputs = n.outputs)
        }
    }
    
    colnames(sigmasq.opt) <- colnames(outputs)
    rownames(sigmasq.opt) <- colnames(outputs)
    colnames(phi.opt) <- colnames(outputs)
    rownames(phi.opt) <- colnames(inputs)
    
    if (!phi.opt.miss && !sigmasq.opt.miss)
        log.lik <- -negLogLikLMC(phi.opt, sigmasq.opt, inputs, MV.H, stacked.outputs, cor.function, ...)
    
    
    # Calculate various constants in the emulator
    fit <- makeEmulatorConstantsLMC(phi.opt, sigmasq.opt, inputs = inputs, H = MV.H, outputs = stacked.outputs, 
                                    cor.function = cor.function, ...)
    
    class(fit) <- c("emulatorFitLMC", "emulatorFit")
    
    if (phi.opt.miss || sigmasq.opt.miss) {
        fit$opt.convergence <- opt$convergence
        fit$opt.message <- opt$message
        fit$opt.counts <- opt$counts
    }
    
    fit$log.lik <- log.lik
    fit$n.regressors <- n.regressors
    fit$n.train <- n.train
    fit$n.outputs <- n.outputs
    fit$formula <- formula
    fit$sigmasq.param <- param
    # enabling update function to work
    fit$call <- match.call()
    
    fit
}

