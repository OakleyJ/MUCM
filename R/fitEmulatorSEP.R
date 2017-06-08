#' @rdname fitEmulator

#' @export
fitEmulatorSEP <- function(inputs, outputs, prior.mean = "linear",                            
                           cor.function = corGaussian, phi.opt, sigmasq.opt, 
                           MCMC.iterations = 50,   
                           phi.init, sigmasq.init,
                           MC.plot = FALSE, nugget = NULL,
                           optim.method = "BFGS", optim.gr = NULL,
                           optimise.interval, ...) {
    
    inputs <- data.frame(inputs)
    outputs <- data.frame(outputs)
    
    # Modify output column names if equal to input column names
    if (length(intersect(colnames(inputs), colnames(outputs))) > 0) {
      colnames(outputs) <- paste(colnames(outputs),".out", sep = "")
    }
    
    n.outputs <- ncol(outputs)
    n.inputs <-  ncol(inputs)
    n.train <- nrow(inputs)
    
    # Checking if outputs and inputs have the same no of observations
    if (n.train != nrow(outputs))
        stop("number of rows in inputs does not match the number of rows in outputs")
    
    # Checkng if outputs have colnames otherwise improvise
    if (is.null(colnames(outputs))) 
        colnames(outputs) <- paste("Output", 1:n.outputs, sep = "")
    
    # Checking if nugget is correct length or if length 1, then code assumes same uncertainty for all observations
    if (!is.null(nugget)) {
        if (length(nugget) == 1)
            nugget <- rep(nugget, n.train)
        else if (length(nugget) != n.train)
            stop("nugget argument of incorrect length"  )
    }
    
    # Currently, nugget can only be included for univariate output
    if (n.outputs > 1 & !is.null(nugget)) 
        stop("currently, nugget can only be included for univariate output")
    
    # matching "prior.mean" to strings (constant/linear) or allowing user to provide a CORRECT "formula"
    if (is.character(prior.mean)) {
        match <- pmatch(prior.mean, c("constant", "linear"))
        if (is.na(match)) 
            stop("prior.mean should be a string of either 'linear' or 'constant', or an object with class 'formula'")
        if (match == 1) 
            formula <- as.Formula(paste(paste(colnames(outputs), collapse = "+"),  "~1"))
        else 
            formula <- as.Formula(paste(paste(colnames(outputs), collapse = "+"), "~", paste(colnames(inputs), collapse = "+")))
        
    } else if (any(class(prior.mean) == "formula")) {
        prior.mean <- as.Formula(prior.mean)
        # check all RHS variables are provided in inputs
        if (!all(all.vars(prior.mean[[3]]) %in% colnames(inputs))) 
            stop("'prior.mean' contains input variables not found in 'inputs'")
        
        # if no LHS provided, add it
        if (length(prior.mean)[1] == 0)
            formula <- as.Formula(paste(paste(colnames(outputs), collapse = "+"), paste(prior.mean, collapse = " ")))
        else {
            output.names <- all.vars(prior.mean[[2]])
            if (all(output.names %in% colnames(outputs))) 
                formula <- prior.mean
            else 
                stop("'prior.mean' contains output variables not found in 'outputs'")
        }
    } else 
        stop("prior.mean should be a string of either 'linear' or 'constant', or an object with class 'formula'")
    
    H.training <- model.matrix(formula, data = inputs)
    outputs <- as.matrix(outputs)
    
    if (!is.null(rownames(inputs)))
        rownames(outputs) <- rownames(inputs)
    
    phi.opt.miss <- missing(phi.opt) || is.null(phi.opt)
    phi.init.miss <- missing(phi.init) || is.null(phi.init)
    sigmasq.miss <- !is.null(nugget) && (missing(sigmasq.opt) || is.null(sigmasq.opt))
    phi.opt.scaled <- NULL # indicator if phi is scaled or not. Null = no scaling, !NULL = scaling
    
    # checking sigmasq.opt/init is of correct length
    if (is.null(nugget)) {
        if (!missing(sigmasq.opt) && !is.null(sigmasq.opt))
            warning("sigmasq.opt provided, but not used because nugget = NULL")
        if (!missing(sigmasq.init) && !is.null(sigmasq.init))
            warning("sigmasq.init provided, but not used because nugget = NULL")
        
        sigmasq.opt <- NULL
    } else {
        if (!sigmasq.miss && length(sigmasq.opt) != 1)
            stop("sigmasq.opt has incorrect dimensions")
        if ((!missing(sigmasq.init) && !is.null(sigmasq.init)) &&  length(sigmasq.init) != 1)
            stop("sigmasq.init has incorrect dimensions")
    }
    
    # if phi.opt and phi.init are provided, issue a warning 
    if (!phi.opt.miss) { 
        if (!phi.init.miss)
            warning("phi.init value ignored as phi.opt has been provided")
        
        # Checking if phi.opt is correct length or if length 1, then code assumes same phi for all observations
        if (length(phi.opt) == 1)
            phi.opt <- rep(phi.opt, n.inputs)
        else if (length(phi.opt) != n.inputs)
            stop("phi.opt argument of incorrect length")
    }
    
    # if sigmasq.opt provided, and sigmasq.init is not provided, issue a warning 
    if (!sigmasq.miss && !is.null(nugget) && (!missing(sigmasq.init) && !is.null(sigmasq.init)))
        warning("sigmasq.init value ignored as sigmasq.opt has been provided")
    
    # Checking whether phi or sigmasq needs to be estimated
    if (phi.opt.miss || sigmasq.miss) {
        
        # Set likelihood function according to what needs to be estimated - phi, and/or sigma, or none
        if (is.null(nugget)) {
            fn <- negLogLik
        } else {
            if (sigmasq.miss && phi.opt.miss)
                fn <- negLogLikNugget
            else if (!sigmasq.miss && phi.opt.miss) {
                # theta = phi
                fn <- function(theta, nugget, inputs, H, outputs, cor.function,  ...){
                    theta <- c(theta, log(sigmasq.opt))
                    negLogLikNugget(theta, nugget, inputs, H, outputs, cor.function, ...)
                }
            } else if (sigmasq.miss && !phi.opt.miss) {
                # theta = log(sigmasq)
                fn <- function(theta, nugget, inputs, H, outputs, cor.function, ...){
                    theta <- c(phi.opt, theta)
                    negLogLikNugget(theta, nugget, inputs, H, outputs, cor.function, ...)
                }
            }
        }
        
        # if phi.opt is given, then we dont need phi.init, otherwise we use a defult phi.init value, or user provided one
        if (!phi.opt.miss)
            phi.init <- NULL
        else if (phi.init.miss)
            phi.init <- rep(0, ncol(as.matrix(inputs)))
        else if (!phi.init.miss) {# check if you have the correct dimension
            # Checking if phi.init is correct length or if length 1, then code assumes same phi for all observations
            if (length(phi.init) == 1)
                phi.init <- rep(phi.init, n.inputs)
            else if (length(phi.init) != n.inputs)
                stop("phi.init argument of incorrect length"  )
        }
        
        # if sigmasq.opt is given, then we dont need sigmasq.init, otherwise we use a defult sigmasq.init value, or user provided one
        if (!sigmasq.miss)
            log.sigmasq.init <- NULL
        
        else if (!is.null(nugget)) {
            if (missing(sigmasq.init) || is.null(sigmasq.init)) {
                # Get starting value for sigmasq from linear model fit
                lmtemp <- lm(formula, data = data.frame(inputs, outputs))
                log.sigmasq.init <- 2 * log(summary(lmtemp)$sigma)
            } else 
                log.sigmasq.init <- log(sigmasq.init)
        }
        
        scaled.inputs <- as.matrix(inputs)
        scaled.index <- NULL
        
        # Scale inputs to [0, 1] for estimating parameters in the correlation function
        if (identical(cor.function, corGaussian) || identical(cor.function, corMatern2.5)) {
            min.x <- apply(inputs, 2, min)
            max.x <- apply(inputs, 2, max)
            range.x <- matrix(max.x - min.x, n.train, n.inputs, byrow = TRUE)
            min.x.matrix <- matrix(min.x, n.train, n.inputs, byrow = TRUE)
            scaled.inputs <- as.matrix((inputs - min.x.matrix)/range.x)
            scaled.index <- 1:n.inputs
            if (!phi.init.miss)
                phi.init[scaled.index] <- phi.init[scaled.index] - log(max.x - min.x)[scaled.index]
            
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
            scaled.inputs[, scaled.index] <- as.matrix(inputs[, scaled.index])
            if (!phi.init.miss)
                phi.init[scaled.index] <- phi.init[scaled.index] - log(max.x - min.x)[scaled.index]
        }
        
        theta.init <- c(phi.init, log.sigmasq.init)
        
        # estimate correlation parameters - optimisation bit
        # if optim.method is "brent" then MCMC is unnecessary as starting values not used. 
        if (optim.method == "Brent") {
            opt <- optimise(fn, optimise.interval, inputs = scaled.inputs, H = H.training, 
                         outputs = outputs, cor.function = cor.function, nugget = nugget, ...)
            
            theta.opt <- opt$minimum
            log.lik <- -opt$objective 
            
            if (opt$objective == negloglik.lim)
                warning(paste("likelihood equals zero - consider using different optim method or altering search interval"))
            else if (any(theta.opt == optimise.interval))
                warning("optimise.interval limit reached - consider widening interval")
            
        } else {
            # option to run mcmc if nec
            if (MCMC.iterations > 0) {
                # Estimate parameters in the correlation function -MCMC to get the starting values
                MCMCoutput <- MCMCMetropolisGibbs(inputs = scaled.inputs, 
                                                  outputs = outputs, fn = fn, 
                                                  H = H.training, MCMC.iterations = MCMC.iterations,
                                                  starting.values = theta.init, 
                                                  cor.function = cor.function, 
                                                  MC.plot = MC.plot, nugget = nugget, ...)
                index <- which(MCMCoutput$density.sample == max(MCMCoutput$density.sample))[1]
                
                if (all(theta.init == MCMCoutput$theta.sample[index, ]))
                    warning("MCMC output equals default value of theta.init. You may want to increase MCMC.iterations or provide a different theta.init.")
                
                theta.init <- MCMCoutput$theta.sample[index, ]  # initial values of optimisation
            }
            
            opt <- optim(theta.init, fn, method = optim.method, inputs = scaled.inputs, H = H.training, 
                         outputs = outputs, cor.function = cor.function, nugget = nugget,
                         gr = optim.gr, ...)
            theta.opt <- opt$par
            log.lik <- -opt$value
            
            if (opt$value == negloglik.lim)
                warning(paste("likelihood equals zero - consider using different initial values or optim method"))
            else if (opt$convergence != 0) 
                warning("optimisation for correlation parameters not converged - consider using different initial values. ", "Log likelihood value: ", log.lik)
        }
        
        # extract phi.opt and sigmasq.opt, unscaling if neccesary
        if (phi.opt.miss) {
            phi.opt <- theta.opt[1:(length(theta.opt) - length(log.sigmasq.init))] # phi.opt <- theta.opt[1:n.inputs]
            if (!is.null(scaled.index)) {
                phi.opt.scaled <- phi.opt
                phi.opt[scaled.index] <- phi.opt[scaled.index] + log(max.x - min.x)[scaled.index]
            }
        }
        if (sigmasq.miss)
            sigmasq.opt <- exp(tail(theta.opt, 1))
    }
    if (is.null(names(phi.opt))) 
        names(phi.opt) <- colnames(inputs)
    
    # calculate log likelihood if necessary
    if (!phi.opt.miss && !sigmasq.miss)
        log.lik <- -negLogLik(theta = c(phi.opt, sigmasq.opt), inputs = inputs, H = H.training, outputs = outputs, cor.function = cor.function, ...)
    
    # calculating fit using scaled phi and inputs in one case and unscaled in the other to make absolutely certain that the values 
    # are identical to those found using the optim. 
    if (!is.null(phi.opt.scaled)) {
        fit <- makeEmulatorConstants(phi = phi.opt.scaled, inputs = scaled.inputs, H = H.training, outputs = outputs, 
                                     cor.function = cor.function, nugget = nugget, sigmasq.hat = sigmasq.opt, ...)
    } else {
        # Calculate various constants in the emulator
        fit <- makeEmulatorConstants(phi = phi.opt, inputs = as.matrix(inputs), H = H.training, outputs = outputs, 
                                     cor.function = cor.function, nugget = nugget, sigmasq.hat = sigmasq.opt, ...)
    } 
    
    # finish creating fit object
    class(fit) <- "emulatorFit"
    
    if (phi.opt.miss || sigmasq.miss) {
        fit$opt.convergence <- opt$convergence
        fit$opt.message <- opt$message
        fit$opt.counts <- opt$counts
    }
    
    fit$phi.hat = phi.opt
    fit$training.inputs = as.matrix(inputs)
    fit$log.lik <- log.lik
    fit$n.train <- n.train
    fit$n.outputs <- n.outputs
    fit$formula = formula
    fit$call <- match.call()  # enabling update function to work
    
    fit
} 
