#' @rdname fitEmulator

#' @export
fitEmulatorTemp <- function(inputs, outputs, prior.mean = "linear",                            
                            cor.function = corGaussian, phi.opt, sigmasq.opt = NULL, fn.grad = NULL, optim = TRUE,
                            MCMC.iterations = 50,   
                            phi.init, sigmasq.init,
                            MC.plot = FALSE, nugget = NULL,
                            optim.method = "BFGS",
                            optimise.interval, ...) {
    
    inputs <- data.frame(inputs)
    outputs <- data.frame(outputs)
    
    n.outputs <- ncol(outputs)
    n.inputs <-  ncol(inputs)
    n.train <- nrow(inputs)
    
    # matching "prior.mean" to strings (constant/linear) or allowing user to provide a CORRECT "formula"
    if (is.character(prior.mean)) {
        match <- pmatch(prior.mean, c("constant", "linear"))
        if (is.na(match)) 
            stop("prior.mean should be a string of either 'linear' or 'constant', or an object with class 'formula'")
        if (match == 1) 
            formula <- as.Formula(paste(paste(colnames(outputs), collapse = "+"),  "~1"))
        else 
            formula <- as.Formula(paste(paste(colnames(outputs), collapse = "+"), "~", paste(colnames(inputs), collapse = "+")))
    } 
    
    H.training <- model.matrix(formula, data = inputs)
    outputs <- as.matrix(outputs)
    
    if (!is.null(rownames(inputs)))
        rownames(outputs) <- rownames(inputs)
    
    phi.init.miss <- missing(phi.init) || is.null(phi.init)
    
    if (phi.init.miss)
        phi.init <- rep(0, ncol(as.matrix(inputs)))
    
    # Set likelihood function according to what needs to be estimated - phi, and/or sigma, or none
    
    if(optim)
        fn <- negLogLik
    else 
        fn <- negLogLikGrad
    
    
    # lmtemp <- lm(formula, data = data.frame(inputs, outputs))
    # log.sigmasq.init <- 2 * log(summary(lmtemp)$sigma)
    
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
    } 
    log.sigmasq.init <- NULL
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
        theta.MCMC <- theta.init[1:(length(theta.init) - length(log.sigmasq.init))]
        if (!is.null(scaled.index)) {
            theta.MCMC <- theta.MCMC
            theta.MCMC[scaled.index] <- theta.MCMC[scaled.index] + log(max.x - min.x)[scaled.index]
        }
        
        timing <- system.time({
            if (optim) {
                opt <- optim(theta.init, fn, gr = fn.grad, method = optim.method, inputs = scaled.inputs, H = H.training, 
                             outputs = outputs, cor.function = cor.function, nugget = nugget,...)
                theta.opt <- opt$par
                log.lik <- -opt$value
                # if (opt$value == negloglik.lim)
                #     warning(paste("likelihood equals zero - consider using different initial values or optim method"))
                # else if (opt$convergence != 0) 
                #     warning("optimisation for correlation parameters not converged - consider using different initial values. ", "Log likelihood value: ", log.lik)
                
            } else {
                # browser()
                # fn <- negLogLikGrad
                opt <- nlm(fn, p = theta.init, inputs = scaled.inputs, H = H.training, 
                           outputs = outputs, nugget = nugget, ..., print.level = 1)
                theta.opt <- opt$estimate
                log.lik <- -opt$minimum
            }
        })
    }
    
    # extract phi.opt and sigmasq.opt, unscaling if neccesary
    # if (phi.opt.miss) {
    phi.opt <- theta.opt[1:(length(theta.opt) - length(log.sigmasq.init))] # phi.opt <- theta.opt[1:n.inputs]
    if (!is.null(scaled.index)) {
        phi.opt.scaled <- phi.opt
        phi.opt[scaled.index] <- phi.opt[scaled.index] + log(max.x - min.x)[scaled.index]
    }
    # }
    
    if (is.null(names(phi.opt))) 
        names(phi.opt) <- colnames(inputs)
    
    # calculating fit using scaled phi and inputs in one case and unscaled in the other to make absolutely certain that the values 
    # are identical to those found using the optim. 
    fit <- makeEmulatorConstants(phi = phi.opt.scaled, inputs = scaled.inputs, H = H.training, outputs = outputs, 
                                 cor.function = cor.function, nugget = nugget, sigmasq.hat = sigmasq.opt, ...)
    
    # finish creating fit object
    class(fit) <- "emulatorFit"
    
    # if (phi.opt.miss || sigmasq.miss) {
    fit$opt <- opt
    
    fit$opt.convergence <- opt$convergence
    fit$opt.message <- opt$message
    fit$opt.counts <- opt$counts
    fit$time.taken <- timing[1]
    
    fit$phi.MCMC.unscaled <- theta.init
    fit$phi.MCMC.scaled <- theta.MCMC
    fit$phi.hat = phi.opt
    fit$training.inputs = as.matrix(inputs)
    # fit$log.lik <- log.lik
    fit$n.train <- n.train
    fit$n.outputs <- n.outputs
    fit$formula = formula
    fit$call <- match.call()  # enabling update function to work
    
    fit
} 
