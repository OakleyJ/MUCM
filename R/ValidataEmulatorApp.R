#' validating the emulator emulator
#' @param emulator A fit object of class inheriting from \code{'emulatorFit'}.
#' @param new.inputs A data matrix of input(s) at which emulation is desired (new inputs). 
#'        Must contain at least all parameters given in \code{emulator$training.inputs}.
#'        (Ensure they are the same used to generate the \code{emulator.predictions}.)
#' @param new.outputs A matrix outputs used to validate the emulator predictions. This Matrix must have columns equal to the number of outputs. 
#' @param emulator.predictions Output generated from the \code{\link{predict.emulatorFit}} function. It is an optional argument. 
#'        If this is not provided it will be calculated using arguments provided. Note Posterior Variance will be calculated, unless you provide this argument. 
#' @param verbose Defaults to TRUE. If FALSE, text under graphs (explanations) will not appear. 
#' @examples
#' # Fit the emulator
#' fit <- fitEmulator(inputs = surfebm[1:25, 1:2], 
#'                    outputs = surfebm[1:25, 3, drop = FALSE])
#' 
#' # Use fitted emulator to predict posterior means and variances at the new points
#' predictions <- predict(fit, surfebm[26:35, 1:2], sd = FALSE, var.cov = TRUE)
#' 
#' # Compare predictions with true values for the new inputs
#' # Can also compare accuracy of prediction based on posterior variance
#' validateEmulatorApp(fit, surfebm[26:35, 1:2], surfebm[26:35, 3], predictions, verbose = TRUE, launch.browser = TRUE)
#' 
#' @export
validateEmulatorApp <- function(emulator, new.inputs, new.outputs, emulator.predictions = NULL, verbose = TRUE, ...){
    
    if (!require(shiny)) 
        stop("this function requires R package shiny to be installed")
    
    if (!require(plotly)) 
        stop("this function requires R package plotly to be installed")
    
    if (!require(shinythemes, quietly = TRUE)) {
        theme <- NULL
        warning("This app will look much better if you install the package 'shinythemes'")
    } else {
        theme <- shinytheme("united")
    }
    
    # ensure new and old fits are backward compatible
    emulator <- updateFit(emulator)
    
    # if emulator.predictions is not given, calculate it with var.cov = T
    # determine whether it is possible to get post.var, 
    if (is.null(emulator.predictions)) {
        emulator.predictions <- predict(emulator, new.inputs, var.cov = T)
        post.var.given <- TRUE
    } else if (!is.null(emulator.predictions$posterior.variance))
        post.var.given <- TRUE
    else if (!is.null(emulator.predictions$sigmasq) && !is.null(emulator.predictions$correlation.matrix))
        post.var.given <- TRUE
    else 
        post.var.given <- FALSE
    
    # extract output from newdata
    new.outputs <- as.matrix(new.outputs)
    
    n.validation <- nrow(new.outputs)
    
    # check if training.output and validation output are of same size
    if (emulator$n.outputs != ncol(new.outputs)) 
        stop("new.outputs is not the same size as outputs used in fitEmulator pls check new.outputs with emulator.predictions$posterior.mean")
    
    # creating output variable list with names and indices
    y.axis.index <- 1:emulator$n.outputs
    names(y.axis.index) <- colnames(emulator$training.outputs)
    
    x.axis.index <- 1:ncol(emulator$training.inputs)
    names(x.axis.index) <- colnames(emulator$training.inputs)
    
    # creating output variable list with names and indices
    output.names.index <- 1:emulator$n.outputs
    names(output.names.index) <- colnames(emulator$training.outputs)
    
    
    ui <- fluidPage(
        
        # setting theme
        theme = theme,
        
        #include JavaScript
        includeScript(system.file("hover_all.js", package = "MUCM")),
        
        
        fluidRow(
            #header
            headerPanel("Validating the Emulator Fit"),
            
            
            # intro paragraph
            h4("This is where I will give a brief intro about this app. This app works better if you include posterior.variance in the emulator.predictions"),
            
            selectInput("selected.output", "Select output variable of interest", names(output.names.index)),
            
            # div( # or u can use shinydashboard::box()
                textOutput("RMSE"), 
                textOutput("NormRMSE"), 
                textOutput("Coverage"),
            # ),
            
            column(5, plotlyOutput("Approx.vs.Truth")),
            
            column(2, selectInput("y.variable", "Select variable for Y axis", names(y.axis.index)),
                   selectInput("x.variable", "Select variable for X axis", names(x.axis.index))),
            
            column(5, plotlyOutput("training.vs.prediction.plot")),
            
            # explanation of the first 2 plots
            if (verbose) {
                fluidRow(
                    column(5, h4("Plot showing emulator predicted output against the true output. 
                          The error bars are the 95th percent confidence intervals of the predicted mean. 
                          Ideally, all the points should be close to the diagonal line, which represents the one to one line, and the points should have small error bars. 
                          If you select points using either box select or lasso select, you will be able to see this point on the plot to the left highlighted in red.")),
                    
                    # column(2),
                    
                    column(5, offset = 2, h4("Plot showing the training input and output variables based on your choice. 
                         Select points on the right hand plot (validation points) to see where they appear in this space."))
                )
            },
            
            # plot qq, pcd and MD plot
            if (post.var.given) {
                fluidRow(
                    column(4, plotlyOutput("qqplot")),
                    
                    column(4, plotlyOutput("pcd.plot")),
                    
                    column(4, plotlyOutput("MD.plot"))
                )},
            
            # explation of the qq, pcd and MD plot
            if (post.var.given && verbose) {
                
                fluidRow(
                    column(4, h4("This is a QQ plot of Cholesky residuals. If the points lie close to the 45 degree line then the normality assumptions for the simulator output is reasonable. If the gradient of the points is greater than 1 (less than 1), it suggests that the predictive variability was underestimated (over estimated). 
                            Curvature in the plot indicates nonnormality, and outliers at either end of the plot suggests local fitting problems or nonstationarity.")),
                    
                    column(4, h4("This is a Pivoted cholesky prediciton errors against the pivoting index. The pivoting index gives the order of the Pivoted cholesky prediciton errors with the largest conditional predictive variance. 
                                No patterns are expected. Too many large errors indicate an underestimation of variance and vice versa.Both cases can aslo suggest a nonstationary process. 
                                Either large or very small errors at the beginning of the plot (i.e., on the left side) indicates poor estimation of predictive variance or nonstationarity. 
                                However, large (or very small) errors at the end of the plot (i.e., on the right side) indicates overestimation (or underestimation) of the correlation length parameters, or that the chosen correlation structure is unsuitable. ")),
                    
                    column(4, h4("This is a plot of a normal distribution curve. The red line on the plot is the value of the calculated Mahalanobis distance taking into account correlation among outputs. Ideally it should be between the boundaries of the curve. Extreme values (large or small) indicate a conflict between the emulator and simulator. "))
                )
            }
        )
    )
    
    
    server <- function(input, output) {
        
        selectedOutput <- reactive({
            match(input$selected.output, names(output.names.index))
        })
        
        postMean <- reactive({
            if (emulator$n.outputs == 1) {
                post.mean <- emulator.predictions$posterior.mean
            } else { 
                post.mean <- emulator.predictions$posterior.mean[, selectedOutput(), drop = FALSE] 
            }
            post.mean
        })
        
        postVar <- reactive({
            if (emulator$n.outputs == 1) {
                if (!is.null(emulator.predictions$posterior.variance))
                    post.var <- emulator.predictions$posterior.variance
                else 
                    post.var <- NULL 
            } else if (emulator$n.outputs != 1) {
                if (!is.null(emulator.predictions$posterior.variance)) {
                    post.var <- emulator.predictions$posterior.variance[(((selectedOutput() - 1)*emulator$n.train) + 1):(selectedOutput()*emulator$n.train),
                                                                        (((selectedOutput() - 1)*emulator$n.train) + 1):(selectedOutput()*emulator$n.train)]
                } else
                    post.var <- NULL
            }
            post.var
        })
        
        emPredSd <- reactive({
            
            if (!is.null(postVar()))
                em.pred.sd <- matrix(sqrt(diag(postVar())), ncol = 1)
            else {
                if (emulator$n.outputs == 1) {
                    if (!is.null(emulator.predictions$standard.deviation))
                        em.pred.sd <- emulator.predictions$standard.deviation[, 1, drop = FALSE]
                    else 
                        em.pred.sd <- NULL
                } else 
                    if (!is.null(emulator.predictions$standard.deviation))
                        em.pred.sd <- emulator.predictions$standard.deviation[, selectedOutput(), drop = FALSE]
                    else 
                        em.pred.sd <- NULL
            }
            em.pred.sd
        })
        
        newOutputs <- reactive({
            if (emulator$n.outputs != 1) 
                new.outputs <- new.outputs[, selectedOutput(), drop = FALSE]
            else 
                new.outputs <- new.outputs
            
            new.outputs
        })
        
        # calculate RMSE and NormRMSE in all cases
        residuals <- reactive(newOutputs() - postMean())
        
        rmse <- reactive(sqrt(mean((residuals()) ^ 2)))
        normrmse <- reactive(rmse()/(max(newOutputs()) - min(newOutputs())))
        
        coverage <- reactive(
            if (!is.null(emPredSd()))
                mean(abs(residuals()) < qnorm(0.975) * emPredSd(), na.rm = TRUE)
        )
        
        output$RMSE <-     renderText(paste("RMSE: ", format(rmse()), sep = ""))
        output$NormRMSE <- renderText(paste("Normalised RMSE: ", format(normrmse()), sep = ""))
        output$Coverage <- renderText(paste("Coverage: ", format(coverage()), sep = ""))
        
        # approximations against truth plot (with error bars if possible)
        output$Approx.vs.Truth <- renderPlotly({
            data <- data.frame("Approximations" = postMean()[,1], truth = newOutputs()[,1], sd = emPredSd()[,1])
            
            p <- plot_ly(data, source = "Approx.vs.Truth")
            if (!is.null(emPredSd()))
                p <- p %>% add_markers(x = ~truth, y = ~Approximations, type = "scatter", mode = "markers", error_y = ~list(value = ~sd))
            else 
                p <- p %>% add_markers(x = ~truth, y = ~Approximations, type = "scatter", mode = "markers")
            
            p <- p %>% add_trace(x = ~truth, y = ~truth, type = "scatter", mode = "lines", hoverinfo = "skip", color = I("orange")) %>%
                layout(showlegend = FALSE, title = "Emulator Approximations against Truth",
                       xaxis = list(title = "Truth"),
                       yaxis = list(title = "Approximations"))
            p
        })
        
        # Plot for training and prediction points
        output$training.vs.prediction.plot <- renderPlotly({
            eventdata.selected.approx.truth <- event_data("plotly_selected", source = "Approx.vs.Truth")
            eventdata.selected.qqplot <- event_data("plotly_selected", source = "qqplot")
            eventdata.selected.pcd.plot <- event_data("plotly_selected", source = "pcd.plot")
            
            data <- as.data.frame(cbind("x.var" = emulator$training.inputs[, input$x.variable], "y.var" = emulator$training.outputs[, input$y.variable]))
            
            p <- plot_ly(data, x = ~x.var, y = ~y.var, type = "scatter", mode = "markers", name = "Training<br>Points", source = "training.vs.prediction.plot") %>% #, hoverinfo = 'text', text = __)
                # add_data(x = ~truth, y = ~truth, mode = "lines", hoverinfo = "skip") %>%
                layout(showlegend = TRUE, title = "Emulator Training Points" ,
                       xaxis = list(title = input$x.variable),
                       yaxis = list(title = input$y.variable) )
            
            if (!is.null(eventdata.selected.approx.truth)) {
                new.data <- as.data.frame(cbind("x.var" = new.inputs[(1 + eventdata.selected.approx.truth$pointNumber), input$x.variable], 
                                                "y.var" = new.outputs[(1 + eventdata.selected.approx.truth$pointNumber), match(input$y.variable, names(y.axis.index))]))
                p <- p  %>% add_markers(data = new.data, color = I("red"), name = "Validation<br>points<br>(Approx<br>Vs Truth)")#, hoverinfo = "skip") 
            }
            
            if (!is.null(eventdata.selected.qqplot)) {
                new.data <- as.data.frame(cbind("x.var" = new.inputs[(1 + eventdata.selected.qqplot$pointNumber), input$x.variable], 
                                                "y.var" = new.outputs[(1 + eventdata.selected.qqplot$pointNumber), match(input$y.variable, names(y.axis.index))]))
                p <- p  %>% add_markers(data = new.data, color = I("darkgreen"), name = "Validation<br>points<br>(QQ plot)")#, hoverinfo = "skip") 
            }
            
            if (!is.null(eventdata.selected.pcd.plot)) {
                new.data <- as.data.frame(cbind("x.var" = new.inputs[(1 + eventdata.selected.pcd.plot$pointNumber), input$x.variable], 
                                                "y.var" = new.outputs[(1 + eventdata.selected.pcd.plot$pointNumber), match(input$y.variable, names(y.axis.index))]))
                p <- p  %>% add_markers(data = new.data, color = I("yellow"), name = "Validation<br>points<br>(PCPE) ")#, hoverinfo = "skip") 
            }
            
            
            p %>% layout()
            
        })
        
        # Calculate cholecky residuals for next 2 plot
        cholResid <- reactive({
            if (!is.null(postVar())) {
                GP.resid <- newOutputs() - postMean()
                C <- chol(postVar(), pivot = TRUE)
                chol.resid <- solve(C, GP.resid[attr(C, "pivot")], tol = 1e-100)
                chol.resid
            } else  
                chol.resid <- NULL
        })
        
        # QQplot
        output$qqplot <- renderPlotly({
            if (!is.null(postVar())) {
                orig.order <- sapply(strsplit(names(cholResid()), ":", fixed = TRUE), function(x){x[2]})
                n <- length(cholResid())
                data <- data.frame("x" = qnorm(ppoints(n))[order(order(cholResid()))], "y" = cholResid())
                data <- data[match(rownames(postMean()), orig.order),]
                
                p <- plot_ly(data, x = ~x, y = ~y, type = "scatter", mode = "markers", source = "qqplot") %>% 
                    add_lines(x = extendrange(data[,"x"]), y = extendrange(data[,"x"]), color = I("orange"), type = "line") %>% 
                    layout(showlegend = FALSE, title = "Normal Q-Q plot",
                           xaxis = list(title = "Theoretical Quantiles"),
                           yaxis = list(title = "Sample Quantiles"))
                p
            } 
        })
        
        # plot of pivoted Cholesky errors
        output$pcd.plot <- renderPlotly({
            
            orig.order <- sapply(strsplit(names(cholResid()), ":", fixed = TRUE), function(x){x[2]})
            
            data <- data.frame("index" = 1:nrow(emulator.predictions$posterior.mean), "chol" = cholResid())
            data <- data[match(rownames(postMean()), orig.order),]
            
            x.range <- extendrange(data$index, f = 0.1)
            y.range <- c(1,1)
            
            p <- plot_ly(data, x= ~index, y = ~chol, type = "scatter", mode = "markers", source = "pcd.plot") %>% 
                add_lines(x = x.range, y.range*qnorm(0.975), name = "95% CI", linetype = 1, color = I("orange"), hoverinfo = "skip") %>% 
                add_lines(x = x.range, y.range*-qnorm(0.975), name = "95% CI", linetype = 1, color = I("orange"), hoverinfo = "skip") %>% 
                add_lines(x = x.range, y.range*0, linetype = 2, name = "95% CI", color = I("orange"), hoverinfo = "skip") %>% 
                layout(showlegend = FALSE, title = "Pivoted Cholesky Prediction Errors" ,
                       xaxis = list(title = "Pivot Index"),
                       yaxis = list(title = "Pivoted Cholesky Residuals"))
            p            
        })
        
        # plot of Mahalonobis distance
        output$MD.plot <- renderPlotly({
            GP.resid <- newOutputs() - postMean()
            MD.stat <- t(GP.resid) %*% solve(postVar(), GP.resid, tol = 1e-100) * 
                (emulator$n.train - emulator$n.regressors)/n.validation/(emulator$n.train - emulator$n.regressors - 2)
            lower.tail <- seq(0, qf(0.025, n.validation, emulator$n.train - emulator$n.regressors), length = 100) 
            middle <- seq(max(lower.tail), qf(0.975, n.validation, emulator$n.train - emulator$n.regressors) * 1.1, length = 100)
            upper.tail <- seq(max(middle), qf(0.9999, n.validation, emulator$n.train - emulator$n.regressors), length = 100)
            x.plot <- c(lower.tail, middle, upper.tail)
            
            if (!(MD.stat > max(x.plot) && MD.stat < min(x.plot)))
                x.plot <- c(x.plot, seq(max(x.plot), max(MD.stat) * 1.1, length = 100))
            
            data <- data.frame("x.plot" = x.plot, "Density" = df(x.plot, n.validation, emulator$n.train - emulator$n.regressors))
            
            p <- plot_ly(data, x= ~x.plot, y = ~Density, type = "scatter", mode = "lines") %>% 
                add_polygons(x = c(data[1:100, "x.plot"], rev(data[1:100, "x.plot"])), y = c(data[1:100, "Density"], rep(0,100)), color = I("steelblue2")) %>% 
                add_polygons(x = c(data[201:300, "x.plot"], rev(data[201:300, "x.plot"])), y = c(data[201:300, "Density"], rep(0,100)), color = I("steelblue2")) %>% 
                add_lines(x = c(MD.stat, MD.stat), y = range(data$Density), color = I("Orange")) %>% 
                layout(showlegend = FALSE, title = "Mahalanobis Distance Test",
                       xaxis = list(title = "x"),
                       yaxis = list(title = "Density"))
        })
        
    }
    
    runApp(list(ui = ui, server = server), ...)
    
}
