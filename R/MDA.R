#' @title Maximum Dissimilarilty Algorithm
#' 
#' @description  Selects a subset of events by maximizing the dissimilarity between selected samples.
#'               Supports biased selection via event weighting.
#'               Does not require starting events to have been selected from same population or with same algorithm.
#' 
#' @param data.selected  a matrix (or data.frame) of events already selected.
#' @param data.new  a matrix (or data.frame) from which to select new events from.
#' @param n  the number of additional events to be selected from \code{data.new}.
#' @param weight.selected  optional argument (vector) used to weight the events already selected
#'        (paired with events in \code{data.selected}, and should have length = \code{nrow(data.selected)}).
#' @param weight.new  optional argument (vector) used to weight the events to be selected
#'        (paired with events in \code{data.new}, and should have length = \code{nrow(data.new)}).
#' @param normalise  whether to first normalise the data so each column is scaled between 0 and 1.
#'        Defaults to \code{TRUE} but use \code{FALSE} if already normalised.
#' @param index  a logical argument. If \code{FALSE} (default), the function returns a matrix (or data.frame) of selected events.
#'        If \code{TRUE} it returns indices of the selected events in \code{data.new}.
#' 
#' @return If \code{index = FALSE} the function returns a matrix (or data.frame) of selected events.
#'         Otherwise it returns the indices of the selected events in \code{data.new} in terms of the row number
#'         or the row name (if available).
#' 
#' @author Sajni Malde
#' 
#' @importFrom fields rdist
#' @export
MDA <- function(data.selected, data.new, n, weight.selected = NULL, weight.new = NULL, normalise = TRUE, index = FALSE) {
    
    # check n is sensible (n <= nrow(data.new))
    if (n > nrow(data.new))
        stop("n must be less than or equal to nrow(data.new)")
    
    if (n == 0)
        stop("n must be 1 or greater")
    
        # check data.selected and data.new have same dimension and colnames
    if (ncol(data.selected) != ncol(data.new))
        stop("data.selected and data.new have to have the same number of columns")
    if (any(colnames(data.selected) != colnames(data.new)))
        stop("data.selected and data.new don't have matching colnames")
    
    if (!is.null(weight.selected)) {
        # check dimensions of weight.selected
        if (is.vector(weight.selected) && length(weight.selected) != nrow(data.selected)) 
            stop("weight.selected should have length = nrow(data.selected)")
        else if (!is.vector(weight.selected))
            stop("weight.selected should be a vector")
    }
    
    if (!is.null(weight.new)) {
        # check dimensions of weight.new
        if (is.vector(weight.new) && length(weight.new) != nrow(data.new)) 
            stop("weight.new should have length = nrow(data.new)")
        else if (!is.vector(weight.new))
            stop("weight.new should be a vector")
    }
    
    # Normalising data
    if (normalise) {
        max.combined <- apply(rbind(apply(data.selected, 2, max), apply(data.new, 2, max)), 2, max)
        min.combined <- apply(rbind(apply(data.selected, 2, min), apply(data.new, 2, min)), 2, min)
        range.selected <- matrix(max.combined - min.combined, nrow(data.selected), ncol(data.selected), byrow = TRUE)
        range.new <- matrix(max.combined - min.combined, nrow(data.new), ncol(data.new), byrow = TRUE)
        min.selected <- matrix(min.combined, nrow(data.selected), ncol(data.selected), byrow = TRUE)
        min.new <- matrix(min.combined, nrow(data.new), ncol(data.new), byrow = TRUE)
        scaled.selected <- (data.selected - min.selected) / range.selected
        scaled.new <- (data.new - min.new) / range.new
    } else {
        scaled.selected <- data.selected
        scaled.new <- data.new
    }
    
    # Initialise output data
    selected.id <- NULL 
    
    data.new.ids <- 1:nrow(data.new)
    
    for (i in 1:n) {

        #calculate the distance
        diss <- rdist(scaled.selected, scaled.new[data.new.ids, ]) 
        
        # incorporate weighting if necessary
        if (!(is.null(weight.selected) && is.null(weight.new))) 
            diss <- tcrossprod(weight.selected, weight.new[data.new.ids]) * diss
        
        # choose maximum dissimilar point 
        chosen <- data.new.ids[which.max(apply(diss, 2, min))]
        
        # add it to selected id and subset (Not scaled!)
        selected.id <- c(selected.id, chosen)
        
        # add selected point to set 1 (and weight.selected) to prepare for next iteration
        scaled.selected <- rbind(scaled.selected, scaled.new[chosen, , drop = FALSE])
        
        if (!(is.null(weight.selected) && is.null(weight.new)))
            weight.selected <- c(weight.selected, weight.new[data.new.ids[data.new.ids %in% chosen]])
        
        data.new.ids <- setdiff(data.new.ids, chosen)
    }
    
    # set object to return
    if (index) {
        
        # return row numbers if data.new is a matrix
        ret <- selected.id
        
        # return names, if data.new is a data.frame
        if (!is.null(rownames(data.new)))
            ret <- rownames(data.new)[selected.id]
        
    } else{
        
        ret <- data.new[selected.id, , drop = FALSE]
        
    }
    
    ret 
}
