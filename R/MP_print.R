## Package: MultiPattern
##
## Print functions for summarizing information in package objects 
##




##' Print a summary of a MP configuration object
##'
##' @param x a MultiPattern configuration object
##' @param ... additional arguments, not used
##' 
##' @export
print.MultiPattern = function(x, ...) {    
    if (class(x) != "MultiPattern") {
        stop("object not of class MultiPattern\n")
    }
    capturex = as.character(deparse(substitute(x)))
    cat("MultiPattern configuration object\n\n")
    cat(sprintf("%-5s", length(x$items)), "observations\n")
    cat(sprintf("%-5s", length(x$data)), "data objects\n")
    for (i in seq_along(x$data)) {
        xd = x$data[[i]]
        cat(paste0("      |- ", sprintf("%-12s", names(x$data)[i]), 
                   " - ", sprintf("%6d", ncol(xd)), " features\n"))   
    }
    cat(sprintf("%-5s", length(x$configs)), "analyses configurations\n")
}




##' Print a summary of a MP distances objects
##'
##' @param x a MultiPatternSimilarities object
##' @param ... additional arguments, not used
##' 
##' @export
print.MultiPatternSimilarities = function(x, ...) {
    if (class(x) != "MultiPatternSimilarities") {
        stop("object not of class MultiPatternSimilarities\n")
    }    
    cat("MultiPattern similarities\n\n")
    cat("Similarities:\t\t", length(x), "\n")
}




##' Print a summary of a MP distances objects
##'
##' @param x a MultiPatternSettings object (list)
##' @param ... additional arguments, not used
##' 
##' @export
print.MultiPatternSettings = function(x, ...) {
    if (class(x) != "MultiPatternSettings") {
        stop("object not of class MultiPatternSettings\n")
    }    
    cat("MultiPattern analysis settings\n\n")
    ## find length of names(x)
    xmax = max(nchar(names(x)))
    for (s in names(x)) {
        cat(sprintf(paste0("%", xmax, "s"),s), ":\t", x[[s]], "\n")
    }
}

