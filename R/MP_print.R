## Print functions for summarizing information in objects from the MultiMetric package
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
    ##cat("\n")
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
    cat("\nMP similarities\n\n")
    cat("Similarities:\t\t", length(x), "\n")
    cat("\n")
}



