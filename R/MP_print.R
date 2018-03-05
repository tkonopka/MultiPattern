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
  checkArgClass(x, "MultiPattern")
  message("MultiPattern configuration object\n")
  message(sprintf("%-5s", length(x$items)), "observations")
  message(sprintf("%-5s", length(x$data)), "data objects")
  for (i in seq_along(x$data)) {
    xd = x$data[[i]]
    message(paste0("      |- ", sprintf("%-12s", names(x$data)[i]), 
               " - ", sprintf("%6d", ncol(xd)), " features"))   
  }
  message(sprintf("%-5s", length(x$configs)), "analyses configurations")
}




##' Print a summary of a MP distances objects
##'
##' @param x a MultiPatternSimilarities object
##' @param ... additional arguments, not used
##' 
##' @export
print.MultiPatternSimilarities = function(x, ...) {  
  checkArgClass(x, "MultiPatternSimilarities")
  message("MultiPattern similarities\n")
  message("Similarities:\t\t", length(x))
}




##' Print a summary of a MP distances objects
##'
##' @param x a MultiPatternSettings object (list)
##' @param ... additional arguments, not used
##' 
##' @export
print.MultiPatternSettings = function(x, ...) {
  checkArgClass(x, "MultiPatternSettings")
  message("MultiPattern analysis settings\n")
  ## find length of names(x)
  xmax = max(nchar(names(x)))
  for (s in names(x)) {
    xs = paste(x[[s]], collapse=" ")
    message(sprintf(paste0("%", xmax, "s"),s), ":\t", xs)
  }
}

