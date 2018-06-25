## Package: MultiPattern



##' check that an object is of the epected class, or stop
##'
##' @keywords internal
##' @param x object to check
##' @param expected character, expected class of object x
##' @param msg character, message printed at beginning of stop message
##'
checkArgClass= function(x, expected, msg="") {
  if (!expected %in% class(x)) {
    xcall = deparse(substitute(x))
    stop(paste0(msg, xcall, " should be of class ", expected, "\n"), call.=FALSE)
  }
}




##' check that object x is not null, or stop
##' 
##' @keywords internal
##' @param x object to check
##' @param msg character, message printed at beginning of stop message
##'
checkNotNull = function(x, msg) {
  if (is.null(x)) {
    stop(paste0(msg, " cannot be null\n"), call.=FALSE)
  }
}
