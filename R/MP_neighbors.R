## Package: MultiPattern
##
## Functions for dealing with neighborhoods 
##




##' Compute non-parametric neighborhoods of data elements using ranks on a per column basis
##'
##' This function takes as input a distance object or distance matrix (values
##' describing distance between pairs of data elements). For each data element,
##' it considers what other points are nearby and ranks them by distance.
##'
##' The output is non-parametric neighborhood object. "x[i,j] = r" implies that
##' data element "j" is the "r-th" one to data element "i". 
##' 
##' @param dd - matrix of distances between neighbors
##' @param ties - character, used to determine how equal ranks are handled by rank function
##' @param ... - additional arguments passed on to rank
##'
##' @export
MPrankNeighbors = function(dd, ties="average", ...) {

  ## perhapss convert fancy objects into a matrix
  if (class(dd)=="data.frame" | class(dd)=="data.table" | class(dd)=="dist") {
    dd = as.matrix(dd)
  }
  ## distance matrices must be square
  if (nrow(dd)<1 | ncol(dd)<1) {
    stop("input object must be non-empty\n")
  }
  if (nrow(dd)!=ncol(dd)) {
    stop("input matrix must be square\n")
  }
  ## check that dd is numeric (distances are numeric or integer)
  if (class(dd[,1])!="numeric" & class(dd[,1])!="integer") {
    stop("input object must contain numeric data\n")
  }
  
  ## convert distances into ranks of neighbors
  ans = apply(dd, 2, rank, ties=ties, ...)
  ans = ans/ncol(dd)
  
  ## symmetrize the matrix, make diagonal always zero
  ans = (ans+t(ans))/2
  diag(ans) = 0
  
  ans
}

