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
MPrankNeighbors = function(dd, ties="min", ...) {
    
    ## perhapss convert fancy objects into a matrix
    if (class(dd)=="data.frame" | class(dd)=="data.table" | class(dd)=="dist") {
        dd = as.matrix(dd)
    }
    ## distance matrices must be square
    if (nrow(dd)<1 | ncol(dd)<1) {
        stop("input object must be non-empty\n")
    }
    if (nrow(dd)!=ncol(dd)) {
        stop("input matrix must be square\n");
    }
    ## check that dd is numeric (distances are numeric or integer)
    if (class(dd[,1])!="numeric" & class(dd[,1])!="integer") {
        stop("input object must contain numeric data\n")
    }
    
    ## convert distances into ranks of neighbors
    ## There is ambiguity in how to encode, this line defines the encoding
    ans = apply(dd, 2, rank, ties=ties, ...)
    ans = (-1+ans)/(ncol(dd)-1)

    ## symmetrize the matrix
    ans = (ans+t(ans))/2
    
    ## make diagonal entries always zero
    diag(ans) = 0
    
    ans
}





##' Transform a dissimilariy matrix using ranks
##'
##' This function takes as input a distance object or distance matrix (values
##' describing distance between pairs of data elements). For each data element,
##' it considers what other points are nearby and ranks them by distance.
##'
##' 
##' @param dd matrix of distances between neighbors
##' @param ties.method character, used to determine how equal ranks are handled by rank function.
##' @param ... additional arguments passed on to rank 
##'
##' @export
MPrankWhole = function(dd, ties.method="min", ...) {
    
    ## perhapss convert fancy objects into a matrix
    if (class(dd)=="data.frame" | class(dd)=="data.table" | class(dd)=="dist") {
        dd = as.matrix(dd)
    }
    ## distance matrices must be square
    if (nrow(dd)<1 | ncol(dd)<1) {
        stop("input object must be non-empty\n")
    }
    if (nrow(dd)!=ncol(dd)) {
        stop("input matrix must be square\n");
    }
    ## check that dd is numeric (distances are numeric or integer)
    if (class(dd[,1])!="numeric" & class(dd[,1])!="integer") {
        stop("input object must contain numeric data\n")
    }
    
    ## convert distances into ranks of neighbors
    ## There is ambiguity in how to encode, this line defines the encoding
    ans = matrix(rank(dd, ties.method=ties.method, ...), nrow=nrow(dd), ncol=ncol(dd))
    ans = (-1+ans)/ncol(dd)
    rownames(ans) = rownames(dd)
    colnames(ans) = colnames(dd)        
    diag(ans) = 0
    
    ans
}





##' Get a set of nearest neighbors (ranked)
##'
##' This function scores graph similarity using rank neighborhoods
##' 
##' @param dd - distance object
##' @param seedsample - name of sample
##' @param maxrank - integer/numeric. Leave 0 to obtain complete ranking. Set to N if
##' interested only in items with rank <=N
##' @param ... - passed on to rank()
##'
##' @return vector with names, values are ranks, names correspond to elements in dd
##' 
##' @export 
MPgetNeighborSet = function(dd, seedsample, maxrank=0, ...) {

    if (class(dd)=="dist") {
        dd = as.matrix(dd)
    }
    
    distto = dd[seedsample, colnames(dd)!=seedsample]
    distto = sort(rank(distto, ...))
    if (maxrank >= 1) {
        distto = distto[distto<= maxrank]
    }
    
    return(distto)    
}







