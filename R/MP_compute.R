## Package: MultiPattern
##
## Functions for computing based on a MultiPattern configuration
##


## ###################################################################################
## Functions for computing based on a MultiPattern configuration object
## ###################################################################################


##' Run MultiPattern analysis.
##'
##' This function runs analyses specified in a MultiPattern configuration object.
##' It returns an object holding similarity matrices for each of the MultiPattern analyses. 
##'
##' @param MP MultiPattern configuration object
##' @param configs string, vector of configuration names to evaluate;
##' leave NULL to evaluate all configurations
##' @param verbose logical, when TRUE function prints a warning message when the input is large
##'
##' @export
MPgetDistances = function(MP, configs=NULL, verbose=TRUE) {
    
  if (class(MP) != "MultiPattern") {
    stop("object not of class MultiPattern\n")
  }
  
  if (verbose & object.size(MP)>1e6) {
    cat("This may take some time. Please wait... ")
  }
  
  ## check that specified configs are defined
  if (is.null(configs)) {
    configs = names(MP$configs)
  }
  badconfigs = configs[!configs %in% names(MP$configs)]
  if (length(badconfigs)) {
    stop("specified configs are not in MP: ",
         paste(badconfigs, collapse=", "))
  }
  
  ## helper function that computes one distance object given one configuration
  adfun = function(NC) {
    nowdata = MP$data[[NC$data]]
    if (!is.null(NC$prep)) {
      nowdata = NC$prep(nowdata)
    }            
    ## compute similarities and store in the output object
    NC$dist.fun(nowdata)  
  }
  
  ## main part of the work - compute all distances
  ans = lapply(MP$configs[configs], adfun)
  
  ## housekeeping for the list
  names(ans) = configs
  class(ans) = "MultiPatternSimilarities"
  
  if (verbose & object.size(MP)>1e6) {
    cat("done\n")
  }
  
  ans
}




##' Compute a meta distance object from an MP object. Uses subsampling and averaging,
##' which provides a cross-validation of sorts and reduces computation load for large
##' datasets.
##'
##' Settings for subsampling and averaging are determined by
##' MP$settings$subsample.N and MP$settings$subsample.R
##' 
##' @param MP MultiPattern configuration object
##' @param standardize function used to transform similarity matrices
##' @param subsample.N integer, number of samples in each subsampling; when NULL,
##' value is extracted from MP$settings
##' @param subsample.R integer, number of repetitions; when NULL, value
##' is extracted from MP$settings
##' @param alpha numeric, used as an alternative way to transform similarity matrices,
##' applied after transform.fun; when null, value is extracted from MP$settings
##' @param beta numeric, sets the Lbeta norm; when NULL, value is extracted from
##' MP$settings
##' @param verbose logical, set TRUE to print progress messages
##' @param ... additional arguments, used when applying standardize function to the MPS
##' 
##' @export
MPgetAverageMetaDistance = function(MP, standardize=MPrankNeighbors,
                                    subsample.N=NULL, subsample.R=NULL,
                                    alpha=NULL, beta=NULL, verbose=TRUE, ...) {
  
  if (class(MP) != "MultiPattern") {
    stop("object not of class MultiPattern\n")
  }
  
  ## determine properties of sub-sampling and meta-distance scaling
  if (is.null(subsample.N)) {
    subsample.N = MP$settings$subsample.N
  }
  if (is.null(subsample.R)) {
    subsample.R = MP$settings$subsample.R
  }
  if (is.null(alpha)) {
    alpha = MP$settings$alpha
  }
  if (is.null(beta)) {        
    beta = MP$settings$beta
  }
  
  ## check sub sampling is compatible with MP, no. iterations is >= 1
  if (subsample.N<1) {
    subsample.N = abs(length(MP$items)*subsample.N)
  }
  subsample.N = max(1, min(length(MP$items)-1, ceiling(subsample.N)))
  subsample.R = max(1, ceiling(subsample.R))
  
  if (subsample.N<2) {
    stop("Subsampling gives too-small a dataset");
  }
  
  if (verbose) {
    cat("This may take some time. Please wait")
  }
  
  ## define a result object
  result = NULL;
  
  ## compute meta-distances in a loop
  for (i in 1:subsample.R) {
    if (verbose) {
      cat(".")
    }
    
    ## get a subset of the items and create a new MP object
    tempMP = MP
    tempMP$items = sample(MP$items, subsample.N, replace=FALSE)
    ## get smaller versions of the data matrices
    ## (rename rows to avoid duplicate rownames due to bootstrap)
    for (nowdata in names(MP$data)) {
      tempMP$data[[nowdata]] = MP$data[[nowdata]][tempMP$items,]
      rownames(tempMP$data[[nowdata]]) = paste0("B", 1:subsample.N)
    }
    ## compute similarities and meta-similarities for the tempMP
    tempMPS = MPgetDistances(tempMP, verbose=FALSE);
    tempmeta = MPgetMetaDistance(tempMPS, standardize=standardize,
                                 alpha=alpha, beta=beta, ...)
    tempmeta = tempmeta/subsample.R
    
    ## update the average meta-distance result
    if (i==1) {
      result = tempmeta
    } else {
      result = result + tempmeta
    }        
  }
  
  if (verbose) {
    cat(" done\n")
  }
  
  result    
}




##' Compare MultiPattern similarity matrices
##'
##' @param MPS a MultiPatternSimilarities object (list of similarity matrices)
##' @param standardize function used to transform similarity matrices
##' @param alpha numeric, used as an alternative way to transform similarity matrices,
##' applied after transform.fun
##' @param beta numeric, sets the Lbeta norm 
##' @param ... additional arguments, used when applying standardize function to the MPS
##' 
##' @export
MPgetMetaDistance = function(MPS, standardize=MPrankNeighbors, alpha=1, beta=2, ...) {
    
  if (class(MPS) != "MultiPatternSimilarities") {
    stop("object not of class MultiPatternSimilarities\n")
  }
  
  nowsize = length(MPS)
  ans = matrix(0, ncol=nowsize, nrow=nowsize)
  colnames(ans) = rownames(ans) = names(MPS)
  
  ## perhaps standardize/transform the input similarity matrices
  if (!is.null(standardize)) {
    MPS = lapply(MPS, standardize, ...)
  }
  
  ## perhaps compute alpha powers of the input matrix
  if (alpha!=1) {
    MPS = lapply(MPS,
                 function(x) {
                   if (alpha<0) {
                    x2 = 1-(x^alpha)
                   } else {
                     x2 = x^alpha
                   }
                   diag(x2)=0
                   x2[!is.finite(x2)] = 0
                   return(x2)
                 })
  }
  
  ## helper function with L-alpha norm
  Ldist = function(a, b, beta) {
    temp = (abs(a-b))^beta        
    sum(temp)^(1/beta)
  }
  simpledist = function(a, b) {
    sqrt(sum((a-b)*(a-b)))
  }
  
  if (beta==2) {
    for (i in seq_len(nowsize)) {
      for (j in seq_len(nowsize)) {
        if (j>i) {
          ans[i,j] = ans[j,i] = simpledist(MPS[[i]], MPS[[j]])
        }
      }
    }        
  } else {        
    for (i in seq_len(nowsize)) {
      for (j in seq_len(nowsize)) {
        if (j>i) {
          ans[i,j] = ans[j,i] = Ldist(MPS[[i]], MPS[[j]], beta)                    
        }
      }
        }        
  }
  
  ans
}




##' Identify representative points from a distance object
##'
##' @param d dist object or dist-like matrix
##' @param subsets if NULL, function considers all rows
##'     if vector of characters, interpreted as regex for subsetting rows
##' @param method character, determines what clustering approach
##'     is used to determine representative items within a clustering
##' @param k numeric or integer. If more than 1, number of points to select
##'     If less than 1, proportion of points to select
##' 
##' @export
MPgetRepresentatives = function(d, subsets=NULL,
                                method=c("complete", "single", "average", "pam", "extreme"),
                                k=0.05) {
  
    ## the default behavior for subsets is to look at all the items
  if (is.null(subsets)) {
    return(MPgetRepresentatives(d, subsets=".", k=k, method=method)[[1]])
  }
  
  if (class(d)!="dist" && class(d) !="matrix") {
    stop("input must be a dist or matrix")
  }
  dmat = as.matrix(d)    
  if (is.null(rownames(dmat))) {
    stop("input must have names\n")
  }
  method = match.arg(method)
    
  ## helper function to compute central point of a cluster
  ## x - set of item names
  ## dd - matrix of coordinates
  getCenter = function(x, dd) {
    ## compute mean position of the group of points
    midpos = apply(dd[x,,drop=FALSE], 2, median)
    tempdist = dd[x,,drop=FALSE] - matrix(rep(midpos, each=length(x)), nrow=length(x))
    tempdist = apply(tempdist^2, 1, sum)
    return(names(which.min(tempdist)))   
  }
    
  ## helper function estimates points that are far apart from each other (greedy)
  ## xd - distance object
  ## n - number of elements to select
  getExtremes = function(xd, n) {
    ## scale distances to prevent multiple extreme items clumping
    xdm = as.matrix(xd)
    ## find the first point that is furthers from all others
    ans = names(which.max(colSums(xdm)))
    ## scale distances to prevent multiple extreme items clumping
    xdm = xdm^(1/n)
    for (i in seq_len(n-1)) {
      xdm2 = xdm[!rownames(xdm) %in% ans, ]
      xdm2sums = apply(xdm2[, ans, drop=FALSE], 1, sum)
      ans = c(ans, names(which.max(xdm2sums)))
    }
    return(ans)
  }
  
  ## for each subset, get a set of representative items
  ans = list()
  for (nowsubset in subsets) {
    ## extract subset of rows
    temp = grep(nowsubset, rownames(dmat))
    nowdata = dmat[temp, temp, drop=FALSE]
    rm(temp)
    nowdist = as.dist(nowdata)
    nowk = ifelse(k<1, round(k*nrow(nowdata)), round(k))        
    nowk = max(nowk, 1)
    
    ## cluster and extract representative points
    if (nrow(nowdata)>nowk & nrow(nowdata)>1) {
      if (method=="extreme") {
        temp = getExtremes(nowdist, nowk)
      } else {
        ## split the data points into classes accorind to hclust or pam
        if (method=="pam") {
          nowclass = cluster::pam(nowdist, k=nowk, cluster.only=TRUE)
        } else {
          nowclass = cutree(hclust(nowdist, method=method), k=nowk)
        }
        
        temp = cbind(idno=seq(1, nrow(nowdata)), classid=nowclass)
        temp = split(temp[,1], temp[,2])    
        ## at this stage temp has a list with items in each group            
        temp = as.character(unlist(lapply(temp, getCenter, nowdata)))
      }
      
      ans[[nowsubset]] = temp
    }
  }
  
  ans
}
