## Package: MultiPattern
##
## Collection of "distance" functions
##
## Some of these functions produce "distance" with metric properties.
## Others produce dissimilarities without following metric properties
##
## Many are just wrappers special cases of functions defined elsewhere
##




##' Impute values for NAs in a matrix or other object
##'
##' (This is a general-purpose function and works on many types of objects with some NAs.
##' The object must support operations: sum, method (e.g. mean), and subsetting via []).
##' 
##' @param x matrix or dist object
##' @param method function, is applied on x to find the imputed values. This
##' function must be able to perform method(x, na.rm=TRUE)
##'
##' @export
MPimputeNAs = function(x, method=mean) {    
  ## handle na values
  if (sum(is.na(x))>0) {
    impute.val = method(x, na.rm=TRUE)
    if (is.na(impute.val)) {
      impute.val = 0
    }
    x[is.na(x)] = impute.val
  }
  ## handle infinite values
  if (sum(!is.finite(x))>0) {
    max.val = max(x[is.finite(x)])
    min.val = min(x[is.finite(x)])
    if (!is.finite(max.val)) {
      max.val = 1
    }
    if (!is.finite(min.val)) {
      min.val = 0
    }
    x[!is.finite(x) & x>0] = max.val
    x[!is.finite(x) & x<0] = min.val
  }
  
  x
}




##' Impute values for NAs in a dist object
##'
##' (This is a specialized function - it converts x into a dist.
##' This is to enforce consistency when applying the 'method' function.
##' e.g. to avoid confusion between mean(matrix) and mean(as.dist(matrix)),
##' which are not the same)
##' 
##' @param x matrix or dist object
##' @param method function, is applied on x to find the imputed values. This
##' function must be able to perform method(x, na.rm=TRUE)
##'
##' @export
MPdistImpute = function(x, method=mean) {
  if (!"dist" %in% class(x)) {
    x = as.dist(x)
  }
  MPimputeNAs(x, method=method)
}




##' Euclidean distance
##'
##' This is a wrapper for stats::dist with method="euclidean".
##' The wrapper replaces NA and Inf values with mean and max guesses.
##'
##' @param x numeric matrix (see stats::dist for details)
##' 
##' @export
dist.euclidean = function(x) {
  MPdistImpute( stats::dist(x, method="euclidean") )
}




##' Canberra distance
##'
##' This is a wrapper for stats::dist with method="canberra".
##' The wrapper replaces NA and Inf values with mean and max guesses.
##' 
##' @param x numeric matrix
##'
##' @export
dist.canberra = function(x) {
  MPdistImpute( stats::dist(x, method="canberra") )
}




##' Hamming distance
##'
##' This counts the number of columns in x that are different between rows.
##' It is a suitable distance function for factor and character data.
##'
##' @param x matrix with numbers, characters, or factors
##'
##' @export
dist.hamming = function(x) {  
  nn = nrow(x)
  ans = matrix(0, ncol=nn, nrow=nn)
  for (i in 1:(nn-1)) {
    ix = x[i,]
    for (j in (i+1):nn) {
      iy = x[j,]
      ans[i,j] = sum(ix!=iy)
    }
  }
  ans = ans + t(ans)
  rownames(ans) = colnames(ans) = rownames(x)
  
  MPdistImpute(stats::as.dist(ans))
}




##' Manhattan distance
##'
##' This is a wrapper for stats::dist with method="manhattan".
##' The wrapper replaces NA and Inf values with mean and max guesses.
##'
##' @param x numeric matrix
##'
##' @export
dist.manhattan = function(x) {
  MPdistImpute( stats::dist(x, method="manhattan") )
}




##' Distance using 1st principal component from PCA
##'
##' The function computes the PCA, then collapses onto PC1.
##' 
##' @param x numeric matrix
##'
##' @export
dist.pc1 = function(x) {
  ## just in case, fill in NA values in x
  x = MPimputeNAs(x)
  x1 = prcomp(x)$x[,1, drop=FALSE]
  stats::dist(x1, method="euclidean")    
}




##' Compute dissimilarity using correlation
##'
##' @param x matrix of values
##' @param directional logical, set TRUE to get dissimilarities as (1-r) and FALSE
##' to compute dissimilarities as (1-abs(r))
##' @param use.ranks logical, set TRUE to use rho correlation for ranks
##'
##' @export
dist.correlation = function(x, directional=TRUE, use.ranks=FALSE) {

  ## helper: computing correlation using centered values
  centered.corr = function (a, b) {
    sum(a*b) / sqrt(sum(a*a)*sum(b*b))
  }

  nn = nrow(x)
  xrownames = rownames(x)
  
  ## transpose the matrix so most computations will be done by column
  x = t(x)
  if (use.ranks) {
    x = apply(x, 2, rank)
  }
  
  ## center the matrix
  x = sweep(x, 2, apply(x, 2, mean))
  
  ## compute correlations
  ans = matrix(0, ncol=nn, nrow=nn)
  for (i in 1:(nn-1)) {
    ix = x[,i]
    for (j in (i+1):nn) {
      iy = x[,j]
      ans[i,j] = centered.corr(ix, iy)
    }
  }
  ans = ans + t(ans)
  rownames(ans) = colnames(ans) = xrownames

  ## apply either directional or non-directional dissimilarity
  if (directional) {
    ans = 1-ans
  } else {
    ans = 1-abs(ans)
  }  
  diag(ans) = 0

  MPdistImpute(ans)
}




##' Compute distance using pearson correlation r
##' 
##' @param x matrix of values
##' @param directional logical, set TRUE to get dissimilarities as (1-r) and FALSE
##' to compute dissimilarities as (1-abs(r))
##' 
##' @export
dist.pearson = function(x, directional=TRUE) {
  dist.correlation(x, directional=directional, use.ranks=FALSE)
}





##' Compute distance using spearman rho
##' 
##' @param x matrix of values
##' @param directional logical, set TRUE to get dissimilarities as (1-r) and FALSE
##' to compute dissimilarities as (1-abs(r))
##' 
##' @export
dist.spearman = function(x, directional=TRUE) {
  dist.correlation(x, directional=directional, use.ranks=TRUE)
}




##' Euclidean distance of log2 transformed values
##'
##' @param x numeric matrix
##' @param shift numeric, function computes log2(x+shift) before
##' applying euclidean distance
##'
##' @export
dist.log2euclidean = function(x, shift=1) {    
  MPdistImpute( stats::dist(log2(x+shift), method="euclidean") )
}




##' Distance that incorporates weights from a clustering (hclust or pam)
##'
##' @param x numeric matrix
##' @param clust.dist character, determines what preliminary distance function is
##' applied prior to clustering
##' @param clust.method character, determines agglomeration method in hclust (e.g. complete,
##' single, average) or "pam"
##' @param clust.k integer, determines depth of clustering weights
##' @param clust.weight numeric, determines weighting of cluster distances relative to
##' neighbor rank distance.
##' Set to 0 to obtain pure neighbor rank distance (clustering is ignored).
##' Set to Inf to obtain pure cluster-based distance.
##; Default is 0.5 to nudge neighbor rank distance in the direction of the cluster distance.
##' @param clust.alt boolean. Set TRUE to obtain an alternative clustering. Default is FALSE
##' which returns distances in which patterns from a usual distance are reinforced by clustering.
##' 
##' @export
dist.clust = function(x,
                      clust.dist="euclidean", clust.method="complete", 
                      clust.k=2, clust.weight=0.5, clust.alt=FALSE) {
  
  ## this function needs sample identifiers, so create some if don't exist
  xdm = x
  if (is.null(rownames(x))) {
    rownames(xdm) = paste0("S", 1:nrow(x))
  }
  
  if (clust.k*2 > nrow(x)+2) {
    stop("dist.clust must act on object with a minimum number of objects clust.k*2")
  }
  
  ## produce a first distance matrix using some method, and make a clustering
  xd = MPdistFactory(clust.dist)(xdm)
  xdm = MPrankNeighbors(as.matrix(xd))
  if (clust.method=="pam") {
    ##require("cluster")
    ## pam clustering does not guarantee that k groups consist of 2k groups merged
    ## together in some pattern. So here obtain 2k groups first and then manually merge
    ## them pairwise
    xpam = cluster::pam(xd, k=clust.k*2, cluster.only=FALSE)
    xcut2 = xpam$clustering
    xcut2 = split(names(xcut2), xpam$medoids[xcut2])
    ## get distances between medoids
    xdmed = as.dist(as.matrix(xd)[xpam$medoids, xpam$medoids])
    xmedpam = cluster::pam(xdmed, k=clust.k, cluster.only=TRUE)
    xcut1 = split(names(xmedpam), xmedpam)
    xcut1 = lapply(xcut1,
                   function(x) {
                     temp = unlist(xcut2[x])
                     names(temp) = NULL
                     temp
                   })        
  } else {
    xdh = hclust(xd, method=clust.method)
    ## produce conventional cuts in the tree at level k
    xcut1 = MPcutree(xdh, k=clust.k)
    xcut2 = MPcutree(xdh, k=clust.k*2)
    xcut1 = split(names(xcut1), xcut1)
    xcut2 = split(names(xcut2), xcut2)
  }
  
  ## adjust the neighbor adjusted distances according to the cutree classes
  if (clust.alt) {
    xdadd = matrix(0, ncol=ncol(xdm), nrow=nrow(xdm),
                   dimnames=list(colnames(xdm), colnames(xdm)))
    ## first penalize within cluster at level clust.k
    for (i in seq_along(xcut1)) {
      temp = xcut1[[i]]
      xdadd[temp, temp] = 1
    }
    ## then un-penalize items within clusters at level clust.k*2
    for (i in seq_along(xcut2)) {
      temp = xcut2[[i]]
      xdadd[temp, temp] = 0
    }        
  } else {    
    ## create a dissimilarity matrix where all items are equally distant
    xdadd = matrix(1, ncol=nrow(xdm), nrow=nrow(xdm),
                   dimnames=list(rownames(xdm), rownames(xdm)))
    
    for (i in seq_along(xcut1)) {
      temp = xcut1[[i]]
      xdadd[temp, temp] = xdadd[temp, temp] - 1; 
    }
    diag(xdadd) = 0                
  }
  
  if (is.finite(clust.weight)) {
    xdm = xdm+(clust.weight*xdadd)
  } else {
    xdm = xdadd
  }
  
  if (is.null(rownames(x))) {
    colnames(xdm) = rownames(xdm) = NULL
  }
  
  ## return a dist object 
  MPdistImpute(xdm)
}



##' Distance that incorporates weights from a dbscan clustering 
##'
##' @param x numeric matrix
##' @param eps numeric, passed to dbscane eps
##' @param clust.weight numeric, determines weighting of cluster distances relative to
##' neighbor rank distance.
##' Set to 0 to obtain pure neighbor rank distance (clustering is ignored).
##' Set to Inf to obtain pure cluster-based distance.
##' Default is 0.5 to nudge neighbor rank distance in the direction of the cluster distance.
##' 
##' @export
dist.dbscan = function(x, eps=1, clust.weight=0.8) {

  x.data = x
  if (class(x.data)!="matrix") {
    x.data = as.matrix(x.data)
  }
  x.data = MPimputeNAs(x.data)
  if (is.null(rownames(x))) {
    rownames(x.data) = paste0("S", 1:nrow(x))
  }

  ## produce a first distance matrix
  xd = dist.euclidean(x.data)
  xdm = MPrankNeighbors(as.matrix(xd))

  ## produce a clustering
  x.dbscan = dbscan::dbscan(x.data, eps=eps)
  xcut = split(rownames(x.data), x.dbscan$cluster)
  
  ## create a dissimilarity matrix where all items are equally distant
  ## then adjust making same-cluster items closer
  xdadd = matrix(1, ncol=nrow(xdm), nrow=nrow(xdm),
                 dimnames=list(rownames(xdm), rownames(xdm)))
  for (i in seq_along(xcut)) {
    if (names(xcut)[i] != "0") {
      temp = xcut[[i]]
      xdadd[temp, temp] = xdadd[temp, temp] - 1; 
    }
  }
  diag(xdadd) = 0                
  
  if (is.finite(clust.weight)) {
    xdm = xdm+(clust.weight*xdadd)
  } else {
    xdm = xdadd
  }
  if (is.null(rownames(x))) {
    colnames(xdm) = rownames(xdm) = NULL
  }
  
  ## return a dist object
  MPdistImpute(xdm)
}




##' Generates a function with syntax f(x) that returns a dist object
##'
##' @param method character, one of the available options. The output function
##' will be able to compute dist using the specified method
##' @param directional logical, used with method=pearson and method=spearman
##' @param log2.shift numeric, Used with method=log2euclidean.
##' @param clust.dist character, used with method=hclust
##' @param clust.method character, used with method=hclust
##' @param clust.weight numeric, used with method=hclust
##' @param clust.k integer, used with method=hclust
##' @param clust.alt boolean, used with method=hclust
##' @param eps numeric used with method=dbscane
##'
##' @export
MPdistFactory = function(method=c("euclidean", "manhattan", "canberra", "pc1",
                                  "log2euclidean", "pearson", "spearman",
                                  "hclust", "pam", "dbscan"),
                         directional = TRUE,
                         log2.shift = 1,   
                         clust.dist="euclidean",
                         clust.method=c("complete", "single", "average"),
                         clust.weight=0.8, clust.k=2, clust.alt=FALSE,
                         eps=1) {
  
  ## collapse the given method into one of the allowed options
  method = match.arg(method)
  
  ## process the available options one by one
  ## (Perhaps switch would be easier, but 
  if (method=="euclidean") {
    return(dist.euclidean)
  } else if (method=="canberra") {
    return(dist.canberra)
  } else if (method=="manhattan") {
    return(dist.manhattan)
  } else if (method=="pc1") {
    return(dist.pc1)
  } else if (method=="pearson") {
    force(directional)
    return(function(x) { dist.pearson(x, directional=directional) })
  } else if (method=="spearman") {
    force(directional)
    return(function(x) { dist.spearman(x, directional=directional) })
  } else if (method=="log2euclidean") {
    force(log2.shift)
    return(function(x) {
      dist.log2euclidean(x, shift=log2.shift)
    })
  }

  ## the following options are sligtly more complex distance functions
  ## that require a bit of extra code    
  if (method=="hclust" | method=="pam") {
    if (method=="pam") {
      clust.method = "pam"
    } 
    clust.method = clust.method[1]
    force(clust.weight)
    force(clust.k)
    force(clust.dist)
    force(clust.method)
    force(clust.alt)
    return(
      function(x) {                
        dist.clust(x, clust.dist=clust.dist, clust.method=clust.method,
                   clust.weight=clust.weight, clust.k=clust.k,
                   clust.alt=clust.alt)
      })
  }

  if (method=="dbscan") {
    force(eps)
    force(clust.weight)
    return(
      function(x) {
        dist.dbscan(x, eps=eps, clust.weight=clust.weight)
      })
  }
  
}

