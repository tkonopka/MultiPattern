## Package: MultiPattern
##
## A set of helper functions used in MultiMetric
##




##' obtain a list of all pairs from a set
##' 
##' @param f vector of features (or numbers for column names)
##'
##' @export
MPgetAll2Subspaces = function(f) {    
  flen = length(f);
  if (flen<2) {
    stop("the feature vector must be at least of length 2\n");
  }
  ans = matrix(NA, ncol=2, nrow=flen*(flen-1)/2)
  kk = 1        
  for (i in 1:(flen-1)) {
    for (j in (i+1):flen) {
      ans[kk,] = c(f[i], f[j])
      kk = kk+1
    }
  }
  ans = split(ans, seq(1, nrow(ans)))
  names(ans) = sapply(ans, paste, collapse=".")    
  ans
}




##' Produce a list of randomly selected feature sets
##'
##' Note: this function is not guaranteed to produce the specified
##' number of features (it might produce less).
##'
##' For example, if a feature set is of length 3 and user requets 5 sets of
##' size 3, the output will contain only 1 (because there are not others).
##'
##' If a feature set is of length 6 and user requests 20 sets of length 3,
##' the output might contain fewer than 20. Although in this case there are exactly 20
##' sets of length 3, the randomization might not find them all. In such cases, try
##' increasing the setting oversample.
##' 
##' @param f character or integer vector, feature set
##' @param n integer, number of required feature sets
##' @param d integer, dimension of the required subspaces
##' @param oversample numeric, used in randomization process (only becomes relevant
##' when function starts to return fewer subspaces than requested, see description above)
##' 
##' @export
MPgetRandomSubspaces = function(f, n, d=4, oversample=4) {
    
  if (length(f)<d) {
    stop("Feature vector is shorter than requested subspace")
  }
  if (n<=0) {
    return(list())
  }
  
  ## try to generate many sets of features
  n2 = ceiling(n*oversample)

  ans = matrix(NA, nrow=n2, ncol=d)
  
  ## create unique combinations of sets
  for (i in seq_len(nrow(ans))) {
    ans[i,] = sort(sample(f, d, replace=FALSE))
  }
  ans = unique(ans)
  
  ## cut the subspaces to the requested n
  if (nrow(ans)>n) {
    ans = ans[1:n,]
  }
  
  ## turn the matrix into a list object with arbitrary names for the features
  split(ans, 1:nrow(ans))
}




##' Make a suggestion for interesting features in a dataset based on PCA weights
##'
##' This function looks at the features that make up individual PCA components. It
##' returns those features that contribute more than others to these componets
##' (according to the weight of the feature compared to the median and IQR interval)
##' 
##' @param x dataset or prcomp object
##' @param ncomp function, integer, or numeric.
##'              If integer, number of principal components to consider.
##'              If number less than 1, fraction of components to consider.
##'              If function, a transformation on the number of features that returns an integer.
##' @param iqrfactor numeric. threshold to consider a feature in a PC as worthy of selection
##' 
##' @export
MPsuggestTopFeatures = function(x, ncomp=function(y) { 1+sqrt(y) },
                                iqrfactor=function(y) { 0.5+log(y)} ) {
  
  ## transform the data using PCA
  if (class(x)!="prcomp") {
    ## eliminate constant columns in x (constant columns make prcomp non-deterministic)
    x = x[, apply(x, 2, sd)>0, drop=FALSE]
    if (ncol(x)<2) {
      return(colnames(x))
    }
    x = prcomp(x, scale.=FALSE)
  }
  
  ## convert the ncomp into an integer
  if (class(ncomp)=="function") {
    ncomp = max(1, ceiling(ncomp(ncol(x$rotation))))
  } else {
    if (as.integer(ncomp)<1) {
      ncomp = max(1, ceiling(ncol(x$rotation)*ncomp))
    }
    ncomp = ceiling(ncomp)
  }    
  ncomp = max(1, min(ncomp, ncol(x$rotation)))
  
  ## convert the iqrfactor into a numeric value
  if (class(iqrfactor)=="function") {
    iqrfactor = iqrfactor(nrow(x$rotation))
  }

  ## look at the top ncomp components and find features that contribute more than others
  ans = list()
  for (i in 1:ncomp) {
    temp = x$rotation[,i]
    temp.qs = as.numeric(quantile(abs(temp), p=c(0.25, 0.5, 0.75)))
    temp.iqr = temp.qs[3]-temp.qs[1]
    ## get at least one feature
    temp.max = (abs(temp)==max(abs(temp)))
    ## get features that contribute more than the rest
    temp.iqr =
      (temp > temp.qs[2] + (iqrfactor*temp.iqr)) |
      (temp < temp.qs[2]-(iqrfactor*temp.iqr))
    ans[[i]] =sort( temp[temp.max | temp.iqr])
  }
  
  ## format the output answer into a vector 
  ans = sort(unique(names(unlist(ans))))
  
  ans
}




##' Create a function for preprocessing a dataset before computing a similarity
##'
##' @param a various types of inputs are accepted.
##' If a is a function, output will be exactly a.
##' If a is a character vector, output is a function that returns the a-columns in a matrix.
##' If a is an integer vector
##' 
##' @export
MPmakePrepFunction = function(a) {
  
  ## no preprocessing or applying a function
  if (class(a)=="function" | is.null(a)) {
    return(a)
  }
  
  ## other output types use a to change a matrix. Need to remember what a is
  force(a)
  
  ## subsetting by named columns in a matrix, or order of columns in a matrix
  if (is.character(a) | (is.integer(a) & is.null(names(a)))) {
    return(function(x) x[, a, drop=FALSE])        
  }
  
  ## subsetting by named columns with rescaling factors
  if (is.numeric(a) & !is.null(names(a))) {
    return(function(x) {
      rescale = matrix(a, nrow=nrow(x), ncol=length(a), byrow=TRUE)
      x[, names(a), drop=FALSE]*rescale
    })        
  }
  
  ## if reached here, unrecognized type of prep
  stop("Unrecognized type of prep instruction:",
       " class:", class(a), ", names:", !is.null(names(a)), "\n")
}




## compute an appropriate number of components for ica, pca
##
## n.col - integer, number of columns (features) in data matrix
## n.comp - desired number of components (works with n.comp<1)
##
get.n.comp = function(n.col, n.comp) {
  if (n.col < 2 | n.comp<=0) {
    return (NULL)
  }
  if (n.comp<=1) {
    n.comp = n.comp*n.col
  }
  max(2, min(ceiling(n.comp), n.col))
}



## helper function cuts a matrix into a small with at least 2 columns
## and with all columns that contain covthreshold of the total variance
##
## dd - input matrix
## n.comp - number of PCA features
##          if less than one, a proportion of features from dd
##
## outputs a new matrix with the same number of rows as dd, with fewer columns
##
getPCAsubset = function(dd, n.comp=4) {
  n.comp = get.n.comp(ncol(dd), n.comp)
  if (is.null(n.comp)) {
    return(NULL)
  }
  result = rsvd::rpca(dd, n.comp)$x
  colnames(result) = paste0("PC", seq(ncol(result)))
  rownames(result) = rownames(dd)
  result
}




## transform matrix using ICA and obtain a matrix with a
## small number of components
##
## outputs a new matrix with same number of rows as dd, with fewer columns
##
getICAsubset = function(dd, n.comp=2) {
  n.comp = get.n.comp(ncol(dd), n.comp)
  if (is.null(n.comp)) {
    return(NULL)
  }
  ## There is PCA pre-processing here.
  ## Somehow, that makes fastICA itself much faster and use less memeory
  ddpca = rsvd::rpca(dd)$x
  ddica = fastICA::fastICA(ddpca, n.comp, method="C")
  result = ddica$S[, 1:n.comp]
  colnames(result) = paste0("IC", seq(ncol(result)))
  rownames(result) = rownames(dd)
  result
}




##' Create equal size xlim and ylim intervals from x,y values
##'
##' This function is useful for creating MDS diagrams in a square format,
##' which should have x and y axes in the same scale.
##'
##' @param x values for the x axis
##' @param y values for the y axis
##'
##' @export
MPsquarelim = function(x, y) {
  ## get plain xlim/ylim
  xlim = range(x)
  ylim = range(y)
  xsize = xlim[2]-xlim[1]
  ysize = ylim[2]-ylim[1]    
  ## make the limits "square"
  if (xsize>ysize) {
    ylim = ylim + (xsize-ysize)*c(-0.5,0.5)
  } else {
    xlim = xlim - (xsize-ysize)*c(-0.5,0.5)
  }    
  list(xlim=xlim, ylim=ylim)
}




##' Create a map for a distance object/matrix
##'
##' This function uses cmdscale, or tsne.
##' 
##' @param d distance object
##' @param method character, method to prepare map
##' @param ... other parameters can be passed to umap
##'
##' @export
MPgetMap = function(d, method=c("umap", "cmdscale"), ...) {
  
  method = match.arg(method)
  
  ## get a cmdscale map as an initial configuration
  d = as.matrix(d)
  
  # perhaps adjust the cmdscale result using tsne
  if (method=="umap") {
    conf = umap::umap.defaults
    conf$init = "spectral"
    conf$input = "dist"
    conf$min_dist = 0.2
    V = nrow(d)
    conf$n_neighbors = min(V-1, 1+floor(sqrt(V)))
    ans = suppressWarnings(umap::umap(d, conf, ...))
    ans = ans$layout
  } else if (method=="cmdscale") {
    ans = cmdscale(d)
  }
  
  ans
}




##' Permute values in a matrix
##'
##' @param x input matrix
##' @param perm.method character, determines how the values are permuted
##'
##' The randomization methods are as follows.
##' shuffle - re-arranges values in a matrix using sampling *without* replacement
##' bootstrap - re-arranges values in a matrix using sampleing *with* replacement
##' 
##' @export
MPrandomizeMatrix = function(x, perm.method=c("shuffle", "bootstrap")) {
    
  perm.method = match.arg(perm.method)
  
  ## get a vector of the x values
  xv = as.numeric(x)
  xp = NULL
  
  if (perm.method=="shuffle") {        
    xp = matrix(sample(xv, length(xv), replace=F), ncol=ncol(x))
  } else if (perm.method=="bootstrap") {
    xp = matrix(sample(xv, length(xv), replace=T), ncol=ncol(x))
  }
  
  rownames(xp) = rownames(x)
  colnames(xp) = colnames(x)
  
  xp
}




##' Normalize a matrix column-wise by median shift and division by sd
##'
##' @param x numeric matrix or data.frame
##' 
##' @export
MPmatrixColNorm = function(x) {    
  xp = x
  for (i in 1:ncol(x)) {
    medx = median(x[,i])
    sdx = sd(x[,i])
    xp[,i] = (x[,i]-medx)/sdx            
  }
  ## avoid NA or Inf when sdx is zero
  xp[!is.finite(xp) | is.na(xp)] = 0
  xp        
}




##' Provides a cut of hclust into k clusters
##'
##' In contrast to stats::cutree(), this can output more than k clusters. However, the
##' output should contain k clusters that are of nontrivial size. The largest cluster
##' should be labeled as "1", the next biggest as "2" and so on
##'
##' @param tree tree object output by hclust
##' @param k integer, number of clusters to output
##' @param minsize numeric in range [0, 1], minimal size of the clusters
##'
##' @export
MPcutree = function(tree, k=2, minsize=1/(10*k)) {
    
  ## size of the tree
  tree.n = length(tree$order)
  ## required cluster size
  csize = tree.n*minsize
  ## attempted size of cluster
  kattempt = k
  
  ## strategy: try cutting the tree with small ks, until reach a point where
  ## the cuts produce at least k clusters of size csize
  while (kattempt < tree.n/k) {
    nowcut = cutree(tree, k=kattempt)
    nowsizes = table(nowcut)
    if (sum(nowsizes>=csize)>=k) {
      break;
    }
    kattempt = kattempt+1        
  }
  
  ## repeat the cut (this is necessary when k is set too large at the beginning)
  nowcut = cutree(tree, k=kattempt)
  nowsizes = table(nowcut)
  
  ## rename the clusters so that large clusters have low ids, cluster 0 contains outliers
  nowrank = rank(-nowsizes, ties.method="first")
  nowcut2 = nowrank[nowcut]
  names(nowcut2) = names(nowcut)
  nowcut2[nowcut2>k] = 0
  
  nowcut2
}

