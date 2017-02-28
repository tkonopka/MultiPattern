## Collection of "distance" functions
##
## Some of these functions produce "distance" with metric properties (metric axioms).
## Others produce "distance" that are informative as similarities without following metric properties
##
## Many are just wrappers special cases of functions defined elsewhere
##
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
        if (!is.finite(max.val)) {
            max.val = 1
        }
        x[!is.finite(x)] = max.val
    }
    return (x)
}





##' Impute values for NAs in a dist object
##'
##' (This is a specialized function - it converts x into a dist as the first step.
##' This is to enforce consistency when applying the 'method' function.
##' e.g. to avoid confusion between mean(matrix) and mean(as.dist(matrix)), which are not the same)
##' 
##' @param x matrix or dist object
##' @param method function, is applied on x to find the imputed values. This
##' function must be able to perform method(x, na.rm=TRUE)
##'
##' @export
MPdistImpute = function(x, method=mean) {
    x = as.dist(x)
    MPimputeNAs(as.dist(x), method=method)
}




##' Random distance 
##'
##' This produces a distances object in which distances are set to random numbers
##' from the N(0,1) distribution
##'
##' @param x numeric matrix (see stats::dist for details) The object is used to obtain the
##' size of the random distance object
##'
##' @export
dist.rnorm = function(x) {    
    NN = nrow(x)
    ans = matrix(abs(rnorm(NN*NN)), ncol=NN, nrow=NN)
    ans[lower.tri(ans)] = t(ans)[lower.tri(ans)]
    diag(ans) = 0
    rownames(ans) = colnames(ans) = rownames(x)
    return(as.dist(ans))        
}




##' Random distance using beta distribution
##'
##' This function generates points in a d-dimensional unit-space with coordinates
##' chosen at random from a beta distribution. The function then computes the euclidean
##' distance in this space.
##' 
##' @param x numeric matrix (see stats::dist for details) The object is used to obtain the
##' size of the random distance object
##' @param alpha numeric, shape parameter for beta distribution
##' @param beta numeric, shape parameter for beta distribution
##' 
##' @export
dist.rbeta = function(x, alpha=0.5, beta=0.5) {
    ## generate random points using a beta distribution
    NN = nrow(x)
    ans = matrix(abs(rbeta(NN*NN, alpha, beta)), ncol=NN, nrow=NN)
    ans[lower.tri(ans)] = t(ans)[lower.tri(ans)]
    diag(ans) = 0
    rownames(ans) = rownames(x)
    return(as.dist(ans))        
}




##' Random distance using a normal distribution
##'
##' This function generates points in a d-dimensional unit-space with coordinates
##' chosen at random from a normal distribution. The function then computes the euclidean
##' distance in this space.
##' 
##' @param x numeric matrix (see stats::dist for details) The object is used to obtain the
##' size of the random distance object
##' @param d integer, dimension of the underlying space
##' 
##' @export
dist.rnorm0 = function(x, d=round(log2(nrow(x)))) {
    ## make sure d is at least 1, otherwise the matrix ans will be all empty
    d = round(min(1, d))
    ## generate random points using a beta distribution
    NN = nrow(x)
    ans = matrix(rnorm(NN*d), ncol=d, nrow=NN)
    rownames(ans) = rownames(x)
    return(dist(ans))        
}




##' Random distance obtained through permutation of an input matrix
##'
##' @param x numeric matrix
##' @param dist.method character, determines how the distance will be computed,
##' e.g. euclidean, canberra, etc.
##' @param perm.method character, determines how the permutation is carried out,
##' e.g. free, marginals
##' 
##' @export
dist.perm = function(x, dist.method=c("euclidean", "canberra", "spearman"),
    perm.method=c("shuffle", "bootstrap")) {
    
    ## get a function to compute distance
    dist.method = match.arg(dist.method)    
    dist.fun = MPdistFactory(method=dist.method)
    
    ## permute the input data
    perm.method = match.arg(perm.method)
    xp = MPrandomizeMatrix(x, perm.method=perm.method)
    
    return(dist.fun(xp))
}





##' Euclidean distance
##'
##' This is a wrapper for stats::dist with method="euclidean".
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
##' 
##' @param x numeric matrix
##'
##' @export
dist.canberra = function(x) {
    MPdistImpute( stats::dist(x, method="canberra") )
}




##' Manhattan distance
##'
##' This is a wrapper for stats::dist with method="manhattan". 
##'
##' @param x numeric matrix
##'
##' @export
dist.manhattan = function(x) {
    MPdistImpute( stats::dist(x, method="manhattan") )
}




##' Distance using 1st principal component from PCA
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




##' Compute distance using spearman rho
##'
##' This function rows as features and columns as features. 
##' 
##' @param x - matrix of values
##' 
##' @export
dist.spearman = function(x) {

    x = t(x)
    
    ## helper function for computing spearman rho from ranks
    rhoFromRanks = function (xranks, yranks, xmid, ymid) {
        TA = sum((xranks - xmid) * (yranks - ymid))
        TB = sqrt(sum((xranks - xmid)^2) * sum((yranks - ymid)^2))
        return(TA/TB)
    }

    nn = ncol(x)
    
    ## precompute ranks for all columns
    xranks = apply(x, 2, rank)
    rankmeans = apply(xranks, 2, mean)
    
    ## make a matrix with distances
    ans = matrix(0, ncol=nn, nrow=nn)
    for (i in 1:(nn-1)) {
        for (j in (i+1):nn) {
            nowdist = 1-abs(rhoFromRanks(xranks[,i], xranks[,j], rankmeans[i], rankmeans[j]))
            ans[i,j] = nowdist
        }
    }
    ans[lower.tri(ans)] = t(ans)[lower.tri(ans)]
    rownames(ans) = colnames(ans) = colnames(x)
     
    return(MPdistImpute(ans))
}




##' Euclidean distance of log2 transformed values
##'
##' @param x numeric matrix
##' @param shift numeric, function computes log2(x+shift) before applying euclidean distance
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
##' neighbor rank distance. Set to 0 to obtain pure neighbor rank distance (clustering is ignored).
##' Set to Inf to obtain pure cluster-based distance. Default is 0.5 to nudge neighbor rank distance
##' in the direction of the cluster distance.
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
        stop("dist.k2alt must act on object with a minimum number of objects clust.k*2")
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
        xdadd = matrix(0, ncol=ncol(xdm), nrow=nrow(xdm), dimnames=list(colnames(xdm), colnames(xdm)))
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
    
    ## return a dist object 
    if (is.null(rownames(x))) {
        colnames(xdm) = rownames(xdm) = NULL
    }

    ## return a dist object 
    return(MPdistImpute(xdm))
}







##' Generates a function with syntax f(x) that returns a dist object
##'
##' @param method character, one of the available options. The output function
##' will be able to compute dist using the specified method
##' @param log2.shift numeric. Used with method=log2euclidean.
##' @param clust.dist character, used with method=clust
##' @param clust.method character, used with method=clust
##' @param clust.weight numeric, used with method=clust
##' @param clust.k integer, used with method=clust
##' @param clust.alt boolean, used with method=clust
##' 
##' @export
MPdistFactory = function(
    method=c("euclidean", "manhattan", "canberra", "rnorm", "pc1",
        "log2euclidean", "spearman", "hclust", "pam"),
    log2.shift = 1,   
    clust.dist="euclidean", clust.method=c("complete", "single", "average"),
    clust.weight=0.5, clust.k=2, clust.alt=FALSE
    ) {
    
    ## collapse the given method into one of the allowed options
    method = match.arg(method)
    
    ## process the available options one by one
    ## (Perhaps switch would be easier, but 
    if (method=="euclidean") {
        return(dist.euclidean)
    }        
    if (method=="canberra") {
        return(dist.canberra)
    }
    if (method=="manhattan") {
        return(dist.manhattan)
    }
    if (method=="pc1") {
        return(dist.pc1)
    }    
    if (method=="rnorm") {
        return(dist.rnorm)
    }    
    if (method=="spearman") {
        return(dist.spearman)
    }
    if (method=="log2euclidean") {
        force(log2.shift)
        return(
            function(x) {
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
                           clust.weight=clust.weight, clust.k=clust.k, clust.alt=clust.alt)
            })
    }
    
}


