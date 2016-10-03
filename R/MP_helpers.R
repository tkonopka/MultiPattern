## A set of helper functions used in MultiMetric
##
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
    ans = matrix("", ncol=2, nrow=flen*(flen-1)/2)
    kk = 1        
    for (i in 1:(flen-1)) {
        for (j in (i+1):flen) {
            ans[kk,] = c(f[i], f[j])
            kk = kk+1
        }
    }
    ans = split(ans, seq(1, nrow(ans)))
    names(ans) = sapply(ans, paste, collapse=".")    
    return(ans)    
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
    
    if (class(f)=="character") {
        ans = matrix("", nrow=n2, ncol=d)
    } else {
        ans = matrix(0, nrow=n2, ncol=d)
    }
    
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
    ans2 = split(ans, 1:nrow(ans))
    ##names(ans2) = apply(ans, 1, paste, collapse=".")
    
    return(ans2)    
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
MPsuggestTopFeatures = function(x, ncomp=function(x) { 1+sqrt(x) },
    iqrfactor=function(x) { log(x) }) {
    
    ## transform the data using PCA
    if (class(x)!="prcomp") {
        x = prcomp(x)
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
    
    return(ans)    
}




##' Extract the k items that are nearest to specified seeds
##'
##' @param nowdist a distance object or matrix
##' @param seeds names of seed items
##' @param k number of nearest neighbors to select
##' 
##' @export
MPgetKNearest = function(nowdist, seeds, k) {

    ## always use nowdist as a matrix
    if (class(nowdist)=="dist") {
        nowdist = as.matrix(nowdist)
    }
    
    ## make sure that k is not too big
    k = min(k, ncol(nowdist)-1)
    
    ## look-up nearest neighbors
    ans = list()
    for (nowseed in seeds) {
        temp = sort(nowdist[nowseed,])
        temp = names(temp)[1:(k+1)]
        ans[[nowseed]] = temp[temp!=nowseed][1:k]
    }
    
    return(ans)
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




## helper function cuts a matrix into a small with at least min.n columns
## and with all columns that contain covthreshold of the total variance
##
## dd - input matrix
## subset - if less than one, proportion of variance explained by pca
##          if greater than one, number of PCA components
## min.n - minimum number of columns to output
##
## outputs a new matrix with the same number of rows as dd, with fewer columns
##
getPCAsubset = function(dd, subset=0.6, min.n=2) {
    if (ncol(dd)<min.n) {
        warning("input data has ", ncol(dd), "columns and min.n is ", min.n,"")
        return(NULL)
    }
    ddpca = prcomp(dd)$x
    if (subset<1) {
        ddcov = sum(diag(cov(dd)))
        ddpcacov = diag(cov(ddpca))
        ## find which columns in ddpcacov explain the featuers
        nowselect = cumsum(ddpcacov)<=subset*ddcov            
    } else {
        nowselect = rep(FALSE, ncol(ddpca))
        nowselect[1:min(subset, ncol(ddpca))] = TRUE
    }
    nowselect[1:min.n] = TRUE
    ## return the transformed data
    return(ddpca[, nowselect, drop=FALSE])
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
    return(list(xlim=xlim, ylim=ylim))
}




##' Create a map for a distance object/matrix
##'
##' This function uses cmdscale, but additionally adds a certain amount of
##' whitening determined by whiten. The whitening is pseudo-random and deterministic;
##' it always gives the same result. 
##'
##' @param d distance object
##' @param whiten numeric, use low number like 0.01 (default) to slightly jiggle points
##' on the map. This is a percentage of the x-y ranges.
##'
##' @export
MPgetMap = function(d, whiten=0.01) {
    
    ## get a cmdscale map
    ans = cmdscale(d)
    
    ## custom pseudo-random number generator with fixed seeds and settings
    ## this choice of a, b, m produces at least 1M unique random-like numbers
    myrn = function(n, x0=1, a=5862, b=3, m=11930007) {
        x = x0
        nextrn = function() {
            (a*x + b) %% m
        }    
        ans = rep(0, n);
        for (i in seq(n)) {
            x = nextrn()
            ans[i] = x
        }    
        ans/m
    }
    
    ## adjust results from cmdscale slightly
    if (whiten>0 & nrow(ans)>0) {
        ## get pseudo-random numbers and adjust cmdscale coordinates
        rn = qnorm(matrix(myrn(2*nrow(ans)), ncol=2))        
        xlim = range(ans[,1])
        ylim = range(ans[,2])
        zz = min(xlim[2]-xlim[1], ylim[2]-ylim[1])
        ans[,1:2] = ans[,1:2] + (rn*zz*whiten)
    }  
    
    return(ans)
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
    
    return (xp)
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
    
    ## if distance matrix is not specified, abort
    ##if (is.null(dist)) {
    ##    return(cutree(tree, k=k))
    ##}
    
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






##' Get information about top features contributing to PC components 
##'
##' @param x prcomp object
##' @param pc integer, the principal component to "explain", can be a vector of integers
##' @param n integer, number of features to output
##' @param plot logical, when true it automatically plots the data
##' @param mar margin values
##' 
##' @export
MPgetTopPC = function(x, pc=1, n=NULL, plot=FALSE, mar=c(4, 24, 4, 4)) {
    if (class(x)!="prcomp") {
        stop("x must be of class prcomp")
    }
    
    if (is.null(n)) {
        n = nrow(x$rotation)
    }
    
    ans = list()
    for (nowpc in pc) {
        temp = x$rotation[, nowpc]
        temp.abs = abs(temp)
        temp = temp[order(temp.abs, decreasing=T)][1:n]
        ans[[p0("PC", nowpc)]] = temp
    }
    
    if (plot) {
        Rcssplot::Rcsspar(mfrow=c(1, length(pc)), mar=mar)
        for (nowpc in pc) {
            Rcssplot::Rcssbarplot(ans[[p0("PC", nowpc)]])
        }
        Rcssplot::Rcsspar(mfrow=c(1,1))
    }
    
    invisible(ans)
}





##' get a subset of genesets that are small size and non-redundant
##'
##' This is an inelegant implementation, maybe fix and simplify later.
##' 
##' @param sets list with gene sets
##' @param target.genes character vector, names of ok gene names
##' @param target.size integer vector of length 2, indicate minimum and maximum size of geneset
##' @param target.JI numeric, maximum Jaccard overlap between the new genesets
##' 
##' @export
MPgetNonredundantGenesets = function(sets, target.genes=NULL, target.size=c(2, 20), target.JI=0.2) {
    
    ## get sets that overla with target gene set
    ## collect information on how large the sets are
    smallsetinfo = matrix(0, ncol=2, nrow=length(sets))
    rownames(smallsetinfo) = names(sets)
    colnames(smallsetinfo) = c("original", "ontarget")
    smallsets = list()
    for (nowsetname in names(sets)) {
        ## get overlap of this set with annotations
        nowset = sets[[nowsetname]]
        nowset2 = sort(nowset[nowset %in% target.genes])
        nowlen = length(nowset2)
        if (nowlen>=target.size[1] & nowlen<=target.size[2]) {
            smallsets[[nowsetname]] = nowset2          
            smallsetinfo[nowsetname, ] = c(length(nowset), length(nowset2))
        }
    }
    smallsetinfo = smallsetinfo[names(smallsets),]
    
    ## helper function for Jaccard index
    JI = function(x, y) {
        t1 = sum(x %in% y)
        t2 = length(unique(c(x, y)))
        return(t1/t2)
    }

    ## order smallsets from small to large
    smallsets = smallsets[order(sapply(smallsets, length), decreasing=T)]

    setsim = mtrx(0, col.names=names(smallsets), row.names=names(smallsets))
    oksets = rep(TRUE, length(smallsets))
    names(oksets) = names(smallsets)
    for (i in 1:length(smallsets)) {
        cat(".")
        if (i%%50==0) {
            cat(" ", i, "/", length(smallsets), "\n")
        }
        nowset1 = smallsets[[i]]
        nowsetname1 = names(smallsets)[i]
        for (j in 1:length(smallsets)) {
            nowsetname2 = names(smallsets)[j]
            if (j>i & oksets[nowsetname1] & oksets[nowsetname2]) {
                nowset2 = smallsets[[j]]
                nowJI = JI(nowset1, nowset2)                              
                if (nowJI>target.JI) {
                    nowval1 = smallsetinfo[nowsetname1, 1:2]
                    nowval2 = smallsetinfo[nowsetname2, 1:2]                    
                    ## select the bigger set
                    if (nowval1[2] > nowval2[2]) {                           
                        oksets[nowsetname2] = FALSE
                    } else if (nowval2[2] > nowval1[2]) {
                        oksets[nowsetname1] = FALSE
                    } else {
                        ## the sets are of equal size, then select on the original set
                        if (nowval1[1]>nowval2[1]) {
                            oksets[nowsetname2] = FALSE
                        } else {
                            oksets[nowsetname1] = FALSE
                        }
                    }                    
                }
            }            
        }
    }    
    smallsets = smallsets[oksets]
    cat("\n")
    
    ## if there are target genes that are leftover, they are put into separate genesets
    missing.genes = target.genes[!target.genes %in% unlist(smallsets)]
    smallsets[["Other"]] = missing.genes
    
    return(smallsets)
}




##' shows the genesets that contain genes
##'
##' @param sets list of gene sets
##' @param genes vector of genes, all genes must be present in a set to provide a hit
##'
##' @export
MPingenesets = function(sets, genes) {
    sets[sapply(sets, function(x) { sum(genes %in% x)==length(genes)})]         
}




##' Creates a set of configurations using prcomp
##' 
##' @param MP MultiPattern
##' @param dd matrix, dataset in MP to process
##' @param ddname character, name of dataset
##'
MPaddPCAset = function(MP, dd, ddname) {
    
    datapca = getPCAsubset(dd, subset=PCAexplain)

    if (!is.null(datapca)) {
        pcaname = paste0(ddname, ".pca")
        ## add dataset to the MPconfig
        MPaddData(MP, setNames(list(datapca), pcaname)) 
        ## create a series of configurations for pca
        pcasets = list()
        imax = min(ncol(dd)-1, ncol(datapca))
        for (i in 1:imax) {
            if (i==1) {
                pcasets[[colnames(datapca)[1]]] = colnames(datapca)[1]                        
            } else if (i>2) {
                pcasets[[paste(colnames(datapca)[c(1,i)], collapse="..")]] = colnames(datapca)[1:i]
            }
        }
        ## add configuration with increasing PCA details
        MPaddConfig(MP, paste0(config.prefix, ddname, ":", names(pcasets)),
                    data.name=pcaname, preprocess=pcasets)
        ## add configurations with pairwise PCA columns
        twospaces = MPgetAll2Subspaces(head(colnames(datapca)))
        if (length(twospaces)>0) {
            MPaddConfig(MP, paste0(config.prefix, ddname, ":", names(twospaces)),
                        data.name=pcaname, preprocess=twospaces)
        }
    }

    return(MP)
}



