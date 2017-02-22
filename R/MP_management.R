## Functions that hold "management" objects for MultiPattern clustering
## (objects that help to keep track of the clustering definitions and distance matrics
##


## General package maintenance - import functions from packages
#' @import stats
#' @import utils
#' @import cluster
#' @import graphics
#' @import Rcssplot
NULL



##' Create a basic MP object
##'
##' @param items character vector, specifies names for observations in the multipattern analysis
##' @param data named list with data objects. Specifying this parameter
##' is equivalent to calling MPnew() and then adding data manually with MPaddData()
##' 
##' @export 
MPnew = function(items, data=NULL) {
    
    ## create a blank MultiPattern object
    ans = list(items=items, data=list(), configs=list())
    ans$settings = list(
        ## number of random configurations to create by default
        num.random=100,
        ## number of PCA components to keep for automated PCA
        num.PCs=4,
        ## used in rpca analysis (smaller than default for faster exec)
        rpca.term.delta=1e-3, 
        ## max number of clusters to use in easyConfig reg and alt configurations
        clust.k=3,
        ##  used in nmf analysis (background level to avoid all-0 rows)
        nmf.bg = 1e-5,
        ## used in nmf analysis (maximal rank of nmf decomposition, zero is auto)
        nmf.rank = 0,
        ## number of random subspaces
        subspace.num.random=100,
        ## number or proportion of features for subspace analysis
        subspace.d.random=0.5,
        ## exponent for similarity transformation and meta-similarity Lbeta distance
        alpha=1,
        beta=2,
        ## subsampling for meta-distance calculation, number of samples and repetitions
        subsample.N=150,
        subsample.R=30
        )
    class(ans) = "MultiPattern"
    
    ## perhaps add datasets into this object
    if (!is.null(data)) {
        MPaddData(ans, data)
    }
    
    ans
}




##' Add a dataset into a MultiPattern configuration object
##' 
##' @param MP an existing MultiPattern object
##' @param data named list with data objects. 
##'
##' @export
MPaddData = function(MP, data) {
    
    if (class(MP) != "MultiPattern") {
        stop("Argument MP must be of class MultiPattern\n")
    }
    if (class(data) != "list") {
        stop("Argument data must be a list\n")
    }    

    ## capture MP expression for assignment at the end
    captureMP = deparse(substitute(MP))
    
    ## check naming of data objects
    datanames = names(data)
    if (is.null(datanames)) {
        stop("Names of data objects cannot be NULL\n")
    }
    badnames = datanames[datanames %in% names(MP$data)]
    if (length(badnames)>0) {
        stop("Duplicate data object names: ", paste(badnames, collapse=", "), "\n")
    }

    ## update the list of data objects in this MultiPattern configuration
    MP$data = c(MP$data, data)

    ## assign and return an updated MP object
    assign(captureMP, MP, parent.frame())
    invisible(MP)          
}





##' Add one or more cluster configurations into a MultiPattern object
##'
##' The function acceps data.fun or dist.fun as a list. In these cases, the
##' configuration names are given as config.name:[id]. 
##' 
##' @param MP existing MP object
##' @param config.name character, base name for new subspace/distance
##' @param data.name character, name for data component in MP object. Defaults to first dataset
##' declared in MP. 
##' @param preprocess function or feature set, data preprocessing applied before distance function.
##' Can be a single value or a list. Defaults to NULL, which indicates no preprocessing is required.
##' When a function, the function is applied on the dataset prior to computing dissimilarities.
##' When a vector (character or integer), it is interpreted a feature set of the dataset;
##' dissimilarities are then computed on this subspace of the data.
##' @param dist.fun function, distance function. Can be a single value or a list.
##' Defaults to euclidean distance. 
##'
##' @export
MPaddConfig = function(MP, config.name, data.name=names(MP$data)[1],
    preprocess=NULL, dist.fun=dist.euclidean) {
    
    ## Hard checks for object class types
    if (class(MP) != "MultiPattern") {
        stop("Argument MP must be of class MultiPattern\n");
    }
    if (class(data.name) != "character") {
        stop("Argument data.name must be character\n")
    }        
    ## Hard checks for data compatibility
    if (!data.name %in% names(MP$data)) {
        stop("data.name -", data.name, "- not defined in MultiPattern object\n")
    }
    ## Hard check, accept only one of preprocess or dist.fun as a list
    if (class(preprocess)=="list" & class(dist.fun)=="list") {
        stop("only one of preprocess or dist.fun can be a list\n")
    }
    ## Hard check for config.name
    if (sum(is.null(config.name))+sum(is.na(config.name))>0) {
        stop("config.name must be valid string code\n")
    }
    
    ## capture MP expression for assignment at the end
    captureMP = deparse(substitute(MP))

    ## get names of new configurations.
    ## If inputs specify just one configuration, its name is just config.name.
    ## If inputs specify lists, config names are formated as config.name:[identifier]
    ## Identifiers are from names of the data.fun or dist.fun lists
    now.list = NULL
    if (class(preprocess)=="list") {
        now.list = preprocess
    } else if (class(dist.fun)=="list") {
        now.list = dist.fun
    }
    if (!is.null(now.list)) {
        numnew = length(now.list)
        if (length(now.list)==length(config.name)) {
            newnames = config.name
        } else  {
            if (is.null(names(now.list))) {
                newnames = paste0(config.name, ".", seq_along(now.list))
            } else {
                newnames = paste0(config.name, ".", names(now.list))
            }
        }
    } else {
        newnames = config.name
        dist.fun = list(a=dist.fun)
    }
    rm(now.list)
    
    ## check if any of the new names overlap with existing names
    badnames = newnames[newnames %in% names(MP$configs)]
    if (length(badnames)>0) {        
        stop("Duplicate clustering configuration names: ", paste(badnames, collapse=", "), "\n")
    }
    
    ## create a list with new configurations
    newconfigs = vector("list", length(newnames))
    names(newconfigs) = newnames;
    ## helper function to create an object with details for one analysis
    makeOneConf = function(nam, dat, prep, distfun) {
        ## create a list with all the input data
        aa = list()
        aa$name = nam
        aa$data = dat
        aa$prep = MPmakePrepFunction(prep)
        aa$dist.fun = distfun
        if (class(prep)!="function" & class(prep)!="NULL") {
            aa$info = prep
        }
        return(aa)
    }
    
    ## create configuration objects and add to newconfigs list
    if (class(preprocess)=="list") {
        for (i in seq_along(preprocess)) {
            nowname = newnames[i]
            newconfigs[[nowname]] = makeOneConf(nowname, data.name, preprocess[[i]], dist.fun)
        }        
    } else if (class(dist.fun)=="list") {
        for (i in seq_along(dist.fun)) {            
            nowname = newnames[i]
            newconfigs[[nowname]] = makeOneConf(nowname, data.name, preprocess, dist.fun[[i]])
        }        
    }
    
    ## update MP with new configurations
    MP$configs = c(MP$configs, newconfigs)
    
    ## assign and return an updated MP object
    assign(captureMP, MP, parent.frame())
    invisible(MP)            
}





##' Remove data or configuration components from a MultiPattern object
##' 
##' @param MP a MultiPattern configuration objects
##' @param data character or vector of data identifiers 
##' @param config character or vector of config identifiers 
##'
##' @export
MPremove = function(MP, data=NULL, config=NULL) {

    ## Hard checks for object class
    if (class(MP) != "MultiPattern") {
        stop("Argument MP must be of class MultiPattern\n");
    }

    ## capture MP expression for assignment at the end
    captureMP = deparse(substitute(MP))

    ## remove configurations
    if (!is.null(config)) {
        if (class(config) != "character") {
            stop("Argument config must be character\n")
        }
        MP$configs[config] = NULL
    }
    
    ## remove data components
    if (!is.null(data)) {
        ## remove datasets from MP$data
        for (nowd in data) {
            if (nowd %in% names(MP$data)) {
                MP$data[[nowd]] = NULL
            } else {
                warning("Data object", nowd, " is not in MultiPattern configuration\n");
            }
        }
        
        ## remove analysis configurations from MP$configs
        okconfigs = !sapply(MP$configs, function(x) x$data %in% data)
        MP$configs = MP$configs[okconfigs]
    }
        
    ## assign and return an updated MP obejct
    assign(captureMP, MP, parent.frame())
    invisible(MP)
}





##' Change settings encoded in a MultiPattern object
##' 
##' @param MP a MultiPattern configuration object
##' @param settings list with new settings
##'
##' Accepted settings names are 'num.PCs', 'num.random', 'alpha',
##' 'rpca.term.delta', 'clust.k', 'subspace.num.random', 'subspace.d.random'
##' 
##' @export
MPchangeSettings = function(MP, settings = list()) {
    
    ## Hard checks for object class
    if (class(MP) != "MultiPattern") {
        stop("Argument MP must be of class MultiPattern\n");
    }
    if (class(settings) != "list") {
        stop("Argument settings my be of class list\n")
    }
    
    ## capture MP expression for assignment at the end
    captureMP = deparse(substitute(MP))
    
    ## check that components in settings are allowed
    sn = names(settings)
    goodsettings = c("num.PCs", "num.random", "alpha", "rpca.term.delta", "clust.k",
        "subspace.num.random", "subspace.d.random", "subsample.N", "subsample.R")
    badsettings = sn[!sn %in% goodsettings]
    if (length(badsettings)>0) {
        warning("MPchangeSettings: Unrecognized items -", paste(badsettings, collapse=", "), "- will be skipped\n")
    }
    for (nows in goodsettings) {
        if (nows %in% names(settings)) {
            MP$settings[[nows]] = settings[[nows]]
        }
    }
    
    ## assign and return an updated MP obejct
    assign(captureMP, MP, parent.frame())
    invisible(MP)
}
    




## ###################################################################################
## Functions for adding configurations
## ###################################################################################


##' Add predefined configuration types to a MultiPattern configuration object
##'
##' This function modifies the object MP in the first argument. The primary modifications
##' are new items in the MP$config list. When pca or rpca are specified in the type argument,
##' the function also precomputes these transformations and adds these into the MP$data list. 
##' All configurations relying on pca or rpca transformed data refer to these pre-computed objects.
##' 
##' @param MP MultiPattern configuration object
##' @param data character, names of datasets in MP to use to create suggestions. If NULL, function
##' applies the configurations to all datasets (but see type). If not NULL, function only considers
##' the specified data objects.
##' @param config.prefix character, prefix used in all configuration names
##' @param preprocess.prefix character, a middle-fix used when naming subspaceR configurations
##' @param type character or list, codes for what types of configuration plugins to use. To see all
##' available plugins, use MPlistPlugins(). 
##' shows all the supported configuration types. Detailed descriptions will appear elsewhere.
##' @param random logical, set TRUE to include random configurations. 
##' @param preprocess object specifying preprocessing, e.g. vector of features for subspaces
##' 
##' @export
MPeasyConfig = function(MP, data=NULL, config.prefix="",
    preprocess.prefix="", 
    type=c("pca", "euclidean", "spearman", "canberra",
        "manhattan", "hclust", "pam"),
    random=TRUE,
    preprocess=NULL) {
    
    ## capture MP expression for assignment at the end
    captureMP = deparse(substitute(MP))
   
    if (!is.null(data)) {
        data.missing = data[!(data %in% names(MP$data))]
        if (length(data.missing)) {
            stop("Missing data: ", data.missing, "\n")
        }
        rm(data.missing)
    }
    if (class(type)=="list") {
        if (!is.null(data)) {
            stop("type in list form requires data=NULL\n")
        }
        if (length(type)>0) {
            if (is.null(names(type))) {
                stop("type in list form requires names\n")
            } else {
                data.type = names(type)
                data.type = data.type[!(data.type %in% names(MP$data))]
                if (length(data.type)>0) {
                    stop("Unrecognized datasets: ", paste(data.type, collapse=", "), "\n")
                }
            }
        }
    }
    
    if (preprocess.prefix != "") {
        preprocess.prefix = paste0(":", preprocess.prefix)
    }
    
    ## standardize the input - into data=NULL and type=list(data=type)
    typelist = list()
    if (is.null(data)) {
        if (class(type)=="list") {
            typelist = type
        } else {
            typelist = setNames(vector("list", length(MP$data)), names(MP$data))
            typelist = lapply(typelist, function(x) {type})
        }
    } else {
        typelist = setNames(vector("list", length(data)), data)
        typelist = lapply(typelist, function(x) {type})
    }
    
    ## Add configurations by applying plugins    
    for (nowd in names(typelist)) {
        nowtypes = tolower(typelist[[nowd]])        
        for (nowtype in nowtypes) {
            plugin.fun = match.fun(paste0(nowtype, ".MultiPatternPlugin"))
            MP = plugin.fun(MP, nowd, config.prefix,
                preprocess.prefix, preprocess)
        }             
    }    
    ## perhaps add a set of random configurations    
    if (random) {
        plugin.fun = match.fun("random.MultiPatternPlugin")
        MP = plugin.fun(MP, NULL, config.prefix,
            preprocess.prefix, preprocess) 
    }
        
    ## assign and return an updated MP obejct
    assign(captureMP, MP, parent.frame())
    invisible(MP)                
}





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
    
    ## check that the sub sampling is compatible with MP, number of iterations is >= 1
    subsample.N = min(length(MP$items)-1, ceiling(subsample.N))
    subsample.R = max(1, ceiling(subsample.R))

    if (subsample.N<2) {
        stop("Subsampling gives too-small a dataset");
    }

    if (verbose) {
        cat("This may take some time. Please wait... ")
    }
    
    ## define a result object
    result = NULL;
    
    ## compute meta-distances in a loop
    for (i in 1:subsample.R) {
        ## get a subset of the items and create a new MP object
        tempMP = MP
        tempMP$items = sample(MP$items, subsample.N, replace=F)
        ## get smaller versions of the data matrices
        for (nowdata in names(MP$data)) {
            tempMP$data[[nowdata]] = MP$data[[nowdata]][tempMP$items,]
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
        cat("done\n")
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
##' @param k numeric of integer. If more than 1, number of points to select
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



