## Functions that hold "management" objects for MultiMetrics clustering
## (objects that help to keep track of the clustering definitions and distance matrics
##





##' Create a basic MP object
##'
##' @param items character vector, specifies names for observations in the multipatter analysis
##' @param data named list with data objects. Specifying this parameter
##' is equivalent to calling MPnew() and then adding data manually with MPaddData()
##' 
##' @export 
MPnew = function(items, data=NULL) {
    
    ## create a blank MultiPattern object
    ans = list()
    ans$items = items;
    ans$data = list()
    ans$configs = list()
    ans$settings = list(
        ## number of random configurations to create by default
        num.random=100,
        ## number of PCA components to keep for automated PCA
        num.PCs=4,
        ## used in rpca analysis (smaller than default for faster exec)
        rpca.term.delta=1e-3, 
        ## max number of clusters to use in easyConfig reg and alt configurations
        clust.k=3,
        ## number of random subspaces
        subspace.num.random=100,
        ## number or proportion of features for subspace analysis
        subspace.d.random=0.5,
        ## exponent for similarity transformation
        alpha=1        
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
        "subspace.num.random", "subspace.d.random")
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


## This block below was part of MPeasyConfig (the settings are now part of the MPchangeSettings)
##
## @param Nrandom integer, number of random configurations to add
## @param PCAexplain numeric, used for pca configurations. When <1, interpreted as the proportion
## of variance explained by the PCA componets. When >1, interpreted as the number of PCA
## components to use.
## @param term.delta numeric, used during rpca decomposition. See package rpca. This is set to a
## larger number than default in rpca to speed up execution.
## @param clust.k integer, number of cluster groups to consider with type="clust" and "pam".
## When set to clust.k=4, configuration will be used with k=2,3,4.


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
##' @param subspace.prefix character, a middle-fix used when naming subspaceR configurations
##' @param type character or list, codes for what types of configurations to add. The default list
##' shows all the supported configuration types. Detailed descriptions will appear elsewhere.
##' @param preprocess object specifying preprocessing, e.g. vector of features for subspaces
##' 
##' @export
MPeasyConfig = function(MP, data=NULL, config.prefix="", subspace.prefix="", 
    type=c("pca", "rpca", "euclidean", "spearman", "canberra", "manhattan", "hclust", "pam",
        "random", "subspace1", "subspace2", "subspaceR"),
    preprocess=NULL) {
    
    ## get values from MP settings
    Nrandom = MP$settings$num.random
    PCAexplain = MP$settings$num.PCs
    term.delta = MP$settings$rpca.term.delta
    clust.k = MP$settings$clust.k

    ## define types of analysis that are supported by package functions
    oktypes = c("pca", "rpca", "euclidean", "spearman", "canberra", "manhattan", "hclust", "pam",
        "random", "subspace1", "subspace2", "subspacer")
    
    ## capture MP expression for assignment at the end
    captureMP = deparse(substitute(MP))
    
    ## check if all types requested by user are supported 
    nowtypes = tolower(unlist(type))
    type.bad = nowtypes[!(nowtypes %in% oktypes)]
    if (length(type.bad)>0) {
        stop("Unrecognized configuration types: ", type.bad, "\n")
    }
    rm(type.bad)
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
                    stop("Unrecognized data type: ", data.type, "\n")
                }
            }
        }
    }

    if (subspace.prefix != "") {
        subspace.prefix = p0(":", subspace.prefix)
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
    
    clust.k.range = seq(2, abs(clust.k))
    clust.k.range = clust.k.range[clust.k.range>1]
    
    ## ###################################################################################
    ## helper function adds a PCA dataset and a set of configurations to MP
    addPCAset = function(dd, ddname) {
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

    ## ###################################################################################
    ## Here, all relevant info is in a list typelist
    for (nowd in names(typelist)) {
        nowtypes = typelist[[nowd]]
        for (oneliner in c("euclidean", "spearman", "canberra", "manhattan")) {
            if (oneliner %in% nowtypes) {
                MPaddConfig(MP, paste0(config.prefix, nowd,":", oneliner), data.name=nowd,
                            dist.fun=MPdistFactory(method=oneliner))                
            }
        }                
        if ("pca" %in% nowtypes) {
            MP = addPCAset(MP$data[[nowd]], nowd)
        }
        if ("rpca" %in% nowtypes) {
            ##require("rpca")
            datarpca = rpca::rpca(MP$data[[nowd]], term.delta=term.delta);
            rownames(datarpca$L) = rownames(datarpca$S) = rownames(MP$data[[nowd]])
            colnames(datarpca$L) = colnames(datarpca$S) = colnames(MP$data[[nowd]])            
            Lname = paste0(nowd, ".rpcaL")
            Sname = paste0(nowd, ".rpcaS")
            MPaddData(MP, setNames(list(datarpca$L, datarpca$S), c(Lname, Sname)))           
            MP = addPCAset(MP$data[[Lname]], Lname)
            MP = addPCAset(MP$data[[Sname]], Sname)
        }
        if ("hclust" %in% nowtypes) {
            clustconf = paste0(config.prefix, nowd, ":clust.")
            ## add cluster-based distance for each k
            for (nowk in clust.k.range) {
                if (nrow(MP$data[[nowd]])>(2*nowk)+2) {
                    for (ncm in c("complete", "single", "average")) {
                        nowdreg = MPdistFactory(method="hclust",
                            clust.k=nowk, clust.method=ncm, clust.alt=FALSE)
                        nowdalt = MPdistFactory(method="hclust",
                            clust.k=nowk, clust.method=ncm, clust.alt=TRUE)
                        ncm.init = toupper(substring(ncm, 1, 1))                        
                        MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "reg"),
                                    data.name=nowd, dist.fun=nowdreg)
                        MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "alt"),
                                    data.name=nowd, dist.fun=nowdalt)
                        rm(nowdreg, nowdalt, ncm.init)
                    }
                }
            }
            rm(clustconf)
        }
        if ("pam" %in% nowtypes) {
            clustconf = paste0(config.prefix, nowd, ":clust.")
            ## add cluster-based distance for each k
            for (nowk in clust.k.range) {
                if (nrow(MP$data[[nowd]])>(2*nowk)+2) {
                    ncm.init = "P"                    
                    nowdreg = MPdistFactory(method="pam",
                        clust.k=nowk, clust.method="pam", clust.alt=FALSE)
                    nowdalt = MPdistFactory(method="pam",
                        clust.k=nowk, clust.method="pam", clust.alt=TRUE)
                    MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "reg"),
                                data.name=nowd, dist.fun=nowdreg)
                    MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "alt"),
                                data.name=nowd, dist.fun=nowdalt)
                    rm(nowdreg, nowdalt)                    
                }
            }
            rm(clustconf)            
        }
        if ("subspace1" %in% nowtypes) {
            if (is.null(preprocess)) {
                temp = as.list(colnames(MP$data[[nowd]]))
                names(temp) = colnames(MP$data[[nowd]])
            } else {
                temp = as.list(preprocess)
                names(temp) = preprocess
            }
            MPaddConfig(MP, paste0(config.prefix, nowd, subspace.prefix, ":subspace1"), data.name=nowd,
                        preprocess=temp)
            rm(temp)
        }
        if ("subspace2" %in% nowtypes) {
            if (is.null(preprocess)) {
                temp = MPgetAll2Subspaces(colnames(MP$data[[nowd]]))
            } else {
                temp = MPgetAll2Subspaces(preprocess)
            }
            MPaddConfig(MP, paste0(config.prefix, nowd, subspace.prefix, ":subspace2"), data.name=nowd,
                        preprocess=temp)
            rm(temp)
        }
        if ("subspaceR" %in% nowtypes) {
            Nsub = ceiling(MP$settings$subspace.num.random)
            Dsub = MP$settings$subspace.d.random
            if (is.null(preprocess)) {
                temp = colnames(MP$data[[nowd]])
            } else {
                temp = preprocess
            }
            if (Dsub<1) {
                Dsub = abs(Dsub*length(temp))
            }            
            Dsub = min(length(temp), ceiling(Dsub))
            temp.subspaces = MPgetRandomSubspaces(temp, Nsub, Dsub)
            MPaddConfig(MP, paste0(config.prefix, nowd, subspace.prefix, ":subspaceR"), data.name=nowd,
                        preprocess=temp.subspaces, dist.fun=dist.euclidean)
            rm(temp.subspaces, Nsub, Dsub)
        }
    }
    
    ## add random configurations if needed
    if (class(type)!="list" & length(type)==1) {
        if (type=="random" & Nrandom>0) {
            MPaddConfig(MP, paste0(config.prefix, "rnorm.", 1:Nrandom), names(MP$data)[1], 
                        preprocess=vector("list", Nrandom), dist.fun=dist.rnorm)
        }
    }
    
    
    ## ###################################################################################
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
##' @param verbose logical, when TRUE function prints a warning message when the input is large
##'
##' @export
MPgetDistances = function(MP, verbose=TRUE) {
    
    if (class(MP) != "MultiPattern") {
        stop("object not of class MultiPattern\n")
    }

    if (verbose & object.size(MP)>1e6) {
        cat("This may take a little time. Please wait... ")
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
    ans = lapply(MP$configs, adfun)
     
    ## housekeeping for the list
    names(ans) = names(MP$configs)
    class(ans) = "MultiPatternSimilarities"
    
    if (verbose & object.size(MP)>1e6) {
        cat("done\n")
    }
    
    ans
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
        temp = (a-b)^beta        
        sum(temp)^(1/beta)
    }
    simpledist = function(a, b) {
        temp = a-b
        sqrt(sum(temp*temp))
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
##' @param d dist object 
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
    
    if (class(d)!="dist") {
        stop("input must be a distance obejct")
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



