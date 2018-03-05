## Package: MultiPattern
##
## Functions for managing/bookkeeping for MultiPattern clustering
##


## General package maintenance - import functions from packages
#' @import stats
#' @import tsne
#' @import utils
#' @import cluster
#' @import graphics
#' @import Rcssplot
NULL




##' Summary of default settings for new MultiPattern objects
##'
##' See details for a description of each component
##' 
##' num.random: number of random configurations
##' 
##' num.PCs: number of PCA components
##' 
##' rpca.term.delta: used in rpca analysis (smaller than default for speed)
##'
##' clust.k: max number of clusters to use in easyConfig reg and alt
##' configurations
##'
##' nmf.bg: background level in nmf analysis (avoids all-zero rows)
##'
##' nmf.rank: maximal rank of NMF analysis (zero triggers an
##' automatic decision)
##'
##' subspace.num.random: number of random subspaces
##'
##' subspace.d.random: number/proportion of features for subspace analysis
##' 
##' alpha: exponent for similarity transformation for meta-similarities
##'
##' beta: determines Lbeta distance for meta-similarities
##'
##' subsample.N: number of observations to use in subsample/bootstrap
##'
##' subsample.R: number of subsampling/bootstrap repetitions
##' 
##' @export
MPdefaultSettings = list(
  num.random=60,
  num.PCs=4,
  rpca.term.delta=1e-3,
  dbscan.intervals=c(0.1,0.2,0.3,0.4),
  clust.k=3,
  nmf.bg = 1e-5,
  nmf.rank = 0,
  subspace.num.random=100,
  subspace.d.random=0.5,
  alpha=0.5,
  beta=2,
  subsample.N=150,
  subsample.R=20    
)




##' Create a basic MP object
##'
##' @param items character vector, specifies names for observations in the multipattern analysis
##' @param data named list with data objects. Specifying this parameter
##' is equivalent to calling MPnew() and then adding data manually with MPaddData()
##' 
##' @export 
MPnew = function(items, data=NULL) {
    
  ## create a blank MultiPattern object
  MP = list(items=items, data=list(), configs=list())
  MP$settings = MPdefaultSettings
  class(MP) = "MultiPattern"
  class(MP$settings) = "MultiPatternSettings"
  
  ## perhaps add datasets into this object
  if (!is.null(data)) {
    MP = MPaddData(MP, data)
  }
  
  MP
}




##' Add a dataset into a MultiPattern configuration object
##' 
##' @param MP an existing MultiPattern object
##' @param data named list with data objects. 
##'
##' @export
MPaddData = function(MP, data) {

  checkArgClass(MP, "MultiPattern")
  checkArgClass(data, "list")
  
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
  
  MP
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
  
  checkArgClass(MP, "MultiPattern")
  checkArgClass(config.name, "character")
  checkArgClass(data.name, "character")
  
  ## Hard checks for data compatibility
  if (!data.name %in% names(MP$data)) {
    stop("data.name -", data.name, "- not defined in MultiPattern object\n")
  }
  ## Hard check, accept only one of preprocess or dist.fun as a list
  if (class(preprocess)=="list" & class(dist.fun)=="list") {
    stop("only one of preprocess or dist.fun can be a list\n")
  }
  
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
    aa
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
  
  MP
}




##' Remove data or configuration components from a MultiPattern object
##' 
##' @param MP a MultiPattern configuration objects
##' @param data character or vector of data identifiers 
##' @param config character or vector of config identifiers 
##'
##' @export
MPremove = function(MP, data=NULL, config=NULL) {
  
  checkArgClass(MP, "MultiPattern")
  
  ## remove configurations
  if (!is.null(config)) {
    checkArgClass(config, "character")
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
  
  MP
}




##' Change settings encoded in a MultiPattern object
##' 
##' @param MP a MultiPattern configuration object
##' @param settings list with new settings
##' @param warn logical, set TRUE to get warnings if a setting value
##' is not part of the core MultiPattern set (MPdefaultSettings)
##' 
##' Accepted settings names are those defined in MPdefaultSettings.
##' 
##' @export
MPchangeSettings = function(MP, settings = list(), warn=TRUE) {
    
  checkArgClass(MP, "MultiPattern")
  checkArgClass(settings, "list")
  
  ## check that components in settings are allowed
  sn = names(settings)
  noncore = sn[!sn %in% names(MPdefaultSettings)]
  if (length(noncore)>0 & warn) {
    warning("\nnon-core setting(s): ", paste(noncore, collapse=", "), "\n")
  }
  for (nows in names(settings)) {
    MP$settings[[nows]] = settings[[nows]]
  }
  
  MP
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
##' @param preprocess object specifying preprocessing, e.g. vector of features for subspaces
##' 
##' @export
MPeasyConfig = function(MP, data=NULL, config.prefix="",
                        preprocess.prefix="", 
                        type=c("pca", "euclidean", "spearman", "canberra",
                               "manhattan", "hclust", "pam", "dbscan"),
                        preprocess=NULL) {
  
  ## capture MP expression for assignment at the end
  ##captureMP = deparse(substitute(MP))
  
  if (!is.null(data)) {
    data.missing = data[!(data %in% names(MP$data))]
    if (length(data.missing)) {
      stop("Missing data: ", data.missing, "\n")
    }
  }
  if (class(type)=="list") {
    if (!is.null(data)) {
      stop("when type is a list, data must null\n")
    }
    if (length(type)>0) {
      checkNotNull(names(type), "names in a type list")      
      data.type = names(type)
      data.type = data.type[!(data.type %in% names(MP$data))]
      if (length(data.type)>0) {
        stop("Unrecognized datasets: ", paste(data.type, collapse=", "), "\n")
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
      MP = plugin.fun(MP, nowd, config.prefix, preprocess.prefix, preprocess)
    }             
  }    
  
  MP
}

