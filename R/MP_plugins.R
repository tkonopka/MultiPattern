## Package: MultiPattern
##
## Functions that create MP configurations (plug-ins for MPeasyConfig)
## 




##' Add a set of configurations to a MultiPattern object
##'
##' This plugin adds one configuration that uses a traditional canberra distance.
##' 
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
canberra.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                       preprocess.prefix="", preprocess=NULL) {
  
  MPaddConfig(MP, paste0(config.prefix, data.name,":canberra"),
              data.name=data.name,
              dist.fun=MPdistFactory(method="canberra"))    
}




##' Add a set of configurations to a MultiPattern object
##'
##' This plugin adds several configurations that use dbscan-driven dissimilarities.
##' 
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
dbscan.MultiPatternPlugin = function(MP, data.name, config.prefix, preprocess.prefix="",
                                  preprocess=NULL) {

  ## get a range of eps scaling factors
  intervals = MP$settings$dbscan.intervals
  
  ## get some characteristic scales in the data
  data.q = apply(MP$data[[data.name]], 2, quantile, p=c(0.25, 0.5, 0.75))
  rownames(data.q) = c("Q25", "M", "Q75")
  data.qdist = as.matrix(stats::dist(data.q, method="euclidean"))
  data.scale = mean(data.qdist["M", c("Q25", "Q75")])
  eps.vals = data.scale*intervals
  
  dbconf = paste0(config.prefix, data.name, preprocess.prefix, ":dbscan.")
  
  ## add cluster-based distance for each eps multiple
  for (i in seq_along(eps.vals)) {
    noweps = eps.vals[i]
    nowdbscan = MPdistFactory(method="dbscan", eps=noweps, clust.method="dbscan")
    MP = MPaddConfig(MP, paste0(dbconf, i), data.name=data.name, dist.fun=nowdbscan)
  }
  
  MP
}





##' Add a set of configurations to a MultiPattern object
##'
##' This plugin adds one configuration that uses a traditional euclidean distance.
##' 
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
euclidean.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                        preprocess.prefix="", preprocess=NULL) {
  
  MPaddConfig(MP, paste0(config.prefix, data.name,":euclidean"),
              data.name=data.name,
              dist.fun=MPdistFactory(method="euclidean"))    
}




##' Add a set of configurations to a MultiPattern object
##'
##' This plugin adds several configurations that use hamming distance on
##' columns in the data that are factors or characters
##'
##' @param MP MultiPattern	
##' @param data.name character	
##' @param config.prefix character	
##' @param preprocess.prefix character	
##' @param preprocess character, (ignored)
##' 
##' @export
hamming.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                      preprocess.prefix="", preprocess=NULL) {
  
  ## fetch the data and make sure it is positive	
  nowdata = MP$data[[data.name]]
  
  ## identify factor or character columns
  if (is.null(preprocess)) {
    if ("data.frame" %in% class(nowdata)) {
      hamming.columns = sapply(nowdata, class) %in% c("character", "factor")
      hamming.columns = colnames(nowdata)[hamming.columns]
    } else {
      hamming.columns = rep(FALSE, ncol(nowdata))
      for (i in seq_len(ncol(nowdata))) {
        hamming.columns[i] = class(nowdata[, i]) %in% c("character", "factor")
      }
      hamming.columns = colnames(nowdata)[hamming.columns]
    }
  } else {
    hamming.columns = preprocess
  }
  hamming.length = length(hamming.columns)
  
  if (hamming.length==0) {
    return(MP)
  }

  ## always add a configuration based on all the hamming-eligible columns
  MP = MPaddConfig(MP, paste0(config.prefix, data.name, preprocess.prefix, ":hamming"),
                   data.name=data.name, preprocess=hamming.columns, dist.fun=dist.hamming)
  
  ## perhaps also add individual subspaces
  if (hamming.length>1) {
    hamming.singles = setNames(as.list(hamming.columns), hamming.columns)
    MP = MPaddConfig(MP, paste0(config.prefix, data.name, preprocess.prefix, ":hamming"),
                     data.name=data.name, preprocess=hamming.singles, dist.fun=dist.hamming)
  }

  MP
}




##' Add a set of configurations to a MultiPattern object
##'
##' This plugin adds several configurations that use clustering-driven dissimilarities
##' (using hclust). Some configurations (marked reg) create dissimilarities that reinforce
##' traditional clustering. Other configurations (marked alt) use an alternative clustering
##' distance. Settings for this plugin are extracted from MP$settings.
##' 
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
hclust.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                     preprocess.prefix="", preprocess=NULL) {

  ## get a range of k for the clustering
  clust.k = MP$settings$clust.k    
  clust.k.range = seq(2, abs(clust.k))
  clust.k.range = clust.k.range[clust.k.range>1]
  
  clustconf = paste0(config.prefix, data.name, ":clust.")
  ## add cluster-based distance for each k
  for (nowk in clust.k.range) {
    if (nrow(MP$data[[data.name]])>(2*nowk)+2) {
      for (ncm in c("complete", "single", "average")) {
        nowdreg = MPdistFactory(method="hclust",
                                clust.k=nowk, clust.method=ncm, clust.alt=FALSE)
        nowdalt = MPdistFactory(method="hclust",
                                clust.k=nowk, clust.method=ncm, clust.alt=TRUE)
        ncm.init = toupper(substring(ncm, 1, 1))                        
        MP = MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "reg"),
                    data.name=data.name, dist.fun=nowdreg)
        MP = MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "alt"),
                         data.name=data.name, dist.fun=nowdalt)
      }
    }
  }
  
  MP
}




##' Add a set of configurations to a MultiPattern object
##'
##' This plugin uses fastICA to identify interesting components, then
##' adds configurations that uses those components. The configurations are of two types:
##' individual independent components, e.g. only IC1 or only IC2
##' pairs of IC coponents, e.g. IC1, IC2 (IC1.IC2) or IC2, IC4 (IC2.IC4)
##'
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
ica.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                      preprocess.prefix="", preprocess=NULL) {

  ## fetch number of components the user wants to consider
  numICs = MP$settings$num.ICs
  
  ## get matrices with the all the raw data and with PCA transformed data
  dataica = getICAsubset(MP$data[[data.name]], numICs)
  if (is.null(dataica)) {
    return(MP)
  }
  
  ## add transformed data and configurations
  pca.ica.helper(MP, data.name, config.prefix,
                 preprocess.prefix=preprocess.prefix, preprocess=preprocess,
                 data=dataica, type="ica")
}




##' Add a set of configurations to a MultiPattern object
##'
##' Adds one configuration that uses a traditional manhattan distance.
##' 
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
manhattan.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                        preprocess.prefix="", preprocess=NULL) {
  
  MPaddConfig(MP, paste0(config.prefix, data.name,":manhattan"),
              data.name=data.name,
              dist.fun=MPdistFactory(method="manhattan"))    
}




##' Add a set of configurations to a MultiPattern object	
##'	
##' This plugin adds several configurations that use non-negative matrix factorization.	
##' This function will not change MP if the dataset is not suitable	
##' for NMF, for example if the dataset contains missing values	or too-few columns.	
##'
##' @param MP MultiPattern	
##' @param data.name character	
##' @param config.prefix character	
##' @param preprocess.prefix character	
##' @param preprocess character	
##'	
##' @export	
nmf.MultiPatternPlugin = function(MP, data.name, config.prefix,	
                                  preprocess.prefix="", preprocess=NULL) {	
  
  ## fetch the data and make sure it is positive	
  nowdata = as.matrix(MP$data[[data.name]])	
  
  ## avoid work if the input data contains non-finite elements	
  if (sum(!is.finite(nowdata))>0) {	
    return (MP)	
  }	
  ## shift the data so that it is all > 0	
  if (min(nowdata)<=0) {	
    nowdata = nowdata + abs(min(nowdata)) + MP$settings$nmf.bg	
  }	
  
  ## determine a maximal rank for the nmf (from settings or from data)	
  maxr = min(MP$settings$nmf.rank, ncol(nowdata)-1)	
  if (maxr==0) {        	
    maxr = floor(2*log2(ncol(nowdata)))	
    maxr = min(maxr, ncol(nowdata)-1)	
    if (maxr<2) {	
      return (MP)	
    }	
  } 	
  	
  ## create nmf representations of the data	
  datanmf = list()    	
  for (nowr in seq(2, maxr, by=2)) {	
    nownmf = NMF::nmf(nowdata, nowr)	
    datanmf[[paste0(data.name, ".nmf", nowr)]] = NMF::basis(nownmf)	
  }	
  ## add datasets into MP	
  MP = MPaddData(MP, datanmf)	
  ## add euclidean distance configurations for each dataset	
  for (nowr in seq(2, maxr, by=2)) {	
    MP = MPaddConfig(MP, paste0(config.prefix, data.name,":nmf", nowr),	
                     data.name=paste0(data.name, ".nmf", nowr),	
                     dist.fun=dist.euclidean)            	
  }	
  
  MP	
}




##' Add a set of configurations to a MultiPattern object
##'
##' This plugin transforms an input dataset using PCA, then adds configurations
##' that use the transformed data. The configurations will of two types:
##' individual PC components, e.g. only PC1 or only PC2;
##' pairs of PC components, e.g. PC1, PC2 (PC1.PC2) or PC2, PC4 (PC2.PC4)
##'
##' The maximal PC component for the configurations is determined from MP$settings.
##' 
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
pca.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                  preprocess.prefix="", preprocess=NULL) {
  
  ## fetch number of PCs the user wants to consider
  numPCs = MP$settings$num.PCs
  
  ## get matrices with the all the raw data and with PCA transformed data
  datapca = getPCAsubset(MP$data[[data.name]], subset=numPCs)
  if (is.null(datapca)) {
    return(MP)
  }
  
  ## add transformed data and configurations
  pca.ica.helper(MP, data.name, config.prefix,
                 preprocess.prefix=preprocess.prefix, preprocess=preprocess,
                 data=datapca, type="pca")
}




##' Add a set of configurations to a MultiPattern object
##'
##' This plugin adds configurations that use a pearson-correlation dissimilarity
##' 
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
pearson.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                       preprocess.prefix="", preprocess=NULL) {
  
  MP = MPaddConfig(MP, paste0(config.prefix, data.name,":pearson.D"),
                   data.name=data.name,
                   dist.fun=MPdistFactory(method="pearson", directional=TRUE))
  MP = MPaddConfig(MP, paste0(config.prefix, data.name,":pearson.U"),
                   data.name=data.name,
                   dist.fun=MPdistFactory(method="pearson", directional=FALSE))
  MP
}




##' Add a set of configurations to a MultiPattern object
##'
##' This plugin adds several configurations that use clustering-driven dissimilarities
##' (using pam clustering). Some configurations (marked reg) create dissimilarities
##' that reinforce traditional clustering. Other configuration (marked alt) create
##' dissimilarities that promote alternative clustering. Settings for this plugin
##' are extracted from MP$settings. 
##' 
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
pam.MultiPatternPlugin = function(MP, data.name, config.prefix, preprocess.prefix="",
                                  preprocess=NULL) {
  
  ## get a range of k for the clustering
  clust.k = MP$settings$clust.k    
  clust.k.range = seq(2, abs(clust.k))
  clust.k.range = clust.k.range[clust.k.range>1]
  
  clustconf = paste0(config.prefix, data.name, ":clust.")
  ## add cluster-based distance for each k
  for (nowk in clust.k.range) {
    if (nrow(MP$data[[data.name]])>(2*nowk)+2) {
      ncm.init = "P"                    
      nowdreg = MPdistFactory(method="pam",
                              clust.k=nowk, clust.method="pam", clust.alt=FALSE)
      nowdalt = MPdistFactory(method="pam",
                              clust.k=nowk, clust.method="pam", clust.alt=TRUE)
      MP = MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "reg"),
                  data.name=data.name, dist.fun=nowdreg)
      MP = MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "alt"),
                  data.name=data.name, dist.fun=nowdalt)
    }
  }
  
  MP
}




##' Add a set of configurations to a MultiPattern object
##'
##' This plugin transforms an input dataset using rpca, then defines configurations
##' on the transformed data. The plugin does not do anything if the input data
##' contains non-finite elements. Settings for the plugin are extracted from MP$settings.
##'
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
rpca.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                   preprocess.prefix="", preprocess=NULL) {
  
  ## make sure the input data is all finite
  nowdata = as.matrix(MP$data[[data.name]])
  if (sum(!is.finite(nowdata))>0) {
    return (MP)
  }
  
  term.delta = MP$settings$rpca.term.delta
  
  datarpca = rpca::rpca(nowdata, term.delta=term.delta);
  rownames(datarpca$L) = rownames(datarpca$S) = rownames(MP$data[[data.name]])
  colnames(datarpca$L) = colnames(datarpca$S) = colnames(MP$data[[data.name]])            
  Lname = paste0(data.name, ".L")
  Sname = paste0(data.name, ".S")
  MP = MPaddData(MP, setNames(list(datarpca$L, datarpca$S), c(Lname, Sname)))           
  MP = pca.MultiPatternPlugin(MP, Lname, config.prefix,
                              preprocess.prefix, preprocess)
  MP = pca.MultiPatternPlugin(MP, Sname, config.prefix,
                              preprocess.prefix, preprocess)
  
  MP
}




##' Add a set of configurations to a MultiPattern object
##'
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
subspace1.MultiPatternPlugin = function(MP, data.name, config.prefix, preprocess.prefix="",
                                        preprocess=NULL) {
  
  if (is.null(preprocess)) {
    subspaces = as.list(colnames(MP$data[[data.name]]))
    names(subspaces) = colnames(MP$data[[data.name]])
  } else {
    subspaces = as.list(preprocess)
    names(subspaces) = preprocess
  }
  MP = MPaddConfig(MP, paste0(config.prefix, data.name, preprocess.prefix, ":subspace1"),
                   data.name=data.name, preprocess=subspaces)
  
  MP    
}




##' Add a set of configurations to a MultiPattern object
##'
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
subspace2.MultiPatternPlugin = function(MP, data.name, config.prefix, preprocess.prefix="",
                                        preprocess=NULL) {
  
  if (is.null(preprocess)) {
    subspaces = MPgetAll2Subspaces(colnames(MP$data[[data.name]]))
  } else {
    subspaces = MPgetAll2Subspaces(preprocess)
  }
  MP = MPaddConfig(MP, paste0(config.prefix, data.name, preprocess.prefix, ":subspace2"),
                   data.name=data.name, preprocess=subspaces)
  
  MP
}




##' Add a set of configurations to a MultiPattern object
##'
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
subspacer.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                        preprocess.prefix="", preprocess=NULL) {
  
  Nsub = ceiling(MP$settings$subspace.num.random)
  if (Nsub==0) {
    return (MP)
  }
  Dsub = MP$settings$subspace.d.random
  if (is.null(preprocess)) {
    temp = colnames(MP$data[[data.name]])
  } else {
    temp = preprocess
  }
  if (Dsub<1) {
    Dsub = abs(Dsub*length(temp))
  }            
  Dsub = min(length(temp), ceiling(Dsub))
  subspaces = MPgetRandomSubspaces(temp, Nsub, Dsub)
  MP = MPaddConfig(MP, paste0(config.prefix, data.name,
                              preprocess.prefix, ":subspaceR"),
                   data.name=data.name, preprocess=subspaces,
                   dist.fun=dist.euclidean)
  
  MP
}




##' Add a set of configurations to a MultiPattern object
##'
##' This plugin adds configurations that use a spearman-correlation dissimilarity
##' 
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
spearman.MultiPatternPlugin = function(MP, data.name, config.prefix,
                                       preprocess.prefix="", preprocess=NULL) {
  
  MP = MPaddConfig(MP, paste0(config.prefix, data.name,":spearman.D"),
                   data.name=data.name,
                   dist.fun=MPdistFactory(method="spearman", directional=TRUE))
  MP = MPaddConfig(MP, paste0(config.prefix, data.name,":spearman.U"),
                   data.name=data.name,
                   dist.fun=MPdistFactory(method="spearman", directional=FALSE))
  MP
}




###############################################################################
## Other functions related to plugins


##' helper function to fastica and pca plugins; first few arguments are same
##' as for MultiPatternPlugin functions. Last two 
##'
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##' @param data matrix, data
##' @param type character, one of pca or ica
##' 
pca.ica.helper = function(MP, data.name, config.prefix,
                          preprocess.prefix="", preprocess=NULL,
                          data=NULL, type=c("pca", "ica")) {
  
  type = match.arg(type)
  data.newname = paste0(data.name, ".", type)
  components = colnames(data)
  
  ## add dataset
  MP = MPaddData(MP, setNames(list(data), data.newname))
  
  ## add configurations (these are all euclidean-distance configs)
  all = setNames(list(components), paste0(toupper(type),".",ncol(data)))
  singles = setNames(as.list(components), components)
  twospaces = MPgetAll2Subspaces(components)  
  sets = c(all, singles, twospaces)
  MPaddConfig(MP, paste0(config.prefix, data.name, ":", names(sets)),
              data.name=data.newname, preprocess=sets)
}




##' Get of MultiPattern plugins available 
##'
##' This returns a vector of strings. A string x indicates that there exist
##' an accessible function called x.MultiPatternPlugin.
##'
##' @export
MPlistPlugins = function() {

  ## helper function to extract object that end in .MultiPatternPlugin from a package
  getPlugins = function(pkg=NULL) {
    if (is.null(pkg)) {
      objs = ls(envir=parent.frame(n=2));
    } else {
      objs = ls(paste0("package:", pkg))
    }
    objs[grep(".MultiPatternPlugin", objs)]
  }
  
  ## get a list of plugins in all loaded packages
  ans = list()
  si = sessionInfo()
  for (nowpkg in names(si$otherPkgs)) {
    ans[[nowpkg]] = getPlugins(nowpkg)
  }
  ## also add plugins in non-package
  ans[[length(ans)+1]] = getPlugins()
  
  ## cut out the .MultiPatternPlugin from the function names
  ans = unlist(ans)
  names(ans) = NULL
  ans = sapply(strsplit(ans, "\\."),
               function(x) {
                 paste(x[1:(length(x)-1)], collapse=".")
               })    
  
  sort(ans)
}

