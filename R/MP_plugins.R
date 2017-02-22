## Functions that create MP configurations (plug-ins for MPeasyConfig)
## 
##




##' Add a set of configurations to an MultiPattern object
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
            MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "reg"),
                        data.name=data.name, dist.fun=nowdreg)
            MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "alt"),
                        data.name=data.name, dist.fun=nowdalt)
            rm(nowdreg, nowdalt, ncm.init)                    
        }
    }
    
    MP
}




##' Add a set of configurations to an MultiPattern object
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
hclust.MultiPatternPlugin = function(MP, data.name, config.prefix, preprocess.prefix="",
    preprocess=NULL) {

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
                MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "reg"),
                            data.name=data.name, dist.fun=nowdreg)
                MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "alt"),
                            data.name=data.name, dist.fun=nowdalt)
                rm(nowdreg, nowdalt, ncm.init)
            }
        }
    }
             
    MP    
}



##' Add a set of configurations to an MultiPattern object
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
        temp = as.list(colnames(MP$data[[data.name]]))
        names(temp) = colnames(MP$data[[data.name]])
    } else {
        temp = as.list(preprocess)
        names(temp) = preprocess
    }
    MPaddConfig(MP, paste0(config.prefix, data.name, preprocess.prefix, ":subspace1"),
                data.name=data.name, preprocess=temp)
    
    MP    
}



##' Add a set of configurations to an MultiPattern object
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
        temp = MPgetAll2Subspaces(colnames(MP$data[[data.name]]))
    } else {
        temp = MPgetAll2Subspaces(preprocess)
    }
    MPaddConfig(MP, paste0(config.prefix, data.name, preprocess.prefix, ":subspace2"),
                data.name=data.name, preprocess=temp)
    
    MP
}




##' Add a set of configurations to an MultiPattern object
##'
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
subspacer.MultiPatternPlugin = function(MP, data.name, config.prefix, preprocess.prefix="",
    preprocess=NULL) {

    Nsub = ceiling(MP$settings$subspace.num.random)
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
    temp.subspaces = MPgetRandomSubspaces(temp, Nsub, Dsub)
    MPaddConfig(MP, paste0(config.prefix, data.name, preprocess.prefix, ":subspaceR"),
                data.name=data.name, preprocess=temp.subspaces, dist.fun=dist.euclidean)
    rm(temp.subspaces, Nsub, Dsub)
    
    MP
}




##' Add a set of configurations to an MultiPattern object
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
                data.name=data.name, dist.fun=MPdistFactory(method="euclidean"))    
    MP
}




##' Add a set of configurations to an MultiPattern object
##'
##' This plugin adds one configuration that uses a spearman-correlation dissimilarity
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

    MPaddConfig(MP, paste0(config.prefix, data.name,":spearman"),
                data.name=data.name, dist.fun=MPdistFactory(method="spearman"))    
    MP
}




##' Add a set of configurations to an MultiPattern object
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
                data.name=data.name, dist.fun=MPdistFactory(method="canberra"))    
    MP
}




##' Add a set of configurations to an MultiPattern object
##'
##' This plugin adds one configuration that uses a traditional manhattan distance.
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
                data.name=data.name, dist.fun=MPdistFactory(method="manhattan"))    
    MP
}




##' Add a set of configurations to an MultiPattern object
##'
##' This plugin transforms an input dataset using PCA, then adds configurations
##' that use the transformed data. The configurations will of three types:
##' (i) individual PC components, e.g. only PC1 or only PC2;
##' (ii) series of PC components, e.g. PC1 to PC2 (PC1..PC2) or PC1 to PC3 (PC1..PC3);
##' (iii) pairs of PC components, e.g. PC1, PC2 (PC1.PC2) or PC2, PC4 (PC2.PC4)
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
    
    if (!is.null(datapca)) {
        pcaname = paste0(data.name, ".pca")
        ## add dataset to the MPconfig
        MPaddData(MP, setNames(list(datapca), pcaname)) 
        ## create a series of configurations for pca
        pcasets = list()
        imax = min(ncol(MP$data[[data.name]])-1, ncol(datapca))
        for (i in 1:imax) {
            if (i==1) {
                pcasets[[colnames(datapca)[1]]] = colnames(datapca)[1]                        
            } else if (i>2) {
                pcasets[[paste(colnames(datapca)[c(1,i)], collapse="..")]] = colnames(datapca)[1:i]
            }
        }
        ## add configuration with increasing PCA details
        ## (these are all euclidean-distance configs based on PCA-trasnsformed data)
        MPaddConfig(MP, paste0(config.prefix, data.name, ":", names(pcasets)),
                    data.name=pcaname, preprocess=pcasets)
        ## add configurations with pairwise PCA columns
        twospaces = MPgetAll2Subspaces(head(colnames(datapca)))
        if (length(twospaces)>0) {
            MPaddConfig(MP, paste0(config.prefix, data.name, ":", names(twospaces)),
                        data.name=pcaname, preprocess=twospaces)
        }
    }
    
    MP
}





##' Add a set of configurations to an MultiPattern object
##'
##' This plugin defined random dissimilarity matrices.
##'
##' @param MP MultiPattern
##' @param data.name character
##' @param config.prefix character
##' @param preprocess.prefix character
##' @param preprocess character
##'
##' @export
random.MultiPatternPlugin = function(MP, data.name=NULL, config.prefix="",
    preprocess.prefix="", preprocess=NULL) {
    
    ## find the number of random configurations in the MP settings
    Nrandom = max(0, MP$settings$num.random)
    
    if (Nrandom>0) {
        ## only add if the same configurations have not been added before
        rannames = paste0(config.prefix, "rnorm.", 1:Nrandom)
        if (sum(rannames %in% names(MP$configs))==0) {                                
            MPaddConfig(MP, paste0(config.prefix, "rnorm.", 1:Nrandom), names(MP$data)[1], 
                        preprocess=vector("list", Nrandom), dist.fun=dist.rnorm)
        }
    }
    
    MP
}




##' Add a set of configurations to an MultiPattern object
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
    Lname = paste0(data.name, ".rpcaL")
    Sname = paste0(data.name, ".rpcaS")
    MPaddData(MP, setNames(list(datarpca$L, datarpca$S), c(Lname, Sname)))           
    MP = pca.MultiPatternPlugin(MP, Lname, config.prefix, preprocess.prefix, preprocess)
    MP = pca.MultiPatternPlugin(MP, Sname, config.prefix, preprocess.prefix, preprocess)
    
    MP
}




##' Add a set of configurations to an MultiPattern object
##'
##' This plugin transforms an input dataset using rpca, then defines configurations
##' on the transformed data. The plugin does not do anything if the input data
##' contains non-finite elements. Settings for the plugin are extracted from MP$settings.
##'
##' Note: this function will not change MP if the dataset is not suitable
##' for NMF. This includes situations when the dataset contains missing values
##' or contains too-few columns.
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
    MPaddData(MP, datanmf)
    ## add euclidean distance configurations for each dataset
    for (nowr in seq(2, maxr, by=2)) {
        MPaddConfig(MP, paste0(config.prefix, data.name,":nmf", nowr),
                    data.name=paste0(data.name, ".nmf", nowr),
                    dist.fun=dist.euclidean)            
    }
    
    MP
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
            objs = ls();
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
