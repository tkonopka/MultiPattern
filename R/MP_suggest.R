## Package: MultiPattern
##
## Function for suggesting configurations for a MultiMetric analysis
##




##' Suggest a set of clustering analyses for a dataset
##'
##' This function can take some time to execute for large datasets. This is
##' because the function applies principal component decompositions to parts of the
##' input data and stores the intermediate results in the MP object. This computation
##' speeds up subsequent steps of the MultiMetric analysis.
##' 
##' @param MP MultiMetric configuration object
##' @param data character. Name of dataset defined in MP.
##' @param verbose logical. Set FALSE to make the function silent;
##' set TRUE to see updates (via message).
##' 
##' @export
MPsuggestConfig = function(MP, data, verbose=TRUE) {
  
  checkArgClass(MP, "MultiPattern")
  checkArgClass(data, "character")
  
  if (length(data)!=1) {
    stop("data must have length 1\n")
  }
  if (!(data %in% names(MP$data))) {
    stop("data is not defined in MP object\n")
  }       
  
  ## get the actual data matrix
  dd = MP$data[[data]]

  if (nrow(dd)<10) {
    stop("data has too-few rows for MPsuggestConfig\n")
  }
  
  if (verbose & object.size(dd)>1e6) {
    message("This may take a little time. Please wait... ")
  }
  
  num.start = length(MP$configs)
  config.prefix = "AUTO:"
  
  
  ## ###############################################################################
  ## Start modifications of MP object with new configurations
  
  ## find which features are real, which are binary, etc.
  feature.class = MPPartitionFeatures(dd)
  MP$auto = feature.class
  
  ## for characters, add hamming
  feature.char = feature.class[["char"]]
  if (length(feature.char)>0) {
    MP = MPeasyConfig(MP, data=data, config.prefix=config.prefix,
                      preprocess.prefix="char",
                      preprocess=feature.char,
                      type="hamming")
  }

  ## remaining features are all numeric
  feature.class[["char"]] = NULL
  
  if (length(feature.class)>0) {
    ## add configurations to MP based on these categories of variables
    distconf = paste0(config.prefix, data, ":", names(feature.class))
    MP = MPaddConfig(MP, paste0(distconf, ":euclidean"),
                     data, preprocess=feature.class, dist.fun=dist.euclidean)
    MP = MPaddConfig(MP, paste0(distconf, ":canberra"),
                     data, preprocess=feature.class, dist.fun=dist.canberra)
    rm(distconf)
  }
  
  ## for real-valued data, add pca 
  for (nowf in intersect(c("real", "realskew", "bin", "binskew"), names(feature.class))) {
    ## for pca and rpca, avoid features with NAs
    nowdd = dd[, feature.class[[nowf]], drop=FALSE]
    nowdd.ok = apply(nowdd, 2, function(x) {sum(!is.finite(x))==0})
    nowdd = nowdd[, nowdd.ok, drop=FALSE]
    if (ncol(nowdd)>0) {            
      data.name = paste0(data, ":", nowf)
      MP = MPaddData(MP, setNames(list(nowdd), data.name))        
      MP = MPeasyConfig(MP, data=data.name,
                        config.prefix=config.prefix, type=c("pca", "ica"))
      ## here can remove the temporary datasets (pca plugin
      ## would have by now created their own helper tables)
      MP = MPremove(MP, data=data.name)                
    }
  }
  
  ## for all the feature types, add clust-based distances
  okf = intersect(c("bin", "binskew", "multi", "multiskew", "real", "realskew"),
                  names(feature.class))
  for (nowf in okf) {
    nowfeatures = feature.class[[nowf]]
    ## add hclust distances
    clustconf = paste0(config.prefix, data, ":", nowf, ":clust.")
    for (nowk in c(2, 3)) {
      ## add clust-based distances for this feature class
      for (ncm in c("average", "single", "complete", "pam")) {
        if (ncm %in% c("average", "single", "complete")) {
          nowdreg = MPdistFactory(method="hclust",
                                  clust.k=nowk, clust.method=ncm, clust.alt=FALSE)
          nowdalt = MPdistFactory(method="hclust",
                                  clust.k=nowk, clust.method=ncm, clust.alt=TRUE)
          ncm.init = toupper(substring(ncm, 1, 1))
        } else if (ncm == "pam") {
          nowdreg = MPdistFactory(method="pam",
                                  clust.k=nowk, clust.method="pam", clust.alt=FALSE)
          nowdalt = MPdistFactory(method="pam",
                                  clust.k=nowk, clust.method="pam", clust.alt=TRUE)
          ncm.init = "P"
        }
        MP = MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "reg"),
                         data.name=data, preprocess=nowfeatures,
                         dist.fun=nowdreg)
        MP = MPaddConfig(MP, paste0(clustconf, ncm.init, nowk, "alt"),
                         data.name=data, preprocess=nowfeatures,
                         dist.fun=nowdalt)
      }                                
    }
    
    if (length(nowfeatures)>2) {
      MP = MPeasyConfig(MP, data=data, config.prefix=config.prefix,
                        preprocess.prefix=nowf, preprocess=nowfeatures,
                        type="subspaceR")
    }
    
  }
  
  ## finish with some updates to the user
  if (verbose & object.size(dd)>1e6) {
    message("done")
  }    
  if (verbose) {
    num.end = length(MP$configs)
    message(paste0("MPsuggestConfig created ", num.end-num.start,
                   " configurations"))
  }
  
  MP
}



##' Partition features in a data matrix into categories like "real",
##' "realskew", "bin", "multi", etc.
##'
##' @param data matrix of data
##'
##' @return named list, names indicate feature class, contents are features
##'
MPPartitionFeatures = function(data) {

  ## helper function to compute skew and kurtosis from vectors
  mystats = function(x) {
    x = x[is.finite(x)]
    n = length(x)
    xmean = mean(x)
    mu2 = sum((x-xmean)^2)
    mu3 = sum((x-xmean)^3)
    mu4 = sum((x-xmean)^4)
    c(skew=n^(1/2)*mu3/mu2^(3/2), kurtosis=n*mu4/(mu2^2))
  } 

  ## identify all character/factor columns
  result = list()
  if ("data.frame" %in% class(data)) {
    result$char = colnames(data)[sapply(data, class) %in% c("character", "factor")]
  } else {
    if (class(data[,1]) %in% c("character", "factor")) {
      result$char = colnames(data)
    }
  }
  
  if (length(result[["char"]]) == ncol(data)) {
    return(result)
  }
  
  ## create an ad-hoc classification of features by data range (integer/real, etc)
  feature.types = data.frame(Column=setdiff(colnames(data), result$character),
                             skew=0, ex.kurtosis=0, n.unique=0, stringsAsFactors=F)
  n.u = "n.unique"
  ft.columns = c("skew", "ex.kurtosis", "n.unique")
  rownames(feature.types) = feature.types[, "Column"]
  for (nowcol in rownames(feature.types)) {
    nowdata = as.numeric(data[,nowcol])
    nowu = length(unique(nowdata))
    nowstats = mystats(nowdata)
    feature.types[nowcol, ft.columns] = c(nowstats[1], nowstats[2]-3, nowu)
  }
  
  skewth = max(1, 0.5*log2(ncol(data)))
  ## make feature.class a vector with a word description of the feature type
  feature.class = rep("real", nrow(feature.types))
  names(feature.class) = rownames(feature.types)
  feature.class[feature.types[, n.u]==1] = "single"
  feature.class[feature.types[, n.u]==2] = "bin"
  feature.class[abs(feature.types[, "skew"])>skewth &
                  feature.types[, n.u]==2] = "binskew"
  feature.class[feature.types[, n.u]<nrow(data)/2 &
                  feature.types[, n.u]>2] = "multi"
  feature.class[feature.types[, n.u]<nrow(data)/2 &
                  feature.types[, n.u]>2 &
                  abs(feature.types[, "skew"])>skewth] = "multiskew"
  feature.class[feature.types[, n.u]>nrow(data)/2 &
                  abs(feature.types[, "skew"])>skewth] = "realskew"     

  ## avoid splits of features with just a few features
  feature.table = table(feature.class)
  for (ii in c("bin", "multi", "real")) {
    iiskew = paste0(ii, "skew")
    if (ii %in% names(feature.table) & iiskew %in% names(feature.table)) {
      iifrac = feature.table[ii]/(feature.table[ii]+feature.table[iiskew])
      if (iifrac<0.01 | feature.table[ii]<3) {
        feature.class[feature.class %in% c(ii, iiskew)] = iiskew
      } else if (iifrac>0.99 | feature.table[iiskew]<3) {
        feature.class[feature.class %in% c(ii, iiskew)] = ii
      }
    }
  }    
  
  ## organize the features into a list by data type
  result = c(result, split(names(feature.class), feature.class))
  ## avoid working with features with a single value
  result[["single"]] = NULL
  if (length(result)==0) {
    stop("data does not appear to have distinct values\n")
  }
  
  result
}


