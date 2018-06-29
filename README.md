# MultiPattern

![Status](https://travis-ci.org/tkonopka/MultiPattern.svg?branch=master)
[![codecov](https://codecov.io/gh/tkonopka/MultiPattern/branch/master/graph/badge.svg)](https://codecov.io/gh/tkonopka/MultiPattern)


MultiPattern is an R package for discovery of multiple patterns in data.


## Overview

Unsupervised exploration of data is an open-ended task. Consider, for example, the task of clustering the following two toy datasets in two dimensions. 

<img src="https://github.com/tkonopka/MultiPattern/blob/master/figures/readme_raw-1.png?raw=true" alt="two toy datasets in two dimensions" width="210px"></img>

In each dataset, there exists a natural grouping for the points. There are four natural clusters in the first case and six in the second. Various machine learning algorithms, for example using hierarchical approaches, would be able to identify these groups.

But let's suppose that we need to assign the points to a maximum of two or three categories. Although this is also a clustering task, it does not have a unique or a natural solution. Machine learning algorithms like hierarchical clustering that produce a single output can provide one suggestion but inherently cannot provide a complete assessment to the task at hand. 

The MultiPattern package provides a framework for exploring data in a systematic manner and is particularly suited to describing ambiguous situations as in these examples. For the first dataset, the output of a multi-pattern analysis may be series of clusterings below (arranged using two colors). 

<img src="https://github.com/tkonopka/MultiPattern/blob/master/figures/readme_k4-1.png?raw=true" alt="multiple patterns in a four-group dataset" width="315px">
</img>

For the second dataset, the output may be as follows (arranged into three colors).

<img src="https://github.com/tkonopka/MultiPattern/blob/master/figures/readme_k6-1.png?raw=true" alt="multiple patterns in a six-group dataset" width="525px">
</img>

Each of the suggested partitions reveals a reasonable pattern. Follow-up analyses can then either exploit one of these patterns or integrate information from several of them. 


## Documentation

The package documentation and vignette contains details about the package usage. 



## Development

The package code is available under a GPL-2 license.

Parts of the package use third-party components. These pacakges should be installed before running MultiPattern.

- [Rcssplot](https://github.com/tkonopka/Rcssplot) - styling for graphics
- [cluster] - methods for cluster analysis
- [dbscan] - density-based clustering
- [fastICA] - independent-component analysis
- [NMF] - non-negative matrix factorization
- [rpca] - robust principal components analysis
- [knitr] - report generation (vignettes)
- [rmarkdown] - report generation (vignettes)
- [umap] - visualization of meta-map


