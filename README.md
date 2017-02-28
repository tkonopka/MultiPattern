# MultiPattern

**MultiPattern** is an R package for discovery of multiple patterns in data.


## Overview

Unsupervised exploration of data is, in general, an open-ended task. 

Consider the task of clustering the following two toy datasets with points arranged in two dimensions. 

<img src="https://github.com/tkonopka/MultiPattern/blob/master/figures/readme_raw-1.png?raw=true" alt="two toy datasets in two dimensions" width="210px"></img>

In each dataset, there exists a natural grouping for the points. In the first case, there are four natural clusters; in the second case, six clusters. Various machine learning algorithms, for example using hierarchical approaches, would be able to identify these groups.

But let's suppose that we need to assign the points into two or three categories. Although this new task is similar to the first, it does not have a unique solution. Therefore, naive machine learning algorithms like hierarchical clustering that produce a single output cannot provide a complete understanding of the underlying data. 

The **MultiPattern** package provides a framework for exploring data in a systematic manner, including in ambiguous situations as in these examples. For the first dataset, the output of a multi-pattern analysis may be as follows (arranged in two colors).

<img src="https://github.com/tkonopka/MultiPattern/blob/master/figures/readme_k4-1.png?raw=true" alt="multiple patterns in a four-group dataset" width="315px">
</img>

For the second dataset, the output may be as follows (arranged into three colors).

<img src="https://github.com/tkonopka/MultiPattern/blob/master/figures/readme_k6-1.png?raw=true" alt="multiple patterns in a six-group dataset" width="525px">
</img>

Each of the suggested partitions reveals a reasonable pattern. Follow-up analyses can then either exploit one of these patterns, or combine information from many of them at once. 

The package documentation and vignette contains details about the package usage. The vignette also describes applications of alternative clusterings and provides links to relevant literature.



## Development

The package code is available under a GPL-2 license.

The package uses some third-party packages. You should install these on your
R platform before running MultiPattern.

- [Rcssplot](https://github.com/tkonopka/Rcssplot) - customization of graphics functions
- [cluster] - methods for cluster analysis
- [NMF] - non-negative matrix factorization
- [rpca] - robust principal components analysis
- [knitr] - report generation (vignettes)
- [rmarkdown] - report generation (vignettes)

(None of these packages are actually required by the core package. However, these packages may be invoked in customized analyses, during visualization, or when generating vignettes)







