# MultiPattern

**MultiPattern** is an R package for discovery of multiple patterns in data.


## Overview

Unsupervised exploration of data is, in general, an open-ended task. Take clustering, for example. In some situations, the natural grouping of elements in a dataset may be straightforward. In such cases, clustering can be seen as a computational task aiming to identify this particular natural pattern. However, there are situations when the natural grouping either does not exists or it is ambiguous. The role of unsupervised analysis in these cases is to provide hints for further investigation. 

The **MultiPattern** package provides a framework to perform such exploration in a systematic manner. As an example, consider two toy datasets with points arranged in two dimensions. 

<img src="" alt="two toy datasets in two dimensions" width="210px"></img>

The points in the left and right panels are arranged in four and six groups, respectively. But let's suppose that we need to group them in a smaller number of categories. Intuitively, there is not a unique or best way to perform this partitioning. 

A multi-pattern analysis yields a set of alternative, distinct, candidate clusterings. For the first dataset, suggested candidate groupings may be as follows (arranged into two colors). 

<img src="" alt="multiple patterns in a four-group dataset" width="315px">
</img>

For the second dataset, candidates may be as follows (arranged into three colors).

<img src="" alt="multiple patterns in a six-group dataset" width="525px">
</img>

Each of these suggestions reveals a reasonable pattern. Follow-up analyses can then either exploit one of these patterns, or combine information from many of them at once. 

The package documentation and vignette contains details about the package scope, usage, and links to literature on alternative clusterings.



## Development

The package code is available under a GPL-2 license.

The package uses some third-party packages. You should install these on your
R platform before running MultiPattern.

- [Rcssplot](https://github.com/tkonopka/Rcssplot) - customization of graphics functions. 









