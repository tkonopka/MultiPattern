% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MP_suggest.R
\name{MPsuggestConfig}
\alias{MPsuggestConfig}
\title{Suggest a set of clustering analyses for a dataset}
\usage{
MPsuggestConfig(MP, data, verbose = TRUE)
}
\arguments{
\item{MP}{MultiMetric configuration object}

\item{data}{character. Name of dataset defined in MP.}

\item{verbose}{logical. Set FALSE to make the function silent;
set TRUE to see updates (via message).}
}
\description{
This function can take some time to execute for large datasets. This is
because the function applies principal component decompositions to parts of the
input data and stores the intermediate results in the MP object. This computation
speeds up subsequent steps of the MultiMetric analysis.
}
