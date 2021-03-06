% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MP_helpers.R
\name{MPrandomizeMatrix}
\alias{MPrandomizeMatrix}
\title{Permute values in a matrix}
\usage{
MPrandomizeMatrix(x, perm.method = c("shuffle", "bootstrap"))
}
\arguments{
\item{x}{input matrix}

\item{perm.method}{character, determines how the values are permuted

The randomization methods are as follows.
shuffle - re-arranges values in a matrix using sampling *without* replacement
bootstrap - re-arranges values in a matrix using sampleing *with* replacement}
}
\description{
Permute values in a matrix
}
