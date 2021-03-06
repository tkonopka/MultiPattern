% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MP_distances.R
\name{dist.euclidean}
\alias{dist.euclidean}
\title{Euclidean distance}
\usage{
dist.euclidean(x)
}
\arguments{
\item{x}{numeric matrix (see stats::dist for details)}
}
\description{
This is a wrapper for stats::dist with method="euclidean".
The wrapper replaces NA and Inf values with mean and max guesses.
}
