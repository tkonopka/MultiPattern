% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MP_plugins.R
\name{ica.MultiPatternPlugin}
\alias{ica.MultiPatternPlugin}
\title{Add a set of configurations to a MultiPattern object}
\usage{
ica.MultiPatternPlugin(MP, data.name, config.prefix,
  preprocess.prefix = "", preprocess = NULL)
}
\arguments{
\item{MP}{MultiPattern}

\item{data.name}{character}

\item{config.prefix}{character}

\item{preprocess.prefix}{character}

\item{preprocess}{character}
}
\description{
This plugin uses fastICA to identify interesting components, then
adds configurations that uses those components. The configurations are of two types:
individual independent components, e.g. only IC1 or only IC2
pairs of IC coponents, e.g. IC1, IC2 (IC1.IC2) or IC2, IC4 (IC2.IC4)
}
