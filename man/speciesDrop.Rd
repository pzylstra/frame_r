% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ensembles.R
\name{speciesDrop}
\alias{speciesDrop}
\title{Removes identified species from a param file and tidies the file}
\usage{
speciesDrop(base.params)
}
\arguments{
\item{base.params}{The params file to be adjusted}
}
\description{
Species where liveLeafMoisture == 100 are removed, then
empty strata removed and species & strata renumbered consecutively
}
