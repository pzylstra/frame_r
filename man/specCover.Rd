% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flora_dynamics.R
\name{specCover}
\alias{specCover}
\title{Finds % cover of surveyed Species and groups minor Species}
\usage{
specCover(dat, thres = 5, pnts = 10)
}
\arguments{
\item{dat}{The dataframe containing the input data,}

\item{thres}{The minimum percent cover (0-100) of a Species that will be kept single}

\item{pnts}{The number of points measured in a transect}
}
\value{
dataframe
}
\description{
Input table requires the following fields:
Point - numbered point in a transect
Species - name of the surveyed Species
Age - age of the site since the triggering disturbance
}
\details{
Species that are less common than the set threshold are combined as "Minor Species"
}
