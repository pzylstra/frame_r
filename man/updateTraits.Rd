% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Building.R
\name{updateTraits}
\alias{updateTraits}
\title{Checks traits table for column names and adds data from a new table}
\usage{
updateTraits(traits, traitsNew, deleteReplicates = TRUE, fill = TRUE)
}
\arguments{
\item{traits}{Existing trait table}

\item{traitsNew}{Table containing new trait data}

\item{deleteReplicates}{If set to TRUE, retains the first mention of a species, then removes all following mentions}

\item{fill}{If set to TRUE, fills empty numeric values with mean values of the genus}
}
\value{
table
}
\description{
Checks traits table for column names and adds data from a new table
}
\examples{
Traits <- ausTraitTable(version = "3.0.2", path = "data/austraits")
TraitsNew <- read.csv("Traits.csv")
Tr <- updateTraits(traits = Traits, traitsNew = TraitsNew, deleteReplicates = TRUE, fill = TRUE)

}
