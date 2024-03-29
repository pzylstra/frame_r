% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ensembles.R
\name{weatherSet_Frame}
\alias{weatherSet_Frame}
\title{Updates parameter files with weather from a dataset,
then models fire from non-deterministic plant parameters
using plantVarFrame to modify individual species with their own measured variability
Differs from weatherSetS by using  FRaME formatted tables}
\usage{
weatherSet_Frame(
  base.params,
  weather,
  Structure,
  Flora,
  a,
  db.path = "out_mc.db",
  jitters = 10,
  l = 0.1,
  Ms = 0.01,
  Pm = 1,
  Mr = 1.001,
  updateProgress = NULL
)
}
\arguments{
\item{base.params}{Parameter input table}

\item{weather}{A dataframe with the four fields:
tm - Sequential numbering of the records
T - Air temperature (deg C)
W - Wind velocity (km/h)
DFMC - Dead fuel moisture content (proportion ODW)}

\item{Structure}{A dataframe with the fields:
record - a unique, consecutively numbered identifier per site
site - a unique identifier per site
NS, El, Mid & Can - the mean separation between plants (m) per stratum
ns_e, ns_m, e_m, e_c, m_c - Logical field indicating whether plants in the stratum
on the left grow directly beneath those in the stratum on the right. Acceptable values
are t, f, or blank, where the outcome will be decided by the relative stratum heights.
nsR, eR, mR, cR - maximum species richness recorded for each stratum}

\item{Flora}{A dataframe with the fields:
record - a unique, consecutively numbered identifier per site
species - the name of the species, which will call trait data from 'default.species.params'
moisture - the moisture content of the species in whole numbers (eg 1 for 100% ODW)
stratum - numeric value from 1 to 4, counting from lowest stratum
comp - % composition or count of that species in the stratum. If absent, all species will be considered equally
hc, he, ht, hp & w - canopy dimensions for that species (m)
clump - mean ratio of clump diameter to crown diameter
openness - proportion of plant canopy occupied by gaps between clumps}

\item{a}{A unique identifier for the record being run}

\item{db.path}{Name of the exported database}

\item{jitters}{Number of repetitions}

\item{l}{Variation around input leaf dimensions}

\item{Ms}{Standard deviation of LFMC}

\item{Pm}{Multiplier of mean LFMC}

\item{Mr}{Truncates LFMC variability by +/- Mr * LFMC}

\item{updateProgress}{Progress bar for use in the dashboard}
}
\value{
dataframe
}
\description{
Updates parameter files with weather from a dataset,
then models fire from non-deterministic plant parameters
using plantVarFrame to modify individual species with their own measured variability
Differs from weatherSetS by using  FRaME formatted tables
}
