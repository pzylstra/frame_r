#' Valid parameters and their default units.
#' 
#' This data frame contains all valid parameter labels with their corresponding 
#' section (site, stratum or species) and default units.
#' A value of \code{NA} for units indicates that a parameter is dimensionless 
#' (e.g. a proportional quantity or a text value).
#'
#' @format A data frame with columns: section (values site,
#' stratum, species); param (standard parameter label);
#' units (unit abbreviation or \code{NA} if dimensionless).
#'
"ParamInfo"
