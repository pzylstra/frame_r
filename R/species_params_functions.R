#' Matches whole or partial species names to DefaultSpeciesParams.
#'
#' For each name in the input character vector \code{names}, finds
#' the index or indices of matching rows in the table of default species
#' parameter values.
#'
#' @param names The names, whole or partial, to search for.
#'
#' @return a named \code{list} where names are the search strings and
#'   each value is the index or indices of matching species (NA for no match).
#'
#' @examples
#' \dontrun{
#' # Look for Acacia dealbata and Cassinia aculeata.
#' # The search is case insensitive.
#' index <- findSpeciesByName(c("Ac deal", "cass ac"))
#'
#' # Search for a term that matches multiple species. The
#' # list element for "cass" will be a vector of indices of
#' # all matching taxa:
#' findSpeciesByName("cass")
#'
#' # Search for a name that does not exist in the default
#' # species parameters table. This results in \code{NA} being
#' # returned for the index:
#' findSpeciesByName("aspidistra")
#' }
#'
#' @seealso \code{\link{getSpeciesParams}}, \link{DefaultSpeciesParams}
#'
#' @export
#'
findSpecies <- function(names) {

  mkregex <- function(s) {
    parts <- stringr::str_split(s, "\\s")[[1]]
    parts <- stringr::str_c(parts, ".*", sep="")
    stringr::regex(str_c(parts, collapse=""), ignore_case = TRUE)
  }

  finder <- function(ptn) {
    matches <- stringr::str_detect(DefaultSpeciesParams$name, ptn)
    if (any(matches)) which(matches)
    else NA
  }

  ptns <- lapply(names, mkregex)
  res <- lapply(ptns, finder)
  names(res) <- names

  res
}


#' Checks if species names are defined in the default species parameters table.
#'
#' For each name (whole or partial) checks whether the name matches one, and only
#' one, record in the default species parameters table.
#'
#' @param names The names, whole or partial, to search for.
#'
#' @return A named boolean vector where names are search terms.
#'
#' @export
#'
isSpeciesKnown <- function(names) {
  sapply(findSpecies(names),
         function(indices) length(indices) == 1 && !is.na(indices))
}


#' Gets default parameter values for one or more species.
#'
#' Takes a character vector of one or more species names (whole or partial)
#' and retrieves the default parameter values for each. If one or more
#' names is not found, or matches mutliple species in the default
#' parameters table, an error message is given.
#'
#' @param names The names, whole or partial, to search for.
#'
#' @return A data frame with default parameter values for each matched species.
#'
#' @examples
#' \dontrun{
#' # Get default parameters for selected species:
#' spp <- c("Poa lab", "D repens", "Hyd laxiflora")
#' params <- getSpeciesParams(spp)
#' }
#'
#' @seealso \code{\link{findSpecies}}, \link{DefaultSpeciesParams}
#'
#' @export
#'
getSpeciesParams <- function(names) {
  indices <- findSpecies(names)

  # check for non-matches and multiple matches
  nas <- sapply(indices, anyNA)
  multis <- sapply(indices, function(ii) length(ii) > 1)

  if (any(nas) || any(multis)) {
    msg <- NULL

    if (any(nas)) {
      err.names <- paste(names[nas], collapse = ", ")
      msg <- c(msg, paste("The following do not match any species:", err.names))
    }

    if (any(multis)) {
      err.names <- paste(names[multis], collapse = ", ")
      msg <- c(msg, paste("The following match multiple species:", err.names))
    }

    stop(msg)
  }

  # convert from list of indices to matrix
  indices <- do.call(rbind, indices)

  # return params
  DefaultSpeciesParams[indices[, 1], ]
}
