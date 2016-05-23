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
#' @export
#'
findSpecies <- function(names) {

  mkregex <- function(s) {
    parts <- stringr::str_split(s, "\\s")[[1]]
    parts <- stringr::str_c(parts, ".*", sep="")
    stringr::regex( stringr::str_c(parts, collapse=""), ignore_case = TRUE )
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


#' Completes a table by adding default parameters for each species as required.
#' 
#' For each species in the input parameter table, this function checks which
#' parameters have been provided and retrieves defaults for missing parameters
#' from the default lookup table: \link{DefaultSpeciesParams}.
#' 
#' If a species has missing parameters but is not present in the table of
#' default lookup table, an error is thrown.
#' 
#' @param tbl The input parameter table
#' 
#' @return The completed parameter table
#' 
#' @export
#' 
completeParams <- function(tbl) {
  if (!.is_param_table(tbl))
    stop("Input table must be a character matrix or data frame with 4 or 5 columns")
  
  tbl <- .as_str_data_frame(tbl)
  
  # Retrieve default parameters as required for a given species
  # and return as additional records to add to the table
  do_species <- function(species.id) {
    recs <- dplyr::filter(tbl, species == species.id)
    stratum <- recs$stratum[1]
    species.name <- dplyr::filter(recs, param == "name")$value
    
    provided <- recs$param
    required <- setdiff( colnames(DefaultSpeciesParams)[-1], provided )
    
    if (length(required == 0)) {
      # Don't need to add any parameters, so return NULL
      NULL
    }
    else if (isSpeciesKnown(species.name)) {
      # Parameters required and species is known
      #
      dat <- getSpeciesParams(species.name) %>% select_(.dots = required)
      nrecs <- ncol(dat)
      
      # return required addition parameter table records
      data.frame(
        stratum = rep(stratum, nrecs),
        species = rep(species.id, nrecs),
        param = colnames(dat),
        value = as.character(dat[1, ]),
        stringsAsFactors = FALSE
      )
      
    } else {
      stop(species.id, " is missing required parameters and is not in the defaults table")
    }
    
  }
  
  ids <- .get_species_ids(tbl)
  new.recs <- do.call(rbind, lapply(ids, do_species))
  
  if (nrow(new.recs) == 0) {
    # no params added, return input
    tbl
    
  } else {
    # add a units column if one was present in the input table
    if (ncol(tbl) == 5) {
      units <- ParamInfo %>%
        dplyr::filter(section == "species") %>%
        select(param, units)
      
      new.recs <- left_join(new.recs, units, by = "param")
    }
    
    tbl <- rbind( tbl, new.recs )
    arrange(tbl, stratum, species)
  }
}

#' Checks that a parameter label is valid.
#' 
#' A parameter label is valid if it is present in the \code{\link{ParamInfo}}
#' table. The comparison ignores case and white-space. If a non-NULL value is
#' supplied for \code{section} the check is restricted to parameter labels
#' in that section (site, stratum or species).
#' 
#' @param label The parameter label.
#' 
#' @param section If not \code{NULL}, one of site, stratum or species.
#' 
#' @param single If \code{TRUE} (default) require that label only
#'   matches a single parameter; if \code{FALSE } allow multiple matches.
#' 
#' @return \code{TRUE} is the label is valid.
#' 
#' @export
#' 
ffm_valid_param <- function(label, section=NULL, single) {
  !is.null(.match_param(label, section, no.match.error = FALSE, single = single))
}


#' Retrieves information for a parameter label.
#' 
#' Searches for a row in the \code{\link{ParamInfo}} table with a matching
#' parameter label. The comparison ignores case and white-space. 
#' If a non-NULL value is supplied for \code{section} the check is restricted
#' to parameter labels in that section (site, stratum or species).
#' 
#' @param label The parameter label.
#' 
#' @param section If not \code{NULL}, one of site, stratum or species.
#' 
#' @param no.match.error If \code{TRUE} an error results when a unique
#'   match is not found for the input label; if \code{FALSE} (default)
#'   the function returns NULL when there is no unique match.
#' 
#' @return Parameter info as a single row data frame; or \code{NULL}
#'   if no matching label was found and \code{no.match.error} is
#'   \code{FALSE}.
#'   
#' @export
#'
ffm_param_info <- function(label, section = NULL, no.match.error = FALSE) {
  i <- .match_param(label, section, no.match.error, single = TRUE)

  if (!is.null(i)) ParamInfo[i, ]
  else NULL
}

#' Returns a summary of the given parameter table.
#' 
#' @export
#' 
ffm_param_table_summary <- function(tbl) {
  u <- function(xs) unique(xs[!is.na(xs)])
  nu <- function(xs) length(u(xs))
  
  # Species IDs within each stratum
  strata.spp <- tapply(tbl$species, 
                       tbl$stratum, 
                       function(ids) u(ids))
  
  names(strata.spp) <- u(tbl$stratum)
  
  list(
    nstrata = nu(tbl$stratum),
    nspecies = nu(tbl$species),
    strata.species = strata.spp)
}

#' Sets the value of a site meta-data parameter.
#' 
#' @param tbl The parameter table in which to set the value.
#' @param param (character) Label of the parameter to set.
#' @param value Value of the parameter (will be converted to character).
#' @param units (option) The units of measurement for the supplied value.
#' 
#' @return The updated parameter table.
#' 
#' @export
#' 
ffm_set_site_param <- function(tbl, param, value, units = NA_character_) {
  .set_param(tbl, stratum.id=NA, species.id=NA, param, value, units)
}


#' Sets the value of a stratum meta-data parameter.
#' 
#' @param tbl The parameter table in which to set the value.
#' @param stratum.id (character) Stratum identifier.
#' @param param (character) Label of the parameter to set.
#' @param value Value of the parameter (will be converted to character).
#' @param units (option) The units of measurement for the supplied value.
#' 
#' @return The updated parameter table.
#' 
#' @export
#' 
ffm_set_stratum_param <- function(tbl, stratum.id, 
                                  param, value, units = NA_character_) {
  .set_param(tbl, stratum.id, species.id=NA, param, value, units)
}


#' Sets the value of a species parameter.
#' 
#' @param tbl The parameter table in which to set the value.
#' @param stratum.id (character) Stratum identifier.
#' @param species.id (character) Species identifier.
#' @param param (character) Label of the parameter to set.
#' @param value Value of the parameter (will be converted to character).
#' @param units (option) The units of measurement for the supplied value.
#' 
#' @return The updated parameter table.
#' 
#' @export
#' 
ffm_set_species_param <- function(tbl, stratum.id, species.id, 
                                  param, value, units = NA_character_) {
  .set_param(tbl, stratum.id, species.id, param, value, units)
}


############################################################################
#
# Non-exported helper functions
#
############################################################################

#' Finds parameter information for a given parameter label.
#' 
#' Searches for a row or rows in the \code{\link{ParamInfo}} table with a matching
#' parameter label. The comparison ignores case and white-space. 
#' If a non-NULL value is supplied for \code{section} the check is restricted
#' to parameter labels in that section (site, stratum or species).
#' 
#' @param param The parameter label.
#' 
#' @param section If not \code{NULL}, one of site, stratum or species.
#' 
#' @param no.match.error If \code{TRUE} an error results when no
#'   match is found, or more than one match is found and \code{single}
#'   is \code{TRUE}; if \code{FALSE} (default) the function returns 
#'   \code{NULL} if unsuccessful.
#'   
#' @param single If \code{TRUE} (default) require that label only
#'   matches a single parameter; if \code{FALSE } allow multiple matches.
#' 
#' @return The index or indices (if \code{single} is \code{FALSE}) of the 
#'   matching row(s) in \code{ParamInfo}; otherwise \code{NULL}
#'   if unsuccessful and \code{no.match.error} is \code{FALSE}.
#'
.match_param <- function(param, section, no.match.error = FALSE, single = TRUE) {
  if (is.null(section))
    labels <- ParamInfo$param
  else
    labels <- dplyr::filter(ParamInfo, section == section)$param
  
  ii <- stringr::str_detect(labels, .make_ptn(param))
  n <- sum(ii)
  
  if (single) {
    if (n == 1) which(ii)
    else if (no.match.error) {
      if (n == 0) stop(param, " not recognized as a parameter label")
      else stop(param, " matches more than one parameter")
    }
    else NULL
  } 
  else {  # not limited to single
    if (n == 0 && no.match.error)
      stop(param, " not recognized as a parameter label")
    else if (n > 0)
      which(ii)
    else NULL
  }
}


#' Updates the value of a parameter in a parameter table.
#' 
#' @param tbl The parameter table in which to set the value.
#' @param stratum.id Stratum identifier.
#' @param species.id Species identifier.
#' @param param Label of the parameter to set.
#' @param value Value of the parameter (will be converted to character).
#' @param units (option) The units of measurement for the supplied value.
#' 
#' @return The updated parameter table.
#' 
.set_param <- function(tbl, stratum.id, species.id,
                       param, value, units = NA_character_) {
  
  section <- 
    if (is.na(stratum.id)) "site"
    else if (is.na(species.id)) "stratum"
    else "species"
  
  i <- .match_param(param, section, no.match.error = TRUE, single = TRUE)
  std.param <- ParamInfo[i, "param"]

  in.scope <- .match_ids(tbl, stratum.id, species.id)
    
  rows <- which( tbl$param == std.param & in.scope ) 
  if (length(rows) == 0) 
    stop("Input table does not contain parameter ", 
         std.param,
         " (", section, ") ",
         "for stratum=", stratum.id,
         ", species=", species.id)
  
  tbl[rows, "value"] <- as.character(value)
  if ( !(is.na(units) | is.null(units)) ) tbl[rows, "units"] <- units
  tbl
}

.match_ids <- function(tbl, stratum.id, species.id) {
  if (is.na(stratum.id)) {
    if (is.na(species.id)) .match_stratum(tbl, stratum.id)
    else stop("Missing stratum ID with non-missing species ID is not a valid combination")
  }
  else .match_stratum(tbl, stratum.id) & .match_species(tbl, species.id)
}

.match_stratum <- function(tbl, stratum.id) {
  if (is.na(stratum.id)) is.na(tbl$stratum)
  else {
    id <- as.character(stratum.id)
    ids <- as.character(tbl$stratum)
    !is.na(ids) & ids == id
  }
}

.match_species <- function(tbl, species.id) {
  # match species NA to stratum meta-data rows
  if (is.na(species.id)) !is.na(tbl$stratum) & is.na(tbl$species)
  else {
    id <- as.character(species.id)
    ids <- as.character(tbl$species)
    !is.na(ids) & ids == id
  }
}

