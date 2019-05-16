#' Checks that a parameter table is complete
#' 
#' A parameter table is considered complete if it contains all parameters specified
#' in the \code{\link{ParamInfo}} table other than the optional parameter
#' \code{propSilicaFreeAsh}. Parameter name comparisons are done ignoring case.
#' 
#' @param tbl The input simulation parameter table (data frame).
#' 
#' @param quiet Set to TRUE to suppress details of missing parameters being
#'   written to the console.
#' 
#' @return \code{TRUE} is all required parameters are present; \code{FALSE}
#'   otherwise.
#'
#' @importFrom dplyr %>% group_by mutate ungroup
#'  
#' @export
#' 
ffm_check_params <- function(params, quiet = FALSE) {
  # optional species parameter
  Silica <- tolower("propSilicaFreeAsh")
  
  required.cols <- c("stratum", "species", "param", "value")
  nrequired <- length(required.cols)
  
  optional.cols <- "units"
  
  given.cols <- colnames(params)
  missing.cols <- setdiff(required.cols, given.cols)
  
  # Just in case
  params <- ungroup(params)
  
  if (length(missing.cols) > 0) {
    if (!quiet) warning("Missing required column(s): ", missing.cols)
    return(FALSE)
  }
    
  extra.cols <- setdiff(given.cols, c(required.cols, optional.cols))
  if (length(extra.cols) > 0) {
    if (!quiet) warning("Extra column(s) present: ", extra.cols)
    return(FALSE)
  }
  
  # ensure correct column order
  ok <- all(given.cols[1:nrequired] == required.cols)
  if (!ok) {
    if (!quiet) warning("Incorrect column order. Should be ",
                        paste(required.cols, collapse = ", "),
                        " optionally followed by units")
    return(FALSE)
  }
  
  
  # Parameter checks are done in lower case
  ParamInfo <- frame::ParamInfo %>%
    mutate(param = tolower(param))
  
  params <- params %>%
    mutate(param = tolower(param))
  
  
  # Check site-level parameters
  given <- params %>%
    dplyr::filter(is.na(stratum))
  
  expected <- frame::ParamInfo %>%
    dplyr::filter(section == "site")
  
  missing <- setdiff(expected$param, given$param)
  if (length(missing) > 0) {
    if (!quiet) warning("Missing site-level parameter(s): ", missing)
    return(FALSE)
  }
  
  extra <- setdiff(given$param, expected$param)
  if (length(extra) > 0) {
    if (!quiet) warning("Unrecognized site-level parameter(s): ", extra)
    return(FALSE)
  }
  
  # Check strata- and species-level parameters
  strata <- sort(unique(params$stratum), na.last = NA)
  
  strata.expected <- {
    x <- dplyr::filter(ParamInfo, section == "stratum")
    x$param
  }
  
  species.expected <- {
    x <- dplyr::filter(ParamInfo, section == "species")
    x$param
  }
  
  for (stratum in strata) {
    sdat <- dplyr::filter(params, stratum == stratum)

    # stratum level
    given <- dplyr::filter(sdat, is.na(species))    
    missing <- setdiff(strata.expected, given$param)
    if (length(missing) > 0) {
      if (!quiet) warning("Stratum ", stratum, " missing parameter(s): ", missing)
      return(FALSE)
    }
    
    # species level
    given <- dplyr::filter(sdat, !is.na(species))
    spp <- sort(unique(given$species))
    
    for (isp in spp) {
      givensp <- dplyr::filter(given, species == isp)
      missing <- setdiff(species.expected, c(given$param, Silica))
      if (length(missing) > 0) {
        if (!quiet) warning("Species ", isp, " in stratum ", stratum,
                            " missing parameter(s): ", missing)
        return(FALSE)
      }
    }
  }
}


#' Disassembles a parameter table into a list of components.
#' 
#' Given a parameter table, this function returns a list with
#' the following named elements:
#' \describe{
#'   \item{site.meta}{Data frame of rows from the input table for
#'     site attributes.}
#'   \item{strata.meta}{Data frame of rows from the input table for
#'     stratum meta-data.}
#'   \item{species.values}{Data frame of species parameter values in
#'     wide format, ie. species as rows, parameters as columns.}
#'   \item{species.units}{Data frame of species parameter units in
#'     wide format.}
#' }
#' For the \code{species.values} data frame, all columns for numeric
#' parameters are converted from character to numeric.
#' 
#' @param tbl The input simulation parameter table (data frame).
#' 
#' @return A named list of table components.
#' 
#' @seealso \code{\link{ffm_assemble_table}} for the reverse process.
#' 
#' @importFrom dplyr %>%
#' 
#' @export
#' 
ffm_disassemble_table <- function(tbl) {
  # flag rows for species data
  spp <- !is.na( tbl$species )
  
  if (!("units" %in% colnames(tbl))) tbl$units <- NA_character_
  
  sp.values <- tbl %>%
    dplyr::filter(spp) %>%
    reshape2::dcast(stratum + species ~ param, value.var = "value")
  
  # Convert columns to numeric where applicable
  sp.values <- as.data.frame( lapply( sp.values, .optionally_numeric), 
                              stringsAsFactors=FALSE )
  
  sp.units <- tbl %>%
    dplyr::filter(spp) %>%
    reshape2::dcast(stratum + species ~ param, value.var = "units")
  
  list(site.meta = dplyr::filter(tbl, is.na(stratum)),
       strata.meta = dplyr::filter(tbl, !is.na(stratum), is.na(species)),
       species.values = sp.values,
       species.units = sp.units)
}

#' Assembles a parameter table from a list of elements.
#' 
#' This function takes a list of parameter table components as produced
#' by \code{\link{ffm_disassemble_table}} and combines them into a single
#' parameter table.
#' 
#' @param lst A named list of table components. 
#'   See \code{\link{ffm_disassemble_table}} for details of valid names.
#' 
#' @return A standard five column parameter table.
#' 
#' @importFrom dplyr %>% arrange left_join mutate_all
#' 
#' @export
#' 
ffm_assemble_table <- function(lst) {
  vals <- lst$species.values %>%
    reshape2::melt(id.vars = c("stratum", "species"), 
                   variable.name="param", value.name="value") %>%
    
    mutate_all(dplyr::funs(as.character))
  
  units <- lst$species.units %>%
    reshape2::melt(id.vars = c("stratum", "species"), 
                   variable.name="param", value.name="units") %>%
    
    mutate_all(dplyr::funs(as.character))
  
  spp <- left_join(vals, units, by=c("stratum", "species", "param"))
  
  rbind(lst$site.meta, lst$strata.meta, spp) %>%
    arrange(stratum, species)
}




#' Matches whole or partial species names to a table of default parameter values
#' for selected species.
#' 
#' For each name in the input character vector \code{names}, finds
#' the index or indices of matching rows in the table of default species
#' parameter values.
#'
#' @param default.species.params A data frame with 'name' as the first column
#'   (full species name) and subsequent columns giving default values for
#'   selected parameters. Parameter column names should match standard parameter
#'   names in \code{\link{ParamInfo}} (although case is ignored in comparisons).
#'
#' @param names The names, whole or partial, to search for.
#' 
#' @return a named \code{list} where names are the search strings and
#'   each value is the index or indices of matching species (NA for no match).
#'
#' @examples
#' \dontrun{
#' # A user-provided data frame of default parameter values 
#' # for species
#' my.defaults <- ...
#' 
#' # Look for Acacia dealbata and Cassinia aculeata.
#' # The search is case insensitive.
#' index <- ffm_find_species(my.defaults, c("Ac deal", "cass ac"))
#'
#' # Search for a name that (probably) does not exist in the default
#' # species parameters table. This results in \code{NA} being
#' # returned for the index:
#' ffm_find_species(my.defaults, "aspidistra")
#' }
#'
#' @export
#'
ffm_find_species <- function(default.species.parameters, names) {

  mkregex <- function(s) {
    parts <- stringr::str_split(s, "\\s")[[1]]
    parts <- stringr::str_c(parts, ".*", sep="")
    stringr::regex( stringr::str_c(parts, collapse=""), ignore_case = TRUE )
  }

  finder <- function(ptn) {
    matches <- stringr::str_detect(default.species.parameters$name, ptn)
    if (any(matches)) which(matches)
    else NA
  }

  ptns <- lapply(names, mkregex)
  res <- lapply(ptns, finder)
  names(res) <- names

  res
}


#' Checks if given species names have unique matches in table of parameter defaults
#'
#' For each name (whole or partial) checks whether the name matches one, and only
#' one, record in the provided default species parameters table.
#'
#' @param default.species.params A data frame with 'name' as the first column
#'   (full species name) and subsequent columns giving default values for
#'   selected parameters. Parameter column names should match standard parameter
#'   names in \code{\link{ParamInfo}}.
#'   
#' @param names The names, whole or partial, to search for.
#'   
#' @return A named boolean vector where names are search terms.
#' 
#' @seealso \code{\link{ffm_complete_params}}
#'
#' @export
#'
ffm_is_species_known <- function(names, default.species.params) {
  sapply(ffm_find_species(names),
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
#' params <- ffm_get_species_params(spp)
#' }
#'
#' @export
#'
ffm_get_species_params <- function(names) {
  indices <- ffm_find_species(names)

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
#' from the provided lookup table of defaults.
#' 
#' If a species has missing parameters but either it, or one or more of the
#' parameters, are not present in the lookup table of defaults an error is
#' thrown.
#' 
#' @param tbl The input simulation parameter table (data frame).
#' 
#' @param default.species.params A data frame with species as rows and
#'   parameters as columns giving default values to use for species parameters
#'   not specified in the input \code{params} table. Any subset of species
#'   parameters can be specified. Column names must match the parameter names in
#'   the \code{\link{ParamInfo}} table (case is ignored). An \code{NA} value
#'   indicates no default. This allows default values for a particular parameter
#'   to be provided for only a subset of species.
#' 
#' @return The completed parameter table
#' 
#' @importFrom dplyr %>% arrange left_join
#' 
#' @export
#' 
ffm_complete_params <- function(tbl, default.species.params) {
  if (!.is_param_table(tbl))
    stop("Input table must be a validly structured parameters data frame\n",
         "with columns: stratum, species, param, value and (optionally) units.")
  
  tbl <- .as_str_data_frame(tbl)
  
  # Retrieve default parameters as required for a given species
  # and return as additional records to add to the table
  do_species <- function(species.id) {
    recs <- dplyr::filter(tbl, species == species.id)
    stratum <- recs$stratum[1]
    species.name <- dplyr::filter(recs, param == "name")$value
    
    provided <- recs$param
    required <- setdiff( colnames(DefaultSpeciesParams)[-1], provided )
    
    if (length(required) == 0) {
      # Don't need to add any parameters, so return NULL
      NULL
    }
    else if (ffm_is_species_known(species.name)) {
      # Parameters required and species is known
      #
      dat <- ffm_get_species_params(species.name) %>% dplyr::select_(.dots = required)
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
      msg <- sprintf("Species %d (%s) is missing required parameters and is not in the defaults table",
                     species.id, species.name)
      stop(msg)
    }
  }
  
  ids <- .get_species_ids(tbl)
  new.recs <- do.call(rbind, lapply(ids, do_species))
  
  if (is.null(new.recs) || nrow(new.recs) == 0) {
    # no params added, return input
    tbl
    
  } else {
    # add a units column if one was present in the input table
    if (ncol(tbl) == 5) {
      units <- ParamInfo %>%
        dplyr::filter(section == "species") %>%
        dplyr::select(param, units)
      
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

# Finds parameter information for a given parameter label.
# 
# Searches for a row or rows in the \code{\link{ParamInfo}} table with a matching
# parameter label. The comparison ignores case and white-space. 
# If a non-NULL value is supplied for \code{section} the check is restricted
# to parameter labels in that section (site, stratum or species).
# 
# @param param The parameter label.
# 
# @param section If not \code{NULL}, one of site, stratum or species.
# 
# @param no.match.error If \code{TRUE} an error results when no
#   match is found, or more than one match is found and \code{single}
#   is \code{TRUE}; if \code{FALSE} (default) the function returns 
#   \code{NULL} if unsuccessful.
#   
# @param single If \code{TRUE} (default) require that label only
#   matches a single parameter; if \code{FALSE } allow multiple matches.
# 
# @return The index or indices (if \code{single} is \code{FALSE}) of the 
#   matching row(s) in \code{ParamInfo}; otherwise \code{NULL}
#   if unsuccessful and \code{no.match.error} is \code{FALSE}.
#
.match_param <- function(param, section, no.match.error = FALSE, single = TRUE) {
  if (is.null(section)) {
    labels <- ParamInfo$param
  } else {
    labels <- dplyr::filter(ParamInfo, section == section)$param
  }
  
  # Special treatment for the pesky 'w' parameter
  if (stringr::str_trim(tolower(param)) == "w") {
    ptn <- stringr::regex("^w$", ignore_case = TRUE)
  } else {
    ptn <- .make_ptn(param)
  }
  
  ii <- stringr::str_detect(labels, ptn)
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


# Updates the value of a parameter in a parameter table.
# 
# @param tbl The parameter table in which to set the value.
# @param stratum.id Stratum identifier.
# @param species.id Species identifier.
# @param param Label of the parameter to set.
# @param value Value of the parameter (will be converted to character).
# @param units (option) The units of measurement for the supplied value.
# 
# @return The updated parameter table.
# 
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

