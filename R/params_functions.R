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
#' @seealso \code{\link{ffm_complete_params}}
#'
#' @importFrom dplyr %>% group_by mutate ungroup
#'  
#' @export
#' 
ffm_check_params <- function(params, quiet = FALSE) {
  # optional species parameter
  Silica <- tolower("propSilicaFreeAsh")
  
  # Just in case
  params <- ungroup(params)
  
  required.cols <- c("stratum", "species", "param", "value")
  nrequired <- length(required.cols)
  
  optional.cols <- "units"
  
  given.cols <- colnames(params)
  missing.cols <- setdiff(required.cols, given.cols)
  
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
  XParamInfo <- frame::ParamInfo %>%
    mutate(loparam = tolower(param))
  
  params <- params %>%
    mutate(loparam = tolower(param))
  
  # Function to format multiple parameter names for messages
  fconcat <- function(x) paste(x, collapse = ", ")
  
  
  # Check site-level parameters
  given <- params %>%
    dplyr::filter(is.na(stratum))
  
  expected <- XParamInfo %>%
    dplyr::filter(section == "site")
  
  missing <- setdiff(expected$loparam, given$loparam)
  if (length(missing) > 0) {
    if (!quiet) warning("Missing site-level parameter(s): ", 
                        fconcat(missing))
    return(FALSE)
  }
  
  extra <- setdiff(given$loparam, expected$loparam)
  if (length(extra) > 0) {
    if (!quiet) warning("Unrecognized site-level parameter(s): ", 
                        fconcat(extra))
    return(FALSE)
  }
  
  # Check strata- and species-level parameters
  strata <- sort(unique(params$stratum), na.last = NA)
  
  strata.expected <- {
    x <- dplyr::filter(XParamInfo, section == "stratum")
    x$loparam
  }
  
  species.expected <- {
    x <- dplyr::filter(XParamInfo, section == "species")
    x$loparam
  }
  
  for (ist in strata) {
    sdat <- dplyr::filter(params, stratum == ist)

    # stratum level
    given <- dplyr::filter(sdat, is.na(species))    
    missing <- setdiff(strata.expected, given$loparam)
    if (length(missing) > 0) {
      if (!quiet) warning("Stratum ", ist, " missing parameter(s): ", 
                          fconcat(missing))
      return(FALSE)
    }
    
    # species level
    given <- dplyr::filter(sdat, !is.na(species))
    spp <- sort(unique(given$species))
    
    for (isp in spp) {
      givensp <- dplyr::filter(given, species == isp)
      missing <- setdiff(species.expected, c(given$loparam, Silica))
      if (length(missing) > 0) {
        if (!quiet) warning("Species ", isp, " in stratum ", ist,
                            " missing parameter(s): ", 
                            fconcat(missing))
        return(FALSE)
      }
    }
  }
  
  TRUE
}


#' Writes a table of parameters to a CSV format file.
#'
#' This is simply a wrapper that calls the standard
#' \code{write.csv} function with appropriate default arguments.
#'
#' @param tbl The input simulation parameter table (data frame).
#'
#' @param path The path and file name to write to. The recommended file
#'   extension is \code{'.csv'} but this is not enforced by the function.
#'
#' @param overwrite Set to \code{FALSE} (default) to issue an error message if
#'   the output file exists, or \code{TRUE} to overwrite any existing file.
#'   
#' @param must.be.complete Set to \code{FALSE} (default) to enable writing a
#'   partial table of parameters, or to \code{TRUE} to abort writing and issue an
#'   error message if the parameter table is not complete.
#'
#' @export
#' 
ffm_write_params <- function(tbl, path, overwrite = FALSE, must.be.complete = FALSE) {
  if (!overwrite && file.exists(path))
    stop("Output file ", path, " exists. Set overwrite=TRUE if you want to replace it.")
  
  if (must.be.complete) {
    ok <- ffm_check_params(tbl, quiet = TRUE)
    if (!ok) stop("Parameters table is not complete. Set must.be.complete = FALSE to ignore this check.")
  }
  
  write.csv(tbl, file = path, row.names = FALSE)
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
  
  spp <- left_join(vals, units, by=c("stratum", "species", "param")) %>%
    mutate(species = as.integer(species))
  
  rbind(lst$site.meta, lst$strata.meta, spp) %>%
    arrange(stratum, species) %>%
    
    mutate_all(dplyr::funs(as.character))
}


#' Completes a table by adding default parameters for each species as required.
#' 
#' For each species in the input parameter table, this function checks which
#' parameters have been provided and retrieves defaults for missing parameters
#' from the provided lookup table of defaults.
#' 
#' If a species has missing parameters but either it, or one or more of the
#' parameters, are not present in the lookup table of defaults an error is
#' thrown. An exception is made for the parameter 'propSilicaFreeAsh' which
#' is allowed to be missing.
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
#' @seealso \code{\link{ffm_check_params}}
#' 
#' @importFrom dplyr %>% arrange left_join
#' 
#' @export
#' 

ffm_complete_params <- function(tbl, default.species.params) {
  if (!frame:::.is_param_table(tbl))
    stop("Input table must be a validly structured parameters data frame\n",
         "with columns: stratum, species, param, value and (optionally) units.")
  
  if (!is.data.frame(default.species.params) || 
      ncol(default.species.params) < 2 ||
      colnames(default.species.params)[1] != "name") {
    
    stop("default.species.params should be a data frame with 'name' as \n",
         "first column and valid parameter names as further columns")
  }
  
  if ( any(table(tolower(default.species.params$name)) > 1) ) {
#    stop("One or more duplicate species names in default.species.params")
    stop(paste0("These species have duplicates: ", names(which(table(tolower(traits$name)) > 1)),"\n"))
  }
  
  # Convert default to long format
  default.species.params <- default.species.params %>%
    tidyr::gather(param, value, -name)
  
  tbl <- frame:::.as_str_data_frame(tbl)
  
  Silica <- tolower("propSilicaFreeAsh")
  
  RequiredSpeciesParams <- {
    x <- frame::ParamInfo %>% dplyr::filter(section == "species")
    tolower(x$param)
  }
  
  
  # Retrieve default parameters as required for a given species
  # and return as additional records to add to the table
  do_species <- function(species.id) {
    all.recs <- dplyr::filter(tbl, species == species.id)
    species.name <- dplyr::filter(all.recs, param == "name")$value[1]
    
    strata <- sort( unique(all.recs$stratum) )
    if (anyNA(strata)) 
      stop("Parameter record for species ", species.name, " is missing stratum ID")
    
    # Look for missing parameters within each stratum
    extra.recs <- NULL
    for (ist in strata) {
      stratum.recs <- dplyr::filter(all.recs, stratum == ist)
      provided <- tolower(stratum.recs$param)
      required <- setdiff(RequiredSpeciesParams, provided)
      
      if (length(required) == 0) {
        # Don't need to add any parameters, so return NULL
        NULL
      }
      else {
        spdefaults <- dplyr::filter(default.species.params, 
                                    name == species.name,
                                    !is.na(value))
        
        if (nrow(spdefaults) == 0) {
          msg <- glue::glue(
            "Species {species.name} in stratum {ist} is missing required \n",
            "parameters but is not in the table of defaults")
          
          stop(msg)
        }
        
        # Are all required parameters in the defaults?
        # (allow propSilicaFreeAsh to be missing)
        missing <- setdiff(required, tolower(spdefaults$param))
        if (length(setdiff(missing, Silica)) > 0) {
          msg <- glue::glue(
            "Species {species.name} in stratum {ist} is missing required \n",
            "parameters that are not provided by the table of defaults")
          
          stop(msg)
        }
        
        # Get required additional parameter table records
        spextras <- spdefaults %>%
          dplyr::filter(tolower(param) %in% required) %>%
          dplyr::select(param, value)
        
        spextras$stratum <- ist
        spextras$species <- species.id
        
        spextras <- arrange(spextras, stratum, species, param, value)
        
        extra.recs <- rbind(extra.recs, spextras)
      } 
    }
    extra.recs
  }
  
  ids <- frame:::.get_species_ids(tbl)
  new.recs <- do.call(rbind, lapply(ids, do_species))
  
  if (is.null(new.recs) || nrow(new.recs) == 0) {
    # no params added, return input
    tbl
    
  } else {
    # add a units column if one was present in the input table
    if (ncol(tbl) == 5) {
      units <- frame::ParamInfo %>%
        dplyr::filter(section == "species") %>%
        dplyr::select(param, units)
      
      new.recs <- left_join(new.recs, units, by = "param")
    }
    
    tbl <- rbind( tbl, new.recs )
    arrange(tbl, stratum, as.integer(species))
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

  if (!is.null(i)) frame::ParamInfo[i, ]
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
    labels <- frame::ParamInfo$param
  } else {
    labels <- dplyr::filter(frame::ParamInfo, section == section)$param
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
  std.param <- frame::ParamInfo[i, "param"]

  in.scope <- .match_ids(tbl, stratum.id, species.id)
    
  rows <- which( tbl$param == std.param & in.scope ) 
  if (length(rows) == 0) 
    stop("Input table does not contain parameter ", 
         std.param,
         " (", section, ") ",
         "for stratum=", stratum.id,
         ", species=", species.id,
         ". This species is not present.")
  
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

