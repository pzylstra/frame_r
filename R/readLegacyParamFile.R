#' Reads a legacy format parameters file
#' 
#' The legacy parameters file format was designed for use with a C++
#' implementation of the model. It does not conform to any more general format
#' (e.g. json or yaml) and, to be read correctly, relies on implicit assumptions
#' about the relationship between variables. For new parameter sets it is
#' preferable to compose them as tables in R or elsewhere and save them as
#' CSV-format files.
#'
#' @param path The full or relative path to the parameters file.
#'
#' @param raw (logical) Whether to return raw parsed data. This can be useful
#'   for diagnosing problems with reading a particular parameter file. It should
#'   be set to \code{FALSE} (default) to return a standard format table of
#'   parameters.
#'
#' @return If \code{raw} is \code{TRUE}, a data frame of parsed data
#'   in raw form (useful for testing).
#'   If \code{raw} is \code{FALSE} (the default), a data frame in a form 
#'   suitable for use with \code{\link{ffm_run}}.
#'   
#' @examples 
#' \dontrun{
#' # Import a legacy format parameter file and save it as a CSV format file
#' # in standard tabular parameter format.
#' params <- ffm_read_legacy_params("old_params.txt")
#' ffm_write_params(params, file = "beautiful_params.csv")
#' }
#'
#' @importFrom dplyr %>%
#'
#' @export
#'
ffm_read_legacy_params <- function(path, raw = FALSE) {
  txt <- .preprocess_text(path)
  dat <- .parse_text(txt)
  
  dat <- .parse_params(dat)
  
  if (raw) {
    # Return raw parsed input
    dat
  }
  else {
    dat <- dat %>% 
      # Convert parsed input to standard format
      .format_output() %>% 
      
      # Replace any missing ignition temperatures based
      # on function of silica free ash proportion
      .check_ignition_temps()
    
    # Return formatted parameters
    dat
  }
}


#' Reads a legacy format parameters file.
#' 
#' This function has been deprecated and will be removed in a future version of
#' the package. New code should use function
#' \code{\link{ffm_read_legacy_params}}.
#' 
#' @param path The full or relative path to the parameters file.
#'
#' @param raw (logical) Whether to return raw parsed data.
#'
#' @return If \code{raw} is \code{TRUE}, a data frame of parsed data
#'   in raw form (useful for testing).
#'   If \code{raw} is \code{FALSE} (the default), a data frame in a form 
#'   suitable for use with \code{\link{ffm_run}}.
#'   
#' @export
#'
readLegacyParamFile <- function(path, raw = FALSE) {
  .Deprecated("ffm_read_legacy_params")
  
  ffm_import_legacy_params(path, raw)
}


# Removes comments and blank lines.
#
.preprocess_text <- function(path) {
  if (!file.exists(path)) stop("File ", path, " not found")

  txt <- readLines(path)
  txt <- stringr::str_trim(txt)

  # remove blank lines
  ii <- stringr::str_length(txt) == 0
  txt <- txt[!ii]

  # remove whole comment lines
  ii <- stringr::str_detect(txt, "(^\\#)")
  txt <- txt[!ii]

  # remove trailing comments
  txt <- stringr::str_replace(txt, "(^[^#]+)(\\#.*$)", "\\1")
  txt <- stringr::str_trim(txt)

  txt
}


# Parses pre-processed parameter file text
# into a structured form.
#
#' @importFrom dplyr %>% group_by left_join mutate n ungroup
#
.parse_text <- function(txt) {
  # partition by strata
  strata <- .get_chunks(txt, "begin stratum", "end stratum", "stratum")

  # partition by species
  species <- .get_chunks(txt, "begin species", "end species", "species")

  out <- left_join(strata, species)
  out$txt <- txt

  out <- dplyr::select(out, line, txt, stratum, species)

  # flag stratum and species delimiter lines
  out <- out %>%
    group_by(stratum, species) %>%
    mutate(delim = c(1, rep(0, n() - 2), 1)) %>%
    ungroup()

  ii <- is.na(out$stratum)
  out$delim[ii] <- 0

  # flag meta-data lines
  out <- out %>%
    mutate(meta = as.integer(is.na(stratum)))

  # flag stratum meta-data lines
  out <- out %>%
    mutate(stratum.meta = as.integer(!is.na(stratum) & is.na(species) & !delim))

  #
  out
}


.parse_params <- function(dat,
                          in.param.col = "txt", in.delim.col = "delim",
                          out.name.col = "param", out.value.col = "value") {

  # dat is probably a tbl_df object, so to get a simple character
  # vector we have to do some undoing...
  x <- unname( unlist( dat[, in.param.col] ) )

  # overlapping parameters get special treatment
  overlap.ptn <- .make_ptn("overlapping")

  x <- lapply(x, .parse_param_expr, overlap.ptn)
  x <- as.data.frame( do.call(rbind, x), stringsAsFactors = FALSE )
  colnames(x) <- c(out.name.col, out.value.col)

  ii <- dat[, in.delim.col] == 1
  x[ii, c(out.name.col, out.value.col)] <- NA_character_

  cbind(dat, x)
}


.parse_param_expr <- function(expr, exclude) {
  # lhs is parameter text
  # rhs is value string: either a single value or two values delimited by comma
  parts <- stringr::str_trim( stringr::str_split(expr, "=")[[1]] )

  name <- parts[1]
  value <- NA_character_

  if (length(parts) > 1) {
    if (stringr::str_detect(name, exclude)) value <- parts[2]
    else value <- stringr::str_extract(parts[2], "^[^,]+")
  }

  c(name = name, value = value)
}

# Formats parsed data to be suitable for the API functions
#
#' @importFrom dplyr %>% arrange bind_rows distinct left_join mutate
#
.format_output <- function(dat) {
  dat <- dat %>%
    dplyr::filter(delim == 0) %>%
    dplyr::select(stratum, species, param, value) %>%
    .legacy_text_to_params

  # hack for absent dead fuel moisture parameter
  # in species - use surface dead fuel value
  deadMoist <- dat %>%
    dplyr::filter(is.na(stratum),
                  is.na(species),
                  stringr::str_detect(param, "^dead"))

  stopifnot(nrow(deadMoist) == 1)

  extras <- dat %>%
    dplyr::filter(!is.na(species)) %>%
    dplyr::select(stratum, species) %>%
    distinct %>%
    mutate(param = "deadLeafMoisture", value = deadMoist$value)

  dat %>%
    bind_rows(extras) %>%
    arrange(stratum, species) %>%
    .minus99_to_NA %>%
    left_join(dplyr::select(ParamInfo, param, units), by="param")
}


# Replace any missing ignition temperatures with values
# calculated from proportion of silica free ash.
#
#' @importFrom dplyr %>% mutate
#
.check_ignition_temps <- function(tbl) {
  
  x <- tbl %>%
    dplyr::filter(param %in% c("ignitionTemp", "propSilicaFreeAsh")) %>%
    reshape2::dcast(stratum + species ~ param) %>%
    dplyr::filter(is.na(ignitionTemp))
  
  if ( nrow(x) > 0 ) {
    x <- x %>%
      mutate(ignitionTemp = ignitionTempPSFA( as.numeric(propSilicaFreeAsh) ) ) %>%
      reshape2::melt(id.vars = c("stratum", "species"), 
                     variable.name = "param", value.name = "new.value")
    
    tbl <- tbl %>% left_join(x) %>%
      mutate(value = ifelse( !is.na(new.value), new.value, value) ) %>%
      dplyr::select(-new.value)
  }
  
  tbl
}


# Takes a data frame with columns 'param' and 'value' and
# matches the legacy parameter text in 'param' to names in
# the lookup table.
.legacy_text_to_params <- function(param.dat) {
  ii <- match( .munge_string(param.dat$param),
               .munge_string(LegacyParamLookup$text) )

  param.dat$param <- LegacyParamLookup$name[ii]

  dplyr::filter(param.dat, !is.na(param))
}


.minus99_to_NA <- function(dat) {
  ii <- .munge_string(dat$value) == "-99"
  dat$value[ii] <- NA

  dat
}


# Finds start and end positions of text chunks defined by
# the given start and end labels.
#
.get_chunks <- function(txt, start.txt, end.txt,
                       chunk.label = "chunk") {

  start.pos <- which( stringr::str_detect(txt, .make_ptn(start.txt)) )
  end.pos <- which( stringr::str_detect(txt, .make_ptn(end.txt)) )

  if (length(start.pos) != length(end.pos))
    stop("Unmatched labels for chunks defined by ",
         start.txt, " and ", end.txt)

  N <- length(start.pos)
  out <- data.frame(line = 1:length(txt), chunk = NA_integer_)

  if (N > 0)
    for (i in 1:N) out$chunk[ start.pos[i]:end.pos[i] ] <- i

  colnames(out)[2] <- chunk.label

  out
}

