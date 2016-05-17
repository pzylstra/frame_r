#' Reads a legacy format parameters file.
#'
#' @param path The full or relative path to the parameters file.
#'
#' @param raw (logical) Whether to return raw parsed data.
#'
#' @return If \code{raw} is \code{TRUE}, a data frame of parsed data
#'   in raw form (useful for testing); 
#'   if \code{raw} is \code{FALSE} (the default), a data frame in a form 
#'   suitable to create Scala objects (e.g. with \code{\link{ffm_site}}).
#'
#' @export
#'
readLegacyParamFile <- function(path, raw = FALSE) {
  txt <- .preprocess_text(path)
  dat <- .parse_text(txt)

  dat <- .parse_params(dat)

  if (raw) dat
  else .format_output(dat)
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
.parse_text <- function(txt) {
  # partition by strata
  strata <- .get_chunks(txt, "begin stratum", "end stratum", "stratum")

  # partition by species
  species <- .get_chunks(txt, "begin species", "end species", "species")

  out <- left_join(strata, species)
  out$txt <- txt

  out <- select(out, line, txt, stratum, species)

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

# Formats parsed data to be suitable for Scala api functions
# such as ffm_site.
#
.format_output <- function(dat) {
  dat <- dat %>%
    filter(delim == 0) %>%
    select(stratum, species, param, value) %>%
    .legacy_text_to_params

  # hack for absent dead fuel moisture parameter
  # in species - use surface dead fuel value
  deadMoist <- dat %>%
    filter(is.na(stratum),
           is.na(species),
           stringr::str_detect(param, "^dead"))

  stopifnot(nrow(deadMoist) == 1)

  extras <- dat %>%
    filter(!is.na(species)) %>%
    select(stratum, species) %>%
    distinct %>%
    mutate(param = "deadLeafMoisture", value = deadMoist$value)

  dat %>%
    rbind(extras) %>%
    arrange(stratum, species) %>%
    .minus99_to_NA %>%
    left_join(select(ParamInfo, param, units), by="param")
}


# Takes a data frame with columns 'param' and 'value' and
# matches the legacy parameter text in 'param' to names in
# the lookup table.
.legacy_text_to_params <- function(param.dat) {
  ii <- match( .munge_string(param.dat$param),
               .munge_string(LegacyParamLookup$text) )

  param.dat$param <- LegacyParamLookup$name[ii]

  filter(param.dat, !is.na(param))
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

