# Non-explored helper functions for use by other
# parts of the package

# Converts a data.frame or numeric matrix to a data.frame 
# of character values.
.as_str_data_frame <- function(dat) {
  dat <- as.data.frame(dat)
  
  # use dplyr to avoid the data.frame being simplified to a matrix
  dplyr::mutate_all(dat, dplyr::funs(as.character))
}


# Converts a data.frame or numeric matrix to a matrix
# of character values.
.as_str_matrix <- function(dat) apply(dat, MARGIN = 2, FUN = as.character)


# Runs the most basic checks to see if an object can be 
# used as a parameter table
.is_param_table <- function(tbl) {
  # must be a data frame
  ok <- is(tbl, "data.frame")
  
  # must have some rows
  ok <- ok && nrow(tbl) > 0
  
  # must have 4 or 5 columns
  ok <- ok && ncol(tbl) %in% 4:5
  
  # must have mandatory columns correctly named
  ok <- ok && c("stratum", "species", "param", "value") %in% colnames(tbl)
  
  # if there are 5 columns, one should be "units"
  ok <- ok && "units" %in% colnames(tbl)
  
  ok
}


# Gets species IDs in the given parameter table
#
.get_species_ids <- function(tbl) unique( na.omit(tbl$species) )


# Creates a regular expression pattern for the given string.
#
# Spaces in the input string are replaced by `\\s*` (treating
# spaces as optional) and the pattern is set to be case
# insensitive.
#
.make_ptn <- function(ptn) {
  ptn <- stringr::str_replace(ptn, "\\s+", "\\\\s*")
  stringr::regex(ptn, ignore_case = TRUE)
}


# Converts strings to lower case and removes spaces.
#
.munge_string <- function(s) {
  s %>% tolower %>% stringr::str_replace_all("\\s+", "")
}


# Converts a vector of values to numeric if possible.
#
# @param x Input values

# @param partial If \code{TRUE} does conversion even if only
#   possible for some elements, returning NA for others.
#
# @return Numeric values if conversion was possible, otherwise
#   the input values unaltered.
#
.optionally_numeric <- function(x, partial = TRUE) {
  z <- suppressWarnings( as.numeric(x) )
  nas <- is.na(z)
  
  out <- 
    if(all(nas)) x 
  else if (any(nas)) {
    if (partial) z
    else x
  } 
  else z
}
