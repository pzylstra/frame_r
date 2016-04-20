#' Create a Site object from a table of parameters.
#'
#' @export
#'
ffm_site <- function(params) {
  stopifnot(is.data.frame(params) || is.matrix(params))

  params <- .as_str_matrix(params)
  E$s$do('ffm.io.r.ObjectFactory')$createSite(params)
}

#' Run the model with the given site.
#'
#' @param The site as a ScalaInterpreterReference object.
#'
#' @export
#'
ffm_run <- function(site) {
  stopifnot(is(site, "ScalaInterpreterReference"))

  E$s$do('ffm.runner.Runner')$run(site)
}

# Converts a data.frame or numeric matrix to a matrix
# of strings.
.as_str_matrix <- function(dat) apply(dat, MARGIN = 2, FUN = as.character)

