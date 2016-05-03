#' Creates a Site object from a table of parameters.
#' 
#' @param params A data.frame or matrix of character data 
#'   with either four or five columns:
#'   \enumerate{
#'     \item stratum identifier (\code{NA} for rows corresponding to
#'        site meta-data such as weather variables)
#'     \item species identifier (\code{NA} for rows corresponding to
#'        stratum meta-data such as the level name)
#'     \item param standard parameter label
#'     \item value parameter value
#'     \item units (optional) parameter units or \code{NA} for 
#'        dimensionless parameters
#'   }
#'   
#' @return a \code{ScalaInterpreterReference} holding the
#'   created Site object.
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
#' @param x Site data as a ScalaInterpreterReference object
#'   created with \code{\link{ffm_site}} or a matrix or data.frame
#'   valid for input to \code{ffm_site}.
#'
#' @export
#'
ffm_run <- function(x, as.ref) UseMethod("ffm_run")


#' Run the model with the given site object.
#'
#' @param site The site as a ScalaInterpreterReference object created
#'   with \code{\link{ffm_site}}.
#'   
#' @return ???
#'
#' @export
#'
ffm_run.ScalaInterpreterReference <- function(site, as.ref=FALSE) {
  stopifnot(is(site, "ScalaInterpreterReference"))

  res <- E$s$do('ffm.runner.Runner')$run(site)
  
  if (as.ref) res
  else ffm_result(res)
}


#' Run the model with a site based on parameters in the
#' given data frame.
#'
#' @param site.params A data frame of site parameters valid for
#'   input to \code{\link{ffm_site}}
#'
#' @export
#'
ffm_run.data.frame <- function(site.params, as.ref=FALSE) {
  ffm_run( ffm_site(site.params), as.ref )
}


#' Run the model with a site based on parameters in the
#' given character matrix.
#'
#' @param site.params A character matrix of site parameters valid for
#'   input to \code{\link{ffm_site}}
#'
#' @export
#'
ffm_run.matrix <- function(site.params, as.ref=FALSE) {
  ffm_run( ffm_site(site.params), as.ref )
}


#' Extracts data from a flammability model result object
#' 
#' @export
#' 
ffm_result <- function(res) {
  stopifnot(is(res, "ScalaInterpreterReference"))
  
  # ObjF <- E$s$do('ffm.io.r.ObjectFactory')
  
  st.res <- .get_stratum_results(res)
  run.res <- .get_run_results(res)
  
  list(stratum.results = st.res, run.results = run.res)
}


# Converts a data.frame or numeric matrix to a matrix
# of strings.
.as_str_matrix <- function(dat) apply(dat, MARGIN = 2, FUN = as.character)


.get_stratum_results <- function(res) {
  stopifnot(is(res, "ScalaInterpreterReference"))
  
  levels <- res$'stratumResults.map( sr => sr.level.toString ).toArray'()
  flames <- res$'stratumResults.map(sr => Array(sr.flameLength, sr.flameAngle, sr.flameHeight)).toArray'()
  
  data.frame(level = levels,
             flameLength = flames[,1],
             flameAngle = flames[,2] * 180 / pi,
             flameHeight = flames[,3])
}


.get_run_results <- function(res) {
  .get_surface_results(res)
}

.get_surface_results <- function(res) {
  dat <- E$s$do('ffm.io.r.ResultFormatter')$'surfaceResults'(res)
  
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  colnames(dat) <- dat[1,]
  dat <- dat[-1, ]
  
  dat[,1] <- as.integer(dat[,1])
  for (i in 2:ncol(dat)) dat[,i] <- as.numeric(dat[,i])
  
  dat$run <- dat$run + 1
  dat$flameAngle <- dat$flameAngle * 180 / pi
  
  rownames(dat) <- NULL
  dat
}
