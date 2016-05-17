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


#' Runs the model with the given site.
#'
#' @param x Site data (parameter set) as one of the following:
#'   \itemize{
#'     \item a ScalaInterpreterReference object created with \code{\link{ffm_site}}
#'     \item a a matrix or data.frame valid for input to \code{ffm_site}
#'     \item a parameter generating function.
#'   }     
#'   
#' @param db Output database
#'   
#' @param additional arguments (see specific ffm_run functions for details)
#' 
#' @seealso \code{\link{ffm_run.ScalaInterpreterReference}}
#' @seealso \code{\link{ffm_run.data.frame}}
#' @seealso \code{\link{ffm_run.matrix}}
#' @seealso \code{\link{ffm_run.function}}
#'
#'    
#' @export
#'
ffm_run <- function(x, db, ...) UseMethod("ffm_run")


#' Run the model with the given site object.
#'
#' @param site The site as a Scala reference object created
#'   with \code{\link{ffm_site}}.
#'   
#' @param db Output database as a Scala reference object.
#'   
#' @return \code{TRUE} if the run completed and results were written to
#'   the output database successfully; \code{FALSE} otherwise.
#'
#' @seealso \code{\link{ffm_create_database}}
#'
#' @export
#'
ffm_run.ScalaInterpreterReference <- function(site, db, ...) {
  stopifnot(is(site, "ScalaInterpreterReference"))
  
  db <- .check_db(db)

  res <- E$s$do('ffm.runner.Runner')$run(site)
  
  db$insertResult(res)
}


#' Run the model with a site based on parameters in the
#' given data frame.
#'
#' @param site.params A data frame of site parameters valid for
#'   input to \code{\link{ffm_site}}
#'
#' @export
#'
ffm_run.data.frame <- function(site.params, db, ...) {
  ffm_run( ffm_site(site.params), db, ... )
}


#' Run the model with a site based on parameters in the
#' given character matrix.
#'
#' @param site.params A character matrix of site parameters valid for
#'   input to \code{\link{ffm_site}}
#'
#' @export
#'
ffm_run.matrix <- function(site.params, db, ...) {
  ffm_run( ffm_site(site.params), db, ... )
}


#' Run the model with parameters provided by a generating function.
#'
#' @param gen.fn A function which generates a valid parameter matrix.
#'
#' @export
#'
ffm_run.function <- function(gen.fn, db, ...) {
  params <- gen.fn(...)
  ffm_run.ScalaInterpreterReference( ffm_site(params), db, ... )
}

#' Creates a database manager object and file.
#' 
#' This function creates a Scala database manager (ffm.io.r.Database object),
#' together with an associated SQLite database file. A reference is returned
#' which can then be passed to \code{\link{ffm_run}} functions and used to
#' write model results.
#'
#' @param path The path for the database file
#' 
#' @param delete.existing If \code{TRUE} (default), any existing database file
#'   matching \code{path} is deleted if possible; if \code{FALSE}, the database
#'   will only be created if the file does not exist.
#'   
#' @param use.transactions If \code{TRUE} (default), all insertions into the database
#'   will be done within transactions (slow but safe); if \code{FALSE},
#'   insertions will be done directly (fast but less safe).
#'   
#' @return A reference object (ScalaInterpreterReference) which can be 
#'   passed to \code{\link{ffm_run}} functions.
#' 
#' @export
#' 
ffm_create_database <- function(path, delete.existing = TRUE, use.transactions = TRUE) {
  optionDB <- E$s$do('ffm.io.r.Database')$create(path, delete.existing, use.transactions)
  if (optionDB$isDefined()) optionDB$get()
  else stop("Database ", path, " could not be created")
}


#' Locate specified site meta-parameters in a parameter table.
#' 
#' @param params Parameter table as a matrix or data.frame of character values.
#' @param names The name(s) of parameters to find.
#' 
#' @return Row numbers corresponding to the specified parameters.
#' 
#' @export
#' 
ffm_find_site_param <- function(params, names) {
  stopifnot(is.data.frame(params) | is.matrix(params))
  
  params <- as.data.frame(params)
  ii <- stringr::str_detect(params$param, names) & 
    is.na(params$stratum) & is.na(params$species)
  
  which(ii)
}


.check_db <- function(db) {
  if(!is(db, "ScalaInterpreterReference"))
    stop("Output database must be supplied as a Scala reference object")
  
  if (!(db$isOpen(TRUE)))
    stop("Database is not open for writing results")
  
  db
}
