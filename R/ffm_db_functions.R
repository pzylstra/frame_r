#' Opens a connection to a SQLite database of model results.
#' 
#' This function attempts to connect to a SQLite database
#' file. If successful it returns a connection which can be
#' used with functions in the RSQLite or dplyr packages to query data.
#' 
#' @param db.path The path and filename for the database.
#' 
#' @return A connection object.
#' 
#' @examples 
#' \dontrun{
#' # Querying the database directly with RSQLite functions
#' 
#' library(RSQLite)
#' 
#' # Open a connection
#' con <- ffm_db_connect("c:/michael/somewhere/my_results.db")
#' 
#' # Query the FlameSummaries table for maximum flame height
#' # in each run
#' sql <- "SELECT repId, MAX(flameHeight) AS height FROM FlameSummaries GROUP BY repId"
#' res <- dbGetQuery(con, sql)
#' 
#' # Close the connection when you are finished with it
#' dbDisconnect(con)
#' 
#' # TODO - An example of querying the database with
#' # the dplyr package
#' }
#' 
#' @export
#' 
ffm_db_connect <- function(db.path) {
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db.path)
  if (RSQLite::dbIsValid(con)) con
  else stop("Unable to connect to the database ", db.path)
}


#' Attempts to load the entire contents of a database into
#' the R workspace.
#' 
#' Given a path and filename for a database of model results,
#' this function connects to the database and loads the full
#' contents of all tables into a \code{list} of data frames.
#' An error is returned if the database cannot be opened.
#' 
#' @param db.path The path and filename for the database.
#' 
#' @return A named list of data frames corresponding to
#'   database tables.
#'   
#' @examples 
#' \dontrun{
#' # Load the contents of a database and display the start
#' # of the FlameSummaries table
#' tbls <- ffm_db_load("c:/michael/somewhere/my_results.db")
#' head( tbls$FlameSummaries)
#' }
#'
#' @export
ffm_db_load <- function(db.path) {
  con <- ffm_db_connect(db.path)
  tbl.names <- RSQLite::dbListTables(con)
    
  res <- lapply(tbl.names, function(tbl.name) {
    RSQLite::dbReadTable(con, tbl.name)
  })
  
  RSQLite::dbDisconnect(con)
    
  names(res) <- tbl.names
  res
}


#' Summarize the contents of a model results database.
#' 
#' At the moment this function simply prints the number of
#' simulation runs found in the database.
#' 
#' @param db.path Path and filename for the database.
#' 
#' @export
#' 
ffm_db_summary <- function(db.path) {
  con <- ffm_db_connect(db.path)
  sql <- "SELECT DISTINCT(repId) AS repId FROM Runs"
  dat <- RSQLite::dbGetQuery(con, sql)
  nruns <- nrow(dat)
  
  RSQLite::dbDisconnect(con)
  
  msg <- if (nruns == 1) "1 run" else paste(nruns, "runs")
  cat("Database contains results for", msg, "\n")
}

