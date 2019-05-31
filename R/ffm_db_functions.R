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
  con <- .get_sqlite_connection(db.path)
  
  tbl.names <- RSQLite::dbListTables(con)
  
  res <- lapply(tbl.names, function(tbl.name) {
    RSQLite::dbReadTable(con, tbl.name)
  })
  
  DBI::dbDisconnect(con)
    
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
  con <- .get_sqlite_connection(db.path)
  
  sql <- "SELECT DISTINCT(repId) AS repId FROM Runs"
  dat <- RSQLite::dbGetQuery(con, sql)
  nruns <- nrow(dat)
  
  DBI::dbDisconnect(con)
  
  msg <- if (nruns == 1) "1 run" else paste(nruns, "runs")
  cat("Database contains results for", msg, "\n")
}


.get_sqlite_connection <- function(db.path) {
  if (!file.exists(db.path)) stop("Cannot find file: ", db.path)
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db.path)
  if (!DBI::dbIsValid(con)) stop("Unable to connect to database: ", db.path)
  
  con
}

