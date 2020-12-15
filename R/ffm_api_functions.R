#' Run the model with a given table of parameters.
#' 
#' This function requests the Scala model application to run a simulation with
#' the given input parameters and write results to a SQLite database file. For a
#' related set of simulations (e.g. when exploring the effects of variation over
#' a subset of parameters), the function can be called repeatedly (once for each
#' realization of the parameters) with the same database file specified, and each
#' set of results will be appended to the database with an incremented replicate 
#' ID number.
#'
#' To run the external Scala model application, the function first writes the
#' parameters to a temporary CSV-format file. It then issues a system call to
#' use the Java installation on your machine to run the simulation with the
#' parameters file and the output database specified by \code{db.path}. For this
#' to work, you must have a properly installed version of Java (version 1.8.x or
#' higher) that can be run from any directory by issuing the command
#' \code{'java'}. Consult the Oracle Java installation instructions for more
#' details.
#' 
#' 
#' @param params A data frame of parameters
#'   
#' @param db.path The path and filename of a new or existing SQLite database
#'   file to which simulation results will be written.
#'   
#' @param default.species.params An optional data frame with species as rows and
#'   parameters as columns giving default values to use for species parameters
#'   not specified in the input \code{params} table. Any subset of species
#'   parameters can be specified. Column names must match the parameter names in
#'   the \code{\link{ParamInfo}} table (case is ignored). An \code{NA} value
#'   indicates no default. This allows default values for a particular parameter
#'   to be provided for only a subset of species.
#'   
#' @param db.recreate If \code{TRUE} and the file specified by \code{db.path}
#'   exists, the database will be recreated and any prior results will be lost.
#'   If \code{FALSE} (default) and the file exists, the simulation results
#'   will be appended to the database with an incremented replicate ID value
#'   (field \code{repId} in all tables).
#' 
#' @return \code{TRUE} if the run completed and results were written to
#'   the output database successfully; \code{FALSE} otherwise.
#'
#' @export
#'
ffm_run <- function(params, db.path, 
                    default.species.params = NULL, 
                    db.recreate = FALSE) {
  
  db.path <- normalizePath(db.path, mustWork = FALSE)
  
  if (!.check_setting("java", TRUE)) {
    message("Checking for Java...\n")
    
    ok <- ffm_check_java()
    .set_setting("java", ok)
    
    if (!ok) 
      stop("Bummer: Try again when Java 1.8 or higher is properly installed")
  }  
  
  if (!ffm_check_params(params, quiet = TRUE)) {
    # Attempt to complete parameters if default are provided.
    # If this fails, ffm_complete_params will throw an error.
    if (!is.null(default.species.params)) {
      params <- ffm_complete_params(params, default.species.params)
    }
    else {
      # No defaults - give up
      warning("Parameters incomplete and no defaults provided")
      FALSE
    }
  }
  
  # If we got this far we should have a complete set of parameters
  param.path <- tempfile(pattern = "ffm_", fileext = ".csv")
  ffm_write_params(params, path = param.path)
  
  # Compose sytem call and run the simulation
  cmd <- ffm_run_command(param.path = param.path, 
                         db.path = db.path,
                         db.recreate = db.recreate)
  
  res <- system(cmd, intern = TRUE)
    
  # Check for 'success' in Scala output
  success <- any( stringr::str_detect(tolower(res), "success") )
  if (!success) warning(res)

  # Return status as TRUE / FALSE  
  success
}


#' Composes a command line call to run a simulation
#' 
#' This function is used by \code{\link{ffm_run}} when calling the Scala model
#' application via the base R \code{\link[base]{system}} function. The command
#' line string passed to \code{system} invokes the local Java runtime and
#' includes arguments for the necessary Java libraries (jar files included with
#' the package), the CSV-format file of input parameters, the path to the output
#' database, and (optionally) a flag indicating that the database should be
#' recreated if the file exists.
#' 
#' @param param.path Path and name of a CSV-format file of simulation 
#'   parameters.
#'   
#' @param db.path Path and name of the output SQLite database file.
#' 
#' @param db.recreate A logical value indicating, if the file specified 
#'   by \code{db.path} exists, whether to recreating the database (\code{TRUE})
#'   or append to it (\code{FALSE}). Default is to append.
#'   
#' @param runtime.class The name of the runtime class in the Scala model
#'   application. Unless you have modified the Scala application this should
#'   be left as the default value: \code{ffm.runner.CSVRunner}.
#' 
#' @return The command line instruction as a one-element character vector.
#' 
#' @export
#' 
ffm_run_command <- function(param.path, db.path, 
                            db.recreate = FALSE,
                            runtime.class = "ffm.runner.CSVRunner") {
  
  # Guard against paths or filenames with embedded spaces
  safepath <- function(path) {
    ifelse (stringr::str_detect(path, "\\s"), paste0('"', path, '"'), path)
  }

  pkgdir <- safepath( system.file(package = "frame") ) 

  jarfiles <- safepath(
    dir(file.path(pkgdir, "java"), 
        pattern = "\\.jar$", full.names = TRUE)
  )
  
  cp <- paste(jarfiles, collapse = ";")
  
  # convert logical db.recreate to an optional flag
  db.recreate <- if (db.recreate) "-x" else ""
  
  glue::glue('java -cp {cp} {runtime.class} -p "{param.path}"',
             ' -d "{db.path}" {db.recreate}')
}


#' Check for a locally installed Java runtime system
#' 
#' Checks that Java version 1.8 or higher (aka Java 8) is installed 
#' on the local system and available from the command line. This is
#' required to run simulations with \code{\link{ffm_run}}.
#' 
#' @return \code{TRUE} if a compatible version of Java is found; 
#'   \code{FALSE} otherwise.
#'
#' @export
#' 
ffm_check_java <- function(quiet = FALSE) {
  x <- system2("java", args = "-version", stdout = TRUE, stderr=TRUE)
  
  found <- any( grepl("version", x, ignore.case = TRUE) )
  
  ver <- 
    if (!found) NA
  else as.numeric(stringr::str_extract(x[1], "\\d\\.\\d+"))
  
  if (is.na(ver)) {
    if (!quiet)
      warning("Java was not found on the local system path\n",
              "Java version 1.8 or higher is required")
    FALSE
    
  } else if (ver < 1.8) {
    if (!quiet)
      warning("Java version ", ver, " found\n",
              "Version 1.8 or higher is required")
    FALSE
    
  } else {
    if (!quiet)
      message("Compatible Java version found (", ver, ")")
    
    TRUE
  }
}


.check_setting <- function(setting, desired.value) {
  if (!exists("._ffm_settings", where = 1)) {
    assign("._ffm_settings", list(), pos = 1)
    FALSE
  }
  
  S <- get("._ffm_settings", pos = 1)
  
  if (!(setting %in% names(S))) {
    FALSE
  } else {
    S[[setting]] == desired.value
  }
}

.set_setting <- function(setting, value) {
  if (!exists("._ffm_settings", where = 1)) {
    assign("._ffm_settings", list(), pos = 1)
  }
  
  S <- get("._ffm_settings", pos = 1)
  S[[setting]] <- value
  assign("._ffm_settings", S, pos = 1)
}
