# framecpp validation bridge (Phase 6, dev/package-boundary.md in the framecpp repo).
#
# Retargets frame's single seam onto the in-process C++ engine `framecpp` instead
# of the Scala JVM + SQLite path. `ffm_run(..., engine = "framecpp")` builds an
# framecpp `fire_scenario` from the same paramBuilder table `CSVRunner` consumes,
# runs the engine, and writes the identical seven-table SQLite result database at
# `db.path`, so every downstream frame analysis (Ensembles / Risk / Fauna /
# Earth / Summaries, all of which read the db back via `ffm_db_load()`) runs
# unchanged and can be diffed against the JVM.
#
# This is a *validation bridge*, kept behind a flag with the JVM path intact and
# the default. It is local / branch-only and NOT pushed to pzylstra/frame_r
# without Zylstra sign-off.

# paramBuilder levelName (lower-case) -> framecpp StratumLevel code.
.FRAMECPP_LEVEL_NUM <- c("near surface" = 1L, "elevated" = 2L,
                       "midstorey" = 3L, "canopy" = 4L)
# frame level name (as framecpp emits it) -> StratumLevel code.
.FRAMECPP_LEVEL_CODE <- c(Surface = 0L, NearSurface = 1L, Elevated = 2L,
                        MidStorey = 3L, Canopy = 4L)

#' Run a FRaME simulation with the framecpp in-process engine
#'
#' Validation-bridge alternative to \code{\link{ffm_run}}'s JVM path: builds an
#' framecpp \code{fire_scenario} from \code{params}, runs the C++ engine in
#' process, and writes the same seven-table SQLite results database at
#' \code{db.path} that the Scala \code{CSVRunner} would. Requires the
#' \pkg{framecpp} package (>= 0.1.0.9002, for the full \code{fire_result}).
#'
#' Species in the \code{IgnitionPaths} table are resolved back to their names
#' from the parameter table (framecpp's engine identifies species by id).
#'
#' @param params A \code{paramBuilder} parameter table (as for \code{\link{ffm_run}}).
#' @param db.path Output SQLite database path.
#' @param default.species.params Optional defaults used to complete \code{params}
#'   (as in \code{\link{ffm_run}}).
#' @param db.recreate If \code{TRUE}, recreate the database; otherwise append,
#'   assigning the next \code{repId} exactly as the JVM runner does.
#' @return \code{TRUE} on success (invisibly writes the database).
#' @export
ffm_run_framecpp <- function(params, db.path,
                           default.species.params = NULL,
                           db.recreate = FALSE) {
  if (!requireNamespace("framecpp", quietly = TRUE))
    stop("engine = \"framecpp\" requires the 'framecpp' package to be installed")

  db.path <- normalizePath(db.path, mustWork = FALSE)

  if (!ffm_check_params(params, quiet = TRUE)) {
    if (!is.null(default.species.params)) {
      params <- ffm_complete_params(params, default.species.params)
    } else {
      warning("Parameters incomplete and no defaults provided")
      return(FALSE)
    }
  }

  tables <- .framecpp_result_tables(params)
  .framecpp_write_db(tables, db.path, db.recreate)
  TRUE
}

# Ordered species names per stratum level (framecpp build order = strata sorted by
# id, species sorted by id within a stratum). Keyed by StratumLevel code.
.framecpp_species_by_level <- function(params) {
  p <- params
  p$value <- as.character(p$value)
  strataIds <- sort(unique(p$stratum[!is.na(p$stratum)]))
  out <- list()
  for (sid in strataIds) {
    inStr <- !is.na(p$stratum) & p$stratum == sid
    meta  <- inStr & is.na(p$species)
    levelName <- tolower(trimws(p$value[meta & p$param == "levelName"]))
    code <- unname(.FRAMECPP_LEVEL_NUM[levelName])
    spIds <- sort(unique(p$species[inStr & !is.na(p$species)]))
    nms <- vapply(spIds, function(spid) {
      v <- p$value[inStr & !is.na(p$species) & p$species == spid & p$param == "name"]
      if (length(v) == 1L) v else as.character(spid)
    }, character(1))
    out[[as.character(code)]] <- nms
  }
  out
}

# The seven frame result tables, in frame's schema, from a framecpp fire_result.
# Only IgnitionPaths differs from framecpp's native tables: its engine speciesId /
# speciesRank are replaced by frame's species name.
.framecpp_result_tables <- function(params) {
  scn <- framecpp::scenario_from_params(params)
  res <- framecpp::run_fire_model(scn)
  t <- res$tables

  namesByLevel <- .framecpp_species_by_level(params)
  ip <- t$IgnitionPaths
  if (nrow(ip)) {
    codes <- unname(.FRAMECPP_LEVEL_CODE[ip$level])
    spname <- vapply(seq_len(nrow(ip)), function(i)
      namesByLevel[[as.character(codes[i])]][ip$speciesRank[i] + 1L], character(1))
    IgnitionPaths <- data.frame(
      repId = ip$repId, runIndex = ip$runIndex, level = ip$level,
      pathType = ip$pathType, species = spname, segIndex = ip$segIndex,
      x0 = ip$x0, y0 = ip$y0, x1 = ip$x1, y1 = ip$y1,
      length = ip$length, flameLength = ip$flameLength, stringsAsFactors = FALSE)
  } else {
    IgnitionPaths <- data.frame(
      repId = integer(), runIndex = integer(), level = character(),
      pathType = character(), species = character(), segIndex = integer(),
      x0 = double(), y0 = double(), x1 = double(), y1 = double(),
      length = double(), flameLength = double(), stringsAsFactors = FALSE)
  }

  list(Sites = t$Sites, Runs = t$Runs, Strata = t$Strata,
       FlameSummaries = t$FlameSummaries, SurfaceResults = t$SurfaceResults,
       IgnitionPaths = IgnitionPaths, ROS = t$ROS)
}

# Write (or append) the seven tables to a SQLite db, matching the JVM runner's
# repId assignment: a fresh db starts at repId 1; appending uses max(repId) + 1.
.framecpp_write_db <- function(tables, db.path, db.recreate) {
  if (db.recreate && file.exists(db.path)) file.remove(db.path)
  con <- DBI::dbConnect(RSQLite::SQLite(), db.path)
  on.exit(DBI::dbDisconnect(con))

  repId <- 1L
  if ("Sites" %in% RSQLite::dbListTables(con)) {
    mx <- DBI::dbGetQuery(con, "SELECT MAX(repId) AS m FROM Sites")$m
    if (length(mx) && !is.na(mx)) repId <- as.integer(mx) + 1L
  }

  for (nm in names(tables)) {
    df <- tables[[nm]]
    if (nrow(df)) df$repId <- repId
    DBI::dbWriteTable(con, nm, df, append = TRUE)
  }
  invisible(repId)
}
