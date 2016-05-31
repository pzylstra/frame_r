# Apply changes to rows in a data frame which meet some condition.
# 
# Similar to \code{\link[dplyr]{mutate}} but only applies changes changes to
# those rows which satisfy the given condition.
# 
# Code for this function is taken from 
# \href{http://stackoverflow.com/a/34096575/40246}{Stackoverflow}
# and was written by G. Grothendieck.
# 
# @param .data Input data table.
# 
# @param condition Expression returning \code{TRUE} or \code{FALSE} to
#   restrict changes to desired rows.
# 
# @param ... Name-value pairs of expressions.
# 
# @examples 
# \dontrun{
# tbl <- tbl %>% 
#   mutate_cond(species == "1" & param == "leafLength", value = "100", units = "mm")
# }
# 
# 
.mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}