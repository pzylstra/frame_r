#' Calculates ignition temperature from proportion of silica free ash.
#' 
#' @param prop Value(s) in the interval (0, 1] for the proportion of 
#'   silica free ash, 
#'   
#' @return Calculated ignition temperature(s) in degrees Celsius.
#' 
#' @export
#'
ignitionTempPSFA <- function(prop) {
  if (prop <= 0 || prop > 1)
    stop("Values for prop must be proportions in the interval (0,1]")
  
  # natural log of percent silica free ash
  logPc <- log(prop * 100)
  
  # temperature calculation
  354.0 - (13.9 * logPc) - (2.91 * logPc * logPc)
}
