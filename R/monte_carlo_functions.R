#' Randomly perturbs values for selected numeric parameters in a table.
#' 
#' Given an input parameters table, this function adds random noise to 
#' specified numeric parameters and returns an updated copy of the table.
#' Noise can be based on a Normal (Gaussian) or Beta distribution. The 
#' maximum variation is controlled by the \code{max.prop} argument. 
#' When adding Gaussian noise, a random value is drawn from a Normal 
#' distribution with a mean of zero and standard deviation equal to
#' \code{max.prop / 3}. This results in 99% of percent of random values
#' being within plus or minus \code{max.prop} times the input value.
#' Any random perturbations outside this range will be discarded and 
#' re-drawn (so effectively working with a truncated Normal distribution).
#' 
#' @export
#' 
ffm_param_variance <- function(tbl,
                               max.prop = 0.1,
                               labels = c("leafWidth", "leafLength", "leafThickness", "leafSeparation"),
                               dist = c("normal", "beta"),
                               beta.a = 1,
                               beta.b = 1) {
  
  if (!forestflames:::.is_param_table(tbl))
    stop("Input table must be a character matrix or data frame with 4 or 5 columns")
  
  if (max.prop <= 0.0 || max.prop >= 1.0) {
    warning("max.prop should be a value between 0 and 1")
    tbl
  }
  
  ii <- tbl$param %in% labels
  if (!any(ii)) {
    warning("The input table does contain any rows with the specified parameter labels")
    tbl
  }
  else {
    dist = match.arg(dist)
    if (dist == "normal") .add_gaussian_noise(tbl, ii, max.prop, clamp)
    else .add_beta_noise(tbl, ii, max.prop, beta.a, beta.b)
  }
}

.add_gaussian_noise <- function(tbl, row.indices, max.prop, clamp) {
  xs <- as.numeric(tbl$value[row.indices])
  if (any(is.na(xs))) {
    warning("Value(s) with missing data not altered")
  }
  
  ii <- !is.na(xs)
  if (any(ii)) {
    # Set standard deviations to contain maximum scatter
    # 99% of the time
    limits <- xs[ii] * max.prop
    dxs <- rnorm(sum(ii), mean = 0, sd = limits / 3)
    
    out.of.range <- dxs > limits
    while (any(out.of.range)) {
      dxs[out.of.range] <- rnorm(sum(out.of.range), mean = 0, sd = limits / 3)
      out.of.range <- dxs > limits
    }
    
    xs[ii] <- xs[ii] + dxs
    
    tbl$value[row.indices] <- as.character(xs)
  }
  
  tbl
}


.add_beta_noise <- function(tbl, ii, max.prop, beta.a, beta.b) {
  stop("Beta option not working yet")
  
  xs <- as.numeric(tbl$value[row.indices])
  if (any(is.na(xs))) {
    warning("Value(s) with missing data not altered")
  }

  # TODO - Write Me !
    
  tbl
}

