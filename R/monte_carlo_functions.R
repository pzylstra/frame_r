#' Generates example Beta distributions.
#' 
#' This function can be used to explore the shapes and parameters of
#' alternative Beta distributions. Given one or more values for
#' each of the two Beta parameters, it calculates the density
#' curve for each combination plus the mean, mode and variance of
#' the distribution.
#' 
#' Note that this function will draw and describe Beta distributions
#' with any valid parameter values, but the \code{\link{ffm_param_variance}}
#' function only accepts distributions where both parameters are greater
#' than one to ensure that the mode exists.
#' 
#' @param beta.a Vector of one or more values for parameter \code{a}
#'   (called \code{shape1} by R beta functions).
#'   
#' @param beta.b Vector of one or more values for parameter \code{b}
#'   (called \code{shape2} by R beta functions).
#'   
#' @param show.plot If \code{TRUE}, display the resulting distributions.
#' 
#' @return Invisibly returns a list with elements:
#'   \describe{
#'     \item{plot}{A ggplot object showing the distributions.}
#'     \item{summary}{A data frame summarizing the moments of the distributions.}
#'   }
#' 
#' @seealso \code{\link{ffm_param_variance}} which has an option to use 
#'   a Beta distribution to vary parameter values.
#' 
#' @export
#' 
#' @importFrom ggplot2 
#'   ggplot aes geom_line scale_colour_discrete theme theme_bw
#'   facet_wrap label_both
#' 
ffm_describe_betas <- function(beta.a = c(2, 5, 10),
                               beta.b = c(2, 5, 10),
                               show.plot = TRUE) {
    
  ab <- expand.grid(a = beta.a, b = beta.b)
  x <- seq(0, 1, 0.01)
  
  betas <- lapply(1:nrow(ab), function(i) {
    a <- ab[i, 1]
    b <- ab[i, 2]
    data.frame(a, b, x = x, density = dbeta(x, a, b))
  })
  
  betas <- do.call(rbind, betas)
  
  g <-ggplot(data = betas, aes(x = x, y = density, colour = factor(b))) +
    geom_line() +
    scale_colour_discrete(name = "b") +
    facet_wrap(~ a, labeller = label_both) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  if (show.plot) print(g)
  
  # Summarize means, modes etc
  summary <- ab %>% 
    mutate(mean = a / (a + b),
           mode = (a - 1) / (a + b - 2),
           variance = (a * b) / ( (a + b + 1) * (a + b)^2 )
    )
  
  invisible( list(plot = g, summary = summary) )
}

#' Randomly perturbs values for selected numeric parameters in a table.
#' 
#' Given an input parameters table, this function adds random noise to 
#' specified numeric parameters and returns an updated copy of the table.
#' Noise can be based on a Normal (Gaussian), Beta or Uniform distribution.
#' The maximum absolute perturbation is controlled by the \code{max.prop} 
#' argument. 
#' 
#' When \code{method = "normal"} (the default), random deviates are drawn
#' from a Normal distribution with mean zero and standard deviation equal
#' to \code{x * max.prop / 3} where \code{x} is the input parameter value. 
#' This results in 99% of deviates lying within the allowed range. Any deviates
#' beyond this range are discarded and re-drawn, effectively rendering the distribution
#' a (slightly) truncated Normal.
#' 
#' When \code{method = "beta"}, it is possible to generate both symmetric and
#' skewed deviates depending on the values of the distribution parameters
#' \code{beta.a} and \code{beta.b}. The function constrains the choice of
#' parameters to both be greater than 1.0 so that a mode exists for the 
#' distribution, and will issue an error message if this is not the case.
#' Random deviates between 0 and 1 are drawn from the distribution and scaled
#' such that the mode of the distribution corresponds to the input parameter
#' value. For a symmetric distribution (\code{beta.a == beta.b}) a deviate of
#' 0 is scaled to \code{-max.prop * x} and a deviate of 1 is scaled to 
#' \code{max.prop * x}, where \code{x} is the input parameter value. Skewed
#' Beta distributions (\code{beta.a != beta.b}) give correspondingly skewed
#' distributions of perturbations. Random deviates drawn from the
#' distribution are scaled such that those lying on the side of the mode with
#' the longest tail can approach the full range of pertubation, while those
#' lying on the other side have a reduced range.
#' 
#' When \code{method = "uniform"}, perturbations are uniformly distributed
#' either side of the input parameter value. 
#' 
#' @param tbl Input parameter table.
#' 
#' @param max.prop The maximum absolute perturbation to apply, expressed
#'   as a proportion between 0 and 1.
#' 
#' @param labels Labels of the parameters whose values are to be
#'   perturbed. It only makes sense to specify parameters with numeric values.
#'   
#' @param method Method used to generate perturbations.
#'   \describe{
#'     \item{normal}{Draw perturbations from a Normal distribution,
#'       with zero mean and where max.prop corresponds to 3 standard
#'       deviations.}
#'     \item{beta}{Draw random deviates from a modal Beta distribution
#'       (ie. beta.a and beta.b both greater than 1) and scale
#'       them to accord with \code{max.prop} }.
#'     \item{uniform}{Draw perturbations from a Uniform distribution}.
#'   }
#'   
#' @param beta.a First shape parameter for the Beta distribution. Ignored
#'   if \code{method} is not \code{"beta"}.
#' 
#' @param beta.b Second shape parameter for the Beta distribution. Ignored
#'   if \code{method} is not \code{"beta"}.
#'   
#' @return A copy of the input parameter table with selected values updated.
#' 
#' @examples
#' \dontrun{
#' # Vary leaf traits for all species with Uniform perturbations which are,
#' # at most, 20% of the input values
#' tbl.updated <- ffm_param_variance(tbl, max.prop = 0.2, method = "unif")
#' 
#' # Vary ignition temperatures by up to 10% for all species based on Normal 
#' # perturbations
#' tbl.updated <- ffm_param_variance(tbl, labels = "ignitionTemp", max.prop=0.1)
#' 
#' # Use with the pipe operator to vary species and site parameters
#' # separately (requires either the dplyr or maggritr package)
#' 
#' tbl.updated <- tbl %>%
#'   # vary leaf traits uniformly
#'   ffm_param_variance(max.prop = 0.2, method = "uniform") %>%
#'   
#'   # vary fuel load with a skew
#'   ffm_param_variance(labels = "fuelLoad", max.prop = 0.5,
#'                      method = "beta", beta.a = 2, beta.b = 5)
#' }
#' 
#' @export
#' 
ffm_param_variance <- function(tbl,
                               max.prop = 0.1,
                               labels = c("leafWidth", "leafLength", "leafThickness", "leafSeparation"),
                               method = c("normal", "beta", "uniform"),
                               beta.a = 2,
                               beta.b = 2) {
  
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
    method = match.arg(method)
    
    if (method == "normal") 
      .add_gaussian_noise(tbl, ii, max.prop)
    else if (method == "beta")
      .add_beta_noise(tbl, ii, max.prop, beta.a, beta.b)
    else if (method == "uniform")
      .add_uniform_noise(tbl, ii, max.prop)
  }
}

.add_gaussian_noise <- function(tbl, row.flags, max.prop) {
  xs <- as.numeric(tbl$value[row.flags])
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
    
    tbl$value[row.flags] <- as.character(xs)
  }
  
  tbl
}


.add_beta_noise <- function(tbl, row.flags, max.prop, beta.a, beta.b) {
  if (!(beta.a > 1 & beta.b > 1))
    stop("Both Beta parameters must be greater than 1")
  
  # modal value and length of longest tail
  M <- (beta.a - 1) / (beta.a + beta.b - 2)
  Tail <- max(M, 1 - M)
  
  xs <- as.numeric(tbl$value[row.flags])
  if (any(is.na(xs))) {
    warning("Value(s) with missing data not altered")
  }

  ii <- !is.na(xs)
  if (any(ii)) {
    dzs <- rbeta(sum(ii), beta.a, beta.b)
    dxs <- xs[ii] * max.prop * (dzs - M) / Tail
    xs[ii] <- xs[ii] + dxs
    tbl$value[row.flags] <- as.character(xs)
  }
    
  tbl
}


.add_uniform_noise <- function(tbl, row.flags, max.prop) {
  xs <- as.numeric(tbl$value[row.flags])
  if (any(is.na(xs))) {
    warning("Value(s) with missing data not altered")
  }
  
  ii <- !is.na(xs)
  if (any(ii)) {
    dxs <- xs[ii] * runif(sum(ii), -max.prop, max.prop)
    xs[ii] <- xs[ii] + dxs
    tbl$value[row.flags] <- as.character(xs)
  }
  
  tbl
}

