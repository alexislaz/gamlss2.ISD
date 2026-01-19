#' bounded Power-Law
#'
#' Density, Distribution function, quantile function and, random generation.
#' Mostly based on `sizeSpectra` package; they have been re-written here to capture all essentials for `gamlss2` functionality and exported for convenience.
#'
#' @param x Vector of values to compute density
#' @param q Vector of values to compute distribution function
#' @param p Vector of probabilities
#' @param b Exponent of ISD
#' @param lower Lowest size limit
#' @param upper Highest size limit
#' @param log Use `log`?
#'
#' @return dbPL returns pdf. pbPL returns cumulative distribution values. qbPL gives the quantiles. rbPL generates random numbers.
#'
#' @export
#'
#' @examples
#'
#' # none
#'
#' @references Edwards, A.M., Robinson, J.P.W., Plank, M.J., Baum, J.K. and Blanchard, J.L. (2017), Testing and recommending methods for fitting size spectra to data. Methods Ecol Evol, 8: 57-67. https://doi.org/10.1111/2041-210X.12641
#' @references Edwards, A. M., Robinson, J. P. W., Blanchard, J. L., Baum, J. K., & Plank, M. J. (2020). Accounting for the bin structure of data removes bias when fitting size spectra. Marine Ecology Progress Series, 636, 19–33. https://www.jstor.org/stable/26920653
#'


dbPL =function(x = 1, b = -2, lower = 1, upper = 100, log = FALSE)
{
  if(NCOL(x) == 2) {
    if(log) {
      .Call("d_bpl_binned_log",
            as.numeric(x[, 1]), as.numeric(x[, 2]),
            as.numeric(b),
            as.numeric(lower), as.numeric(upper))
    } else {
      .Call("d_bpl_binned",
            as.numeric(x[, 1]), as.numeric(x[, 2]),
            as.numeric(b),
            as.numeric(lower), as.numeric(upper))
    }
  } else {
    if(log) {
      .Call("d_bpl_log",
            as.numeric(x),
            as.numeric(b),
            as.numeric(lower), as.numeric(upper))
    } else {
      .Call("d_bpl",
            as.numeric(x),
            as.numeric(b),
            as.numeric(lower), as.numeric(upper))
    }
  }
}


pbPL = function(q = 1, b = -2, lower = 1, upper = 100, log = FALSE)
{
  if(NCOL(q) == 1) {
    if(log) {
      .Call("p_bpl_log",
            as.numeric(q),
            as.numeric(b),
            as.numeric(lower), as.numeric(upper))
    } else {
      .Call("p_bpl",
            as.numeric(q),
            as.numeric(b),
            as.numeric(lower), as.numeric(upper))
    }
  } else {
    if(log) {
      .Call("p_bpl_binned_log",
            as.numeric(q[, 1]), as.numeric(q[, 2]),
            as.numeric(b),
            as.numeric(lower), as.numeric(upper))
    } else {
      .Call("p_bpl_binned",
            as.numeric(q[, 1]), as.numeric(q[, 2]),
            as.numeric(b),
            as.numeric(lower), as.numeric(upper))
    }
  }
}

qbPL = function(p = 0.1, b = -2, lower = 1, upper = 100)
{
  .Call("q_bpl",
        as.numeric(p),
        as.numeric(b),
        as.numeric(lower), as.numeric(upper))
}

rbPL = function(n = 1, b = -2, lower = 1, upper = 100)
{
  rand = stats::runif(n)
  .Call("r_bpl",
        as.integer(n),
        as.numeric(b),
        as.numeric(lower), as.numeric(upper),
        as.numeric(rand))
}




