#' Individual Size Distribution for fitting a GAMLSS
#'
#' `bPL` returns a distribution function for ISD described as a bounded Power-Law distribution (bPL) in
#' Edwards et al. 2017. `bPL` is defined by only one parameter (defined as `mu` here) which is equivalent to 'b'
#' as described in Edwards et al. 2017. When fitting with `gamlss2`, estimate for `mu` is returned.
#'
#' @param counts Number of observations for a corresponding size
#' @param lower Lowest size limit (does not need to equal `min(size)`)
#' @param upper Highest size limit (does not need to equal `max(size)`)
#'
#' @return A `gamlss2` family object
#'
#' @export
#'
#' @examples
#' NA
#'
#' @references Edwards, A.M., Robinson, J.P.W., Plank, M.J., Baum, J.K. and Blanchard, J.L. (2017), Testing and recommending methods for fitting size spectra to data. Methods Ecol Evol, 8: 57-67. https://doi.org/10.1111/2041-210X.12641
#' @references Edwards, A. M., Robinson, J. P. W., Blanchard, J. L., Baum, J. K., & Plank, M. J. (2020). Accounting for the bin structure of data removes bias when fitting size spectra. Marine Ecology Progress Series, 636, 19–33. https://www.jstor.org/stable/26920653
#'

bPL = function(counts, lower, upper)
{
  fam = list(
    "family" = "bPL",
    "names" = "mu",
    "links" = c("mu" = "identity"),
    "score" = list(
      "mu" = function(y, par, ...) {
        .Call("score_bpl", as.numeric(y), as.numeric(par$mu), as.numeric(counts), as.numeric(lower), as.numeric(upper))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        .Call("hess_bpl", as.numeric(y), as.numeric(par$mu), as.numeric(counts), as.numeric(lower), as.numeric(upper))
      }
    ),
    "loglik" = function(y, par, ...) {
      sum(.Call("loglik_bpl", as.numeric(y), as.numeric(par$mu), as.numeric(counts), as.numeric(lower), as.numeric(upper)))
    },
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      .Call("d_bpl", as.numeric(y), as.numeric(par$mu), as.numeric(lower), as.numeric(upper))
    },
    "p" = function(y, par, ...) {
      .Call("p_bpl", as.numeric(y), as.numeric(par$mu), as.numeric(lower), as.numeric(upper))
    },
    "r" = function(n, par) {
      rand = stats::runif(n)
      .Call("r_bpl", as.integer(n), as.numeric(par$mu), as.numeric(lower), as.numeric(upper), as.numeric(rand))
    },
    "q" = function(p, par) {
      .Call("q_bpl", as.numeric(p), as.numeric(par$mu), as.numeric(lower), as.numeric(upper))
    },
    "initialize" = list(
      "mu" = function(y, ...) {
        .Call("init_bpl", as.numeric(y), as.numeric(counts), as.numeric(lower), as.numeric(upper))
      }
    ),
    "mean" = function(par) par$mu,
    "valid.response" = function(x) {
      if(!is.numeric(x)) stop("the response should be numeric")
      if(any(x == 0)) stop("response equal to 0 detected")
      if((!length(counts)) || (!length(lower)) || (!length(upper)) || anyNA(counts) || anyNA(lower) || anyNA(upper)) stop("invalid 'counts'/'lower'/'upper'")
      if(any(lower > upper)) stop("need lower < upper")
      if(any((x < lower) | (x > upper))) stop("response should be in [lower, upper]")
      if(any(counts == 0)) stop("0 counts detected")
      return(TRUE)
    }
  )

  class(fam) = "gamlss2.family"

  return(fam)
}
