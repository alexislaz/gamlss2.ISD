#' Individual Size Distribution for fitting a GAMLSS
#'
#' `bPL` returns a distribution function for ISD described as a bounded Power-Law distribution (bPL) in
#' Edwards et al. 2017. `bPL` is defined by only one parameter (defined as `b` here) which is equivalent to 'b'
#' as described in Edwards et al. 2017. When fitting with `gamlss2`, an estimate for `b` is returned.
#'
#' @param counts Number of observations for a corresponding size
#' @param lower Lowest size limit (does not need to equal `min(size)`)
#' @param upper Highest size limit (does not need to equal `max(size)`)
#'
#' @return A `gamlss2.family` object
#'
#' @export
#'
#' @examples
#'
#' # For these examples we need the `sizeSpetra` package.
#' # remotes::install_github("andrew-edwards/sizeSpectra")
#' # Note, though, that `sizeSpectra` is not needed for `gamlss2.ISD` in general.
#'
#' # simulate a dataset of 1000 sizes coming from an ISD with
#' # b = -2, and range of sizes bounded between 1 and 1000
#'
#' # 1) either use sizeSpectra::rPLB
#' d = data.frame(size = sizeSpectra::rPLB(1e4, -2, 1, 1000))
#' # 2) or, alternatively, use `bPL()$r` (by accessing the `$r` function of `bPL`)
#' d = data.frame(size = bPL(1, 1, 1000)$r(1e4, list(b = -2)))
#'
#' # estimate the exponent 'b'
#' m = gamlss2(size ~ 1,
#'             data = d,
#'             family = bPL(counts = 1, lower = 1, upper = 1000))
#'
#' summary(m)
#'
#' # If data is in binned form; i.e.:
#' bin.d = sizeSpectra::binData(d$size, binWidth = "2k")$binVals
#'
#' # We can, also, use a 2-column matrix for the response:
#' bin.m = gamlss2(cbind(binMin, binMax) ~ 1,
#'                 data = bin.d,
#'                 family = bPL(counts = bin.d$binCount, lower = min(bin.d$binMin), upper = max(bin.d$binMax)))
#'
#'
#' summary(bin.m)
#'
#' @references Edwards, A.M., Robinson, J.P.W., Plank, M.J., Baum, J.K. and Blanchard, J.L. (2017), Testing and recommending methods for fitting size spectra to data. Methods Ecol Evol, 8: 57-67. https://doi.org/10.1111/2041-210X.12641
#' @references Edwards, A. M., Robinson, J. P. W., Blanchard, J. L., Baum, J. K., & Plank, M. J. (2020). Accounting for the bin structure of data removes bias when fitting size spectra. Marine Ecology Progress Series, 636, 19–33. https://www.jstor.org/stable/26920653
#'

bPL = function(counts, lower, upper)
{
  fam = list(
    "family" = "bPL",
    "names" = "b",
    "links" = c("b" = "identity"),
    "type" = "continuous",
    "score" = list(
      "b" = function(y, par, ...) {
        if(NCOL(y) == 1) {
          .Call("score_bpl",
                as.numeric(y),
                as.numeric(par$b),
                as.numeric(counts),
                as.numeric(lower), as.numeric(upper))
        } else {
          .Call("score_bpl_binned",
                as.numeric(y[, 1]), as.numeric(y[, 2]),
                as.numeric(par$b),
                as.numeric(counts),
                as.numeric(lower), as.numeric(upper))
        }
      }
    ),
    "hess" = list(
      "b" = function(y, par, ...) {
        if(NCOL(y) == 1) {
          .Call("hess_bpl",
                as.numeric(y),
                as.numeric(par$b),
                as.numeric(counts),
                as.numeric(lower), as.numeric(upper))
        } else {
          .Call("hess_bpl_binned",
                as.numeric(y[, 1]), as.numeric(y[, 2]),
                as.numeric(par$b),
                as.numeric(counts),
                as.numeric(lower), as.numeric(upper))
        }
      }
    ),
    "loglik" = function(y, par, ...) {
      if(NCOL(y) == 1) {
        sum(.Call("loglik_bpl",
                  as.numeric(y),
                  as.numeric(par$b),
                  as.numeric(counts),
                  as.numeric(lower), as.numeric(upper)),
            na.rm = TRUE)
      } else {
        sum(.Call("loglik_bpl_binned",
                  as.numeric(y[, 1]), as.numeric(y[, 2]),
                  as.numeric(par$b),
                  as.numeric(counts),
                  as.numeric(lower), as.numeric(upper)),
            na.rm = TRUE)
      }
    },
    # needs "logLik" instead of "loglik" in newer versions of 'gamlss2'
    "logLik" = function(y, par, ...) {
      if(NCOL(y) == 1) {
        sum(.Call("loglik_bpl",
                  as.numeric(y),
                  as.numeric(par$b),
                  as.numeric(counts),
                  as.numeric(lower), as.numeric(upper)),
            na.rm = TRUE)
      } else {
        sum(.Call("loglik_bpl_binned",
                  as.numeric(y[, 1]), as.numeric(y[, 2]),
                  as.numeric(par$b),
                  as.numeric(counts),
                  as.numeric(lower), as.numeric(upper)),
            na.rm = TRUE)
      }
    },
    "b" = function(par, ...) {
      par$b
    },
    # "d" = ...
    "pdf" = function(y, par, log = FALSE, ...) {
      #stop("used 'd'")
      if(NCOL(y) == 2) y = y[, 1] #y = (y[, 1] + y[, 2]) / 2
      .Call("d_bpl",
            as.numeric(y),
            as.numeric(par$b),
            as.numeric(lower), as.numeric(upper),
            as.logical(log))
    },
    # "p" = ...
    "cdf" = function(y, par, ...) {
      if(NCOL(y) == 1) {
        .Call("p_bpl",
              as.numeric(y),
              as.numeric(par$b),
              as.numeric(lower), as.numeric(upper))
      } else {
        .Call("p_bpl_binned",
              as.numeric(y[, 1]), as.numeric(y[, 2]),
              as.numeric(par$b),
              as.numeric(lower), as.numeric(upper))
      }
    },
    # "r" = ...
    "random" = function(n, par, ...) {
      rand = stats::runif(n)
      .Call("r_bpl",
            as.integer(n),
            as.numeric(par$b),
            as.numeric(lower), as.numeric(upper),
            as.numeric(rand))
    },
    # "q" = ...
    "quantile" = function(p, par, ...) {
      .Call("q_bpl",
            as.numeric(p),
            as.numeric(par$b),
            as.numeric(lower), as.numeric(upper))
    },
    "initialize" = list(
      "b" = function(y, ...) {
        if(NCOL(y) == 1) {
          .Call("init_bpl",
                as.numeric(y),
                as.numeric(counts),
                as.numeric(lower), as.numeric(upper))
        } else {
          .Call("init_bpl_binned",
                as.numeric(y[, 1]), as.numeric(y[, 2]),
                as.numeric(counts),
                as.numeric(lower), as.numeric(upper))
        }
      }
    ),
    "mean" = function(par) par$b,
    "valid.response" = function(x) {
      #if(NCOL(x) == 1) cat("`bPL`: response is individual sizes (1-d)\n") else cat("`bPL`: response is binned (2-d)\n")
      if(!is.numeric(x)) stop("the response should be numeric")
      #if(any(x == 0)) stop("response equal to 0 detected")
      if(anyNA(x)) stop("`NA`s in response; better to pass a dataset without `NA`s to 'gamlss2'")
      if((!length(counts)) || (!length(lower)) || (!length(upper)) || anyNA(counts) || anyNA(lower) || anyNA(upper)) stop("invalid 'counts'/'lower'/'upper'")
      if(any(lower > upper)) stop("need lower < upper")
      if(any((x < lower) | (x > upper))) stop("response should be in [lower, upper]")
      #if(any(counts == 0)) stop("0 counts detected")
      return(TRUE)
    }
  )

  class(fam) = "gamlss2.family"

  return(fam)
}
