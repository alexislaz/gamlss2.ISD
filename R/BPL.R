#' Individual Size Distribution for fitting a GAMLSS
#'
#' `BPL()` defines the bounded Power-Law distribution within the `gamlss` environment (as the equivalent of `bPL()` in `gamlss2`)
#'
#' @param counts Number of observations for a corresponding size
#' @param lower Lowest size limit (does not need to equal `min(size)`)
#' @param upper Highest size limit (does not need to equal `max(size)`)
#' @param mu.link The exponent of ISD
#'
#' @return A `gamlss.family` object
#'
#' @rdname BPL-gamlss
#' @export
#'
#' @examples
#'
#' d = data.frame(size = rbPL(1e4))
#'
#' # estimate the exponent 'b'
#' # using 'gamlss2'
#' m2 = gamlss2(size ~ 1,
#'              data = d,
#'              family = bPL(counts = 1, lower = 1, upper = 1000))
#' summary(m2)

#' # using 'gamlss'
#' m = gamlss(size ~ 1,
#'            data = d,
#'            family = BPL(counts = 1, lower = 1, upper = 1000))
#' summary(m)
#'
#' @references Edwards, A.M., Robinson, J.P.W., Plank, M.J., Baum, J.K. and Blanchard, J.L. (2017), Testing and recommending methods for fitting size spectra to data. Methods Ecol Evol, 8: 57-67. https://doi.org/10.1111/2041-210X.12641
#' @references Edwards, A. M., Robinson, J. P. W., Blanchard, J. L., Baum, J. K., & Plank, M. J. (2020). Accounting for the bin structure of data removes bias when fitting size spectra. Marine Ecology Progress Series, 636, 19–33. https://www.jstor.org/stable/26920653
#'
#'
#'
BPL = function (counts, lower, upper, mu.link ="identity")
{
  mstats = checklink("mu.link", "bounded Power Law", substitute(mu.link), c("identity"))

  env = new.env(parent = .GlobalEnv)
  assign("counts", counts, envir = env)
  assign("lower", lower, envir = env)
  assign("upper", upper, envir = env)

  dldm = eval(substitute({
    function(y, mu) {
      if(NCOL(y) == 1) {
        .Call("score_bpl", as.numeric(y), as.numeric(mu),
              as.numeric(counts),
              as.numeric(lower),
              as.numeric(upper),
              PACKAGE = "gamlss2.ISD")
      }
      else {
        .Call("score_bpl_binned", as.numeric(y[, 1]), as.numeric(y[, 2]),
              as.numeric(mu),
              as.numeric(counts),
              as.numeric(lower),
              as.numeric(upper),
              PACKAGE = "gamlss2.ISD")
      }
    }
  }, list(
    counts = get("counts"),
    lower = get("lower"),
    upper = get("upper")
  )))
  #environment(dldm) = env


  d2ldm2 = eval(substitute({
    function(y, mu) {
      if(NCOL(y) == 1) {
        -1 * .Call("hess_bpl", as.numeric(y), as.numeric(mu),
                   as.numeric(counts),
                   as.numeric(lower),
                   as.numeric(upper),
                   PACKAGE = "gamlss2.ISD")
      }
      else {
        -1 * .Call("hess_bpl_binned", as.numeric(y[, 1]), as.numeric(y[, 2]),
                   as.numeric(mu),
                   as.numeric(counts),
                   as.numeric(lower),
                   as.numeric(upper),
                   PACKAGE = "gamlss2.ISD")
      }
    }
  }, list(
    counts = get("counts"),
    lower = get("lower"),
    upper = get("upper")
  )))
  #environment(d2ldm2) = env

  G.dev.incr = eval(substitute({
    function(y, mu) {
      if(NCOL(y) == 1) {
        -2 * .Call("loglik_bpl", as.numeric(y), as.numeric(mu),
                   as.numeric(counts),
                   as.numeric(lower),
                   as.numeric(upper),
                   PACKAGE = "gamlss2.ISD")
      }
      else {
        -2 * .Call("loglik_bpl_binned", as.numeric(y[, 1]), as.numeric(y[, 2]),
                   as.numeric(counts),
                   as.numeric(lower),
                   as.numeric(upper),
                   PACKAGE = "gamlss2.ISD")
      }
    }
  }, list(
    counts = get("counts"),
    lower = get("lower"),
    upper = get("upper")
  )))
  #environment(G.dev.incr) = env

  y.valid = eval(substitute({
    function(y) {
      if(!is.numeric(y))
        stop("the response should be numeric")
      if(anyNA(y))
        stop("`NA`s in response; better to pass a dataset without `NA`s to 'gamlss'")
      if((!length(as.numeric(counts))) || (!length(as.numeric(lower))) || (!length(as.numeric(upper))) ||
         anyNA(as.numeric(counts)) || anyNA(as.numeric(lower)) || anyNA(as.numeric(upper)))
        stop("invalid 'counts'/'lower'/'upper'")
      if(any(as.numeric(lower) > as.numeric(upper)))
        stop("need lower < upper")
      if(any((y < as.numeric(lower)) | (y > as.numeric(upper))))
        stop("response should be in [lower, upper]")
      return(TRUE)
    }
  }, list(
    counts = get("counts"),
    lower = get("lower"),
    upper = get("upper")
  )))
  #environment(y.valid) = env

  structure(
    list(family = c("BPL","bounded Power Law"),
         parameters = list(mu = TRUE),
         nopar = 1,
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         mu.linkfun = mstats$linkfun,
         mu.linkinv = mstats$linkinv,
         mu.dr = mstats$mu.eta,
         #dldm = function(y,mu) ((y-mu)/mu^2),
         dldm = dldm,
         #d2ldm2 = function(mu) (-1/mu^2),
         d2ldm2 = d2ldm2,
         #G.dev.incr  = function(y, mu, ...)  -2*dEXP(x = y, mu = mu, log = TRUE),
         G.dev.incr = G.dev.incr,
         #rqres = expression(rqres(pfun = "pbPL2", type = "Continuous", y = y, mu = mu)),
         rqres = as.expression(substitute({
           rqres(pfun = "pBPL", type = "Continuous",
                 y = y, mu = mu, lower = lower, upper = upper)
         }, list(
           counts = get("counts", envir = env),
           lower = get("lower", envir = env),
           upper = get("upper", envir = env)
         ))),

         mu.initial = as.expression(
           substitute({
             mu = if(NCOL(y) == 1) {
               .Call("init_bpl", as.numeric(y),
                     as.numeric(counts),
                     as.numeric(lower),
                     as.numeric(upper),
                     PACKAGE = "gamlss2.ISD")
             }
             else {
               .Call("init_bpl_binned", as.numeric(y[, 1]), as.numeric(y[, 2]),
                     as.numeric(counts),
                     as.numeric(lower),
                     as.numeric(upper),
                     PACKAGE = "gamlss2.ISD")
             }
           }, list(
             counts = get("counts", envir = env),
             lower = get("lower", envir = env),
             upper = get("upper", envir = env)
           ))
         ),
         #mu.valid = function(mu) all(mu > 0),
         mu.valid = function(mu) {
           return(TRUE)
         },
         #y.valid = function(y) all(y > 0),
         y.valid = y.valid,
         mean = function(mu) mu,
         variance = function(mu) mu ^ 2
    ),
    class = c("gamlss.family", "family"))
}

#' @rdname BPL-gamlss
#' @export
#----------------------------------------------------------------------------------------
dBPL = function(x = 1, mu = -2, lower = 1, upper = 100, log = FALSE)
{
  dbPL(x = x, b = mu, lower = lower, upper = upper, log = log)
}


#' @rdname BPL-gamlss
#' @export
#----------------------------------------------------------------------------------------
pBPL = function(q = 1, mu = -2, lower = 1, upper = 100, log = FALSE)
{
  pbPL(q = q, b = mu, lower = lower, upper = upper, log = log)
}

#' @rdname BPL-gamlss
#' @export
#----------------------------------------------------------------------------------------
qBPL = function(p = 0.1, mu = -2, lower = 1, upper = 100)
{
  qbPL(p = p, b = mu, lower = lower, upper = upper)
}


#' @rdname BPL-gamlss
#' @export
#----------------------------------------------------------------------------------------
rBPL = function(n = 1, mu = -2, lower = 1, upper = 100)
{
  rbPL(n = n, b = mu, lower = lower, upper = upper)
}


