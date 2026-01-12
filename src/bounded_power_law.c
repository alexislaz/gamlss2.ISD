#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>


// macros to find max; used to find 'max(length(a), length(b), ...)' to emulate recycling behaviour
#define MAX2(x1, x2) (((x1) > (x2)) ? (x1) : (x2))
#define MAX3(x1, x2, x3) MAX2(MAX2((x1), (x2)), (x3))
#define MAX4(x1, x2, x3, x4) MAX2(MAX3((x1), (x2), (x3)), (x4))
#define MAX5(x1, x2, x3, x4, x5) MAX2(MAX4((x1), (x2), (x3), (x4)), (x5))
#define CHECK_NUMERIC_INPUT(x, fun, arg) if( !( (Rf_isReal((x)) || Rf_isInteger((x))) && (LENGTH(x) > 0) ) ) Rf_error("(in '%s') argument '%s' is not 'numeric' with 'length' >= 1\n", (fun), (arg))



// adapted from Edwards et al. 2017 (sizeSpectra::negLL.PLB.counts)
SEXP loglik_bpl(SEXP y, SEXP b, SEXP c, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y, "loglik_bpl", "y");
  CHECK_NUMERIC_INPUT(b, "loglik_bpl", "b");
  CHECK_NUMERIC_INPUT(c, "loglik_bpl", "c");
  CHECK_NUMERIC_INPUT(l, "loglik_bpl", "l");
  CHECK_NUMERIC_INPUT(u, "loglik_bpl", "u");

  int ny = LENGTH(y), nb = LENGTH(b), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX5(ny, nb, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pc = REAL(c), *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if(pb[i % nb] != -1) {
      pans[i] = pc[i % nc] * log((pb[i % nb] + 1)/(R_pow(pu[i % nu], pb[i % nb] + 1) - R_pow(pl[i % nl], pb[i % nb] + 1))) + pb[i % nb] * pc[i % nc] * log(py[i % ny]);
    } else {
      pans[i] = -pc[i % nc] * log(log(pu[i % nu]) - log(pl[i % nl])) - pc[i % nc] * log(py[i % ny]);
      //pans[i] = pc[i % nc] * (log(log(pl[i % nl]) - log(pu[i % nu])) + pb[i % nb] * log(py[i % ny]));
    }
  }

  UNPROTECT(1);
  return(ans);
}

// adapted from sizeSpectra::negLL.PLB.bins.species
SEXP loglik_bpl_binned(SEXP y1, SEXP y2, SEXP b, SEXP c, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y1, "loglik_bpl_binned", "y1");
  CHECK_NUMERIC_INPUT(y2, "loglik_bpl_binned", "y2");
  CHECK_NUMERIC_INPUT(b, "loglik_bpl_binned", "b");
  CHECK_NUMERIC_INPUT(c, "loglik_bpl_binned", "c");
  CHECK_NUMERIC_INPUT(l, "loglik_bpl_binned", "l");
  CHECK_NUMERIC_INPUT(u, "loglik_bpl_binned", "u");

  int ny = LENGTH(y1), nb = LENGTH(b), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX5(ny, nb, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py1 = REAL(y1), *py2 = REAL(y2), *pb = REAL(b), *pc = REAL(c), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
  double comp2, comp1;

  // adapted from https://github.com/andrew-edwards/sizeSpectra/blob/master/R/likelihood.R#L443

  for(int i = 0; i < N; i++) {
    if(pb[i % nb] != -1) {
      comp2 = pc[i % nc] * log( fabs( R_pow(py2[i % ny], pb[i % nb] + 1) - R_pow(py1[i % ny],  pb[i % nb] + 1) ) );
      comp1 = - pc[i % nc] * log( fabs( R_pow(pu[i % nu], pb[i % nb] + 1) - R_pow(pl[i % nl], pb[i % nb] + 1) ) );
      pans[i] = comp1 + comp2;
    } else {
      comp2 = pc[i % nc] * log( log(py2[i % ny]) - log(py1[i % ny]) );
      comp1 = - pc[i % nc] * log( log(pu[i % nu]) - log(pl[i % nl]) );
      pans[i] = comp1 + comp2;
    }
  }

  UNPROTECT(1);
  return(ans);
}


/*
 first derivative of log-likelihood for binned data

 ll_not_minus1 = expression(- c * log( abs( u^(b + 1) - l^(b + 1) ) ) + c * log( abs( y2^(b + 1) - y1^(b + 1) ) ))
 Deriv::Deriv(ll_not_minus1, "b")

 ll_minus1 = expression( - c * log( log(u) - log(l) ) + c * log( log(y2) - log(y1) ) )
 Deriv::Deriv(ll_minus1, "b")
 */
SEXP score_bpl_binned(SEXP y1, SEXP y2, SEXP b, SEXP c, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y1, "score_bpl_binned", "y1");
  CHECK_NUMERIC_INPUT(y2, "score_bpl_binned", "y2");
  CHECK_NUMERIC_INPUT(b, "score_bpl_binned", "b");
  CHECK_NUMERIC_INPUT(c, "score_bpl_binned", "c");
  CHECK_NUMERIC_INPUT(l, "score_bpl_binned", "l");
  CHECK_NUMERIC_INPUT(u, "score_bpl_binned", "u");

  int ny = LENGTH(y1), nb = LENGTH(b), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX5(ny, nb, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double e1, e2, e3, e4, e5, e6, e7;
  double *py1 = REAL(y1), *py2 = REAL(y2), *pb = REAL(b), *pc = REAL(c), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if(pb[i % nb] != -1) {
      e1 = pb[i % nb] + 1;
      e2 = R_pow(pl[i % nl], e1);
      e3 = R_pow(pu[i % nu], e1);
      e4 = R_pow(py1[i % ny], e1);
      e5 = R_pow(py2[i % ny], e1);
      e6 = e3 - e2;
      e7 = e5 - e4;
      pans[i] = pc[i % nc] * (sign(e7) * (e5 * log(py2[i % ny]) - e4 * log(py1[i % ny]))/fabs(e7) - sign(e6) * (e3 * log(pu[i % nu]) - e2 * log(pl[i % nl]))/fabs(e6));
    } else {
      pans[i] = 0;
    }
  }

  UNPROTECT(1);
  return(ans);
}





// first derivative of ISD log-likelihood for counts of individual sizes
// computed using Deriv::Deriv(<loglik_bpl> for = -1 and != -1, respectively)
SEXP score_bpl(SEXP y, SEXP b, SEXP c, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y, "score_bpl", "y");
  CHECK_NUMERIC_INPUT(b, "score_bpl", "b");
  CHECK_NUMERIC_INPUT(c, "score_bpl", "c");
  CHECK_NUMERIC_INPUT(l, "score_bpl", "l");
  CHECK_NUMERIC_INPUT(u, "score_bpl", "u");

  int ny = LENGTH(y), nb = LENGTH(b), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX5(ny, nb, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double e1, e2, e3;
  double *py = REAL(y), *pc = REAL(c), *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if(pb[i % nb] != -1) {
      e1 = 1 + pb[i % nb];
      e2 = R_pow(pl[i % nl], e1);
      e3 = R_pow(pu[i % nu], e1);
      pans[i] = pc[i % nc] * ((1 - e1 * (e3 * log(pu[i % nu]) - e2 * log(pl[i % nl]))/(e3 - e2))/e1 + log(py[i % ny]));
    } else {
      pans[i] = pc[i % nc] * log(py[i % ny]);
    }
  }

  UNPROTECT(1);
  return(ans);
}


/*
 * second derivative of ISD log-likelihood for binned data
 * Deriv::Deriv(ll_not_minus1, "b", nderiv = 2)
 * Deriv::Deriv(ll_minus1, "b", nderiv = 2)
 */
SEXP hess_bpl_binned(SEXP y1, SEXP y2, SEXP b, SEXP c, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y1, "hess_bpl_binned", "y1");
  CHECK_NUMERIC_INPUT(y2, "hess_bpl_binned", "y2");
  CHECK_NUMERIC_INPUT(b, "hess_bpl_binned", "b");
  CHECK_NUMERIC_INPUT(c, "hess_bpl_binned", "c");
  CHECK_NUMERIC_INPUT(l, "hess_bpl_binned", "l");
  CHECK_NUMERIC_INPUT(u, "hess_bpl_binned", "u");

  int ny = LENGTH(y1), nb = LENGTH(b), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX5(ny, nb, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13;
  double *py1 = REAL(y1), *py2 = REAL(y2), *pb = REAL(b), *pc = REAL(c), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if(pb[i % nb] != -1) {
      e1 = pb[i % nb] + 1;
      e2 = R_pow(pl[i % nl], e1);
      e3 = R_pow(pu[i % nu], e1);
      e4 = R_pow(py1[i % ny], e1);
      e5 = R_pow(py2[i % ny], e1);
      e6 = e3 - e2;
      e7 = e5 - e4;
      e8 = log(pl[i % nl]);
      e9 = log(pu[i % nu]);
      e10 = log(py1[i % ny]);
      e11 = log(py2[i % ny]);
      e12 = sign(e6);
      e13 = sign(e7);

      pans[i] = pc[i % nc] * (((e5 * R_pow(e11, 2) - e4 * R_pow(e10, 2))/fabs(e7) - e13 * R_pow(e5 * e11 - e4 * e10, 2)/R_pow(e7, 2)) * e13 - ((e3 * R_pow(e9, 2) - e2 * R_pow(e8, 2))/fabs(e6) - e12 * R_pow(e3 * e9 - e2 * e8, 2)/R_pow(e6, 2)) * e12);

      pans[i] = -pans[i]; //return negative second derivative for 'gamlss2'

    } else {
      pans[i] = 0;
    }
  }

  UNPROTECT(1);
  return(ans);
}


// second derivative of ISD log-likelihood for counts of individual sizes
// computed using Deriv::Deriv(Deriv::Deriv(<loglik_bpl> for = -1 and != -1, respectively))
SEXP hess_bpl(SEXP y, SEXP b, SEXP c, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y, "hess_bpl", "y");
  CHECK_NUMERIC_INPUT(b, "hess_bpl", "b");
  CHECK_NUMERIC_INPUT(c, "hess_bpl", "c");
  CHECK_NUMERIC_INPUT(l, "hess_bpl", "l");
  CHECK_NUMERIC_INPUT(u, "hess_bpl", "u");

  int ny = LENGTH(y), nb = LENGTH(b), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX5(ny, nb, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double e1, e2, e3, e4, e5, e6, e7, e8, e9;
  double *py = REAL(y), *pc = REAL(c), *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if(pb[i % nb] != -1) {
      e1 = 1 + pb[i % nb];
      e2 = R_pow(pl[i % nl], e1);
      e3 = R_pow(pu[i % nu], e1);
      e4 = log(pl[i % nl]);
      e5 = log(pu[i % nu]);
      e6 = e2 * e4;
      e7 = e3 - e2;
      e8 = e3 * e5;
      e9 = e8 - e6;
      pans[i] = -(pc[i % nc] * ((e1 * (e3 * R_pow(e5, 2) - (R_pow(e9, 2)/e7 + e2 * R_pow(e4, 2))) + e8 - e6)/e7 + (1 - e1 * e9/e7)/e1)/e1);

      pans[i] = -pans[i]; //return negative second derivative
    } else {
      pans[i] = 0.0;
    }
  }

  UNPROTECT(1);
  return(ans);
}




/*
// PDF of ISD; adapted from 'isdbayes' package
SEXP d_bpl(SEXP y, SEXP b, SEXP l, SEXP u)
{
  int ny = LENGTH(y), nb = LENGTH(b), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nb, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if((py[i % ny] < pl[i % nl]) || (py[i % ny] > pu[i % nu])) {
      pans[i] = 0.0;
      continue;
    }
    if(pb[i % nb] != -1) {
      pans[i] = (pb[i % nb] + 1) * R_pow(py[i % ny], pb[i % nb] + 1) / ( R_pow(pu[i % nu], pb[i % nb] + 1) - R_pow(pl[i % nl], pb[i % nb] + 1) );
    } else {
      pans[i] = R_pow(py[i % ny], -2)/(pl[i % nl] * log(pu[i % nu]/pl[i % nl]));
    }
  }

  UNPROTECT(1);
  return(ans);
}
*/


/*
 dbounded_powerlaw_log <- function(x, b, xmin, xmax) {
 logdens <- rep(-Inf, length(x))
 valid <- (x >= xmin) & (x <= xmax)

 if (!any(valid)) return(logdens)

 xv <- x[valid]

 if (abs(b + 1) > .Machine$double.eps) {
 logC <- log(abs(b + 1)) -
 log(abs(xmax^(b + 1) - xmin^(b + 1)))

 logdens[valid] <- logC + b * log(xv)
 } else {
# b == -1
 logdens[valid] <- -log(xv) -
 log(log(xmax) - log(xmin))
 }

 return(logdens)
 }

 */

SEXP d_bpl_log(SEXP y, SEXP b, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y, "d_bpl", "y");
  CHECK_NUMERIC_INPUT(b, "d_bpl", "b");
  CHECK_NUMERIC_INPUT(l, "d_bpl", "l");
  CHECK_NUMERIC_INPUT(u, "d_bpl", "u");

  int ny = LENGTH(y), nb = LENGTH(b), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nb, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
  double logC;

  for(int i = 0; i < N; i++) {
    if((py[i % ny] < pl[i % nl]) || (py[i % ny] > pu[i % nu])) {
      pans[i] = R_NegInf;
      continue;
    }
    if(pb[i % nb] != -1) {
      logC = log(fabs(pb[i % nb] + 1)) - log( fabs( R_pow(pu[i % nu], pb[i % nb] + 1) - R_pow(pl[i % nl], pb[i % nb] + 1) ) );
      pans[i] = logC + pb[i % nb] * log(py[i % ny]);
    } else {
      pans[i] = -log(py[i % ny]) - log( log(pu[i % nu]) - log(pl[i % nl]) );
    }
  }

  UNPROTECT(1);
  return(ans);
}


// adapted from sizeSpectra::dPLB
SEXP d_bpl(SEXP y, SEXP b, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y, "d_bpl", "y");
  CHECK_NUMERIC_INPUT(b, "d_bpl", "b");
  CHECK_NUMERIC_INPUT(l, "d_bpl", "l");
  CHECK_NUMERIC_INPUT(u, "d_bpl", "u");

  //if( (!Rf_isLogical(use_log)) || (LOGICAL(use_log)[0] == NA_LOGICAL) ) Rf_error("argument 'log' is not of type \"logical\" in 'd_bpl'");
  //int uselog = LOGICAL(use_log)[0];

  int ny = LENGTH(y), nb = LENGTH(b), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nb, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
  double C;

  for(int i = 0; i < N; i++) {
    if((py[i % ny] < pl[i % nl]) || (py[i % ny] > pu[i % nu])) {
      pans[i] = 0.0;

      //if(uselog) pans[i] = R_NegInf;

      continue;
    }
    if(pb[i % nb] != -1) {
      C = (pb[i % nb] + 1) / ( R_pow(pu[i % nu], pb[i % nb] + 1) - R_pow(pl[i % nl], pb[i % nb] + 1) );
    } else {
      C = 1 / (log(pu[i % nu]) - log(pl[i % nl]));
    }

    pans[i] = C * R_pow(py[i % ny], pb[i % nb]);

    //if(uselog) pans[i] = log(pans[i]);
  }

  UNPROTECT(1);
  return(ans);
}






/*
 * initialize for individual sizes data;
 * adapted from https://github.com/andrew-edwards/sizeSpectra/blob/517c18d84f4326b59807de5235ab4cddac74876b/R/fitting.R#L1507
 */
SEXP init_bpl(SEXP y, SEXP c, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y, "init_bpl", "y");
  CHECK_NUMERIC_INPUT(c, "init_bpl", "c");
  CHECK_NUMERIC_INPUT(l, "init_bpl", "l");
  CHECK_NUMERIC_INPUT(u, "init_bpl", "u");

  int ny = LENGTH(y), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pc = REAL(c), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
  double sumc = 0.0, sumclogy = 0.0, l1 = R_PosInf;

  for(int i = 0; i < N; i++) {
    sumc += pc[i % nc];
    sumclogy += pc[i % nc] * log(py[i % ny]);
    l1 = (pl[i % nl] < l1) ? pl[i % nl] : l1;
  }

  for(int i = 0; i < N; i++) {
    //pans[i] = 1 / (log(pl[i % nl]) - sumclogy / sumc) - 1;
    pans[i] = 1 / (log(l1) - sumclogy / sumc) - 1;
  }

  UNPROTECT(1);
  return(ans);
}


/*
 * initialize for binned data
 *
 * #https://github.com/andrew-edwards/sizeSpectra/blob/517c18d84f4326b59807de5235ab4cddac74876b/R/simulating.R#L176

   sumCntLogMids = sum(binCounts * log(binMids))
   #MLEmid (maximum likelihood using midpoints) calculations.
   #Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007) as a starting point for nlm for MLE of b for PLB model.
   PL.bMLE = 1/( log(min(binBreaks)) - sumCntLogMids/sum(binCounts) ) - 1
 */

SEXP init_bpl_binned(SEXP y1, SEXP y2, SEXP c, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y1, "init_bpl_binned", "y1");
  CHECK_NUMERIC_INPUT(y2, "init_bpl_binned", "y2");
  CHECK_NUMERIC_INPUT(c, "init_bpl_binned", "c");
  CHECK_NUMERIC_INPUT(l, "init_bpl_binned", "l");
  CHECK_NUMERIC_INPUT(u, "init_bpl_binned", "u");

  int ny = LENGTH(y1), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py1 = REAL(y1), *py2 = REAL(y2), *pc = REAL(c), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
  double sumc = 0.0, sumclogmids = 0.0, l1 = R_PosInf;

  for(int i = 0; i < N; i++) {
    sumc += pc[i % nc];
    sumclogmids += pc[i % nc] * log( (py1[i % ny] + py2[i % ny]) / 2 );
    l1 = (pl[i % nl] < l1) ? pl[i % nl] : l1;
  }

  for(int i = 0; i < N; i++) {
    //pans[i] = 1 / ( log(pl[i % nl]) - sumclogmids / sumc ) - 1;
    pans[i] = 1 / ( log(l1) - sumclogmids / sumc ) - 1;
  }

  UNPROTECT(1);
  return(ans);
}



// adapted from 'sizeSpectra::rPLB'
// load 'rand' directly from R to avoid messing with random number generation in C
SEXP r_bpl(SEXP n, SEXP b, SEXP l, SEXP u, SEXP rand)
{
  CHECK_NUMERIC_INPUT(n, "r_bpl", "n");
  CHECK_NUMERIC_INPUT(b, "r_bpl", "b");
  CHECK_NUMERIC_INPUT(l, "r_bpl", "l");
  CHECK_NUMERIC_INPUT(u, "r_bpl", "u");
  CHECK_NUMERIC_INPUT(rand, "r_bpl", "rand");

  int N = INTEGER(n)[0];
  int nb = LENGTH(b), nl = LENGTH(l), nu = LENGTH(u), nrand = LENGTH(rand);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *prand = REAL(rand), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if((pl[i % nl] <= 0) || (pl[i % nl] >= pu[i % nu])) {
      pans[i] = NA_REAL;
      continue;
    }
    if(pb[i % nb] != -1) {
      pans[i] = R_pow((prand[i % nrand] * R_pow(pu[i % nu], (pb[i % nb] + 1)) + (1 - prand[i % nrand]) * R_pow(pl[i % nl], (pb[i % nb] + 1))), (1/(pb[i % nb] + 1)));
    } else {
      pans[i] = R_pow(pu[i % nu], prand[i % nrand]) * R_pow(pl[i % nl], (1 - prand[i % nrand]));
    }
  }

  UNPROTECT(1);
  return(ans);
}

// adapted from 'sizeSpectra::qPLB'
SEXP q_bpl(SEXP p, SEXP b, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(p, "q_bpl", "p");
  CHECK_NUMERIC_INPUT(b, "q_bpl", "b");
  CHECK_NUMERIC_INPUT(l, "q_bpl", "l");
  CHECK_NUMERIC_INPUT(u, "q_bpl", "u");

  int np = LENGTH(p), nb = LENGTH(b), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(np, nb, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *pp = REAL(p), *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if((pl[i % nl] <= 0) || (pl[i % nl] >= pu[i % nu]) || (pp[i % np] < 0) || (pp[i % np] > 1)) {
      pans[i] = NA_REAL;
      continue;
    }
    if(pb[i % nb] != -1) {
      pans[i] = R_pow((pp[i % np] * R_pow(pu[i % nu], (pb[i % nb] + 1)) + (1 - pp[i % np]) * R_pow(pl[i % nl], (pb[i % nb] + 1))), (1/(pb[i % nb] + 1)));
    } else {
      pans[i] = R_pow(pu[i % nu], pp[i % np]) * R_pow(pl[i % nl], (1 - pp[i % np]));
    }
  }

  UNPROTECT(1);
  return(ans);
}

// adapted from 'sizeSpectra::pPLB'
SEXP p_bpl(SEXP y, SEXP b, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y, "p_bpl", "y");
  CHECK_NUMERIC_INPUT(b, "p_bpl", "b");
  CHECK_NUMERIC_INPUT(l, "p_bpl", "l");
  CHECK_NUMERIC_INPUT(u, "p_bpl", "u");

  int ny = LENGTH(y), nb = LENGTH(b), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nb, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
  double xmintobplus1, denom, logxmin;

  for(int i = 0; i < N; i++) {
    if((pl[i % nl] <= 0) || (pl[i % nl] >= pu[i % nu])) {
      pans[i] = NA_REAL;
      continue;
    }

    if(py[i % ny] < pl[i % nl]) {
      pans[i] = 0;
      continue;
    }

    if(py[i % ny] > pu[i % nu]) {
      pans[i] = 1;
      continue;
    }

    if(pb[i % nb] != -1) {
      xmintobplus1 = R_pow(pl[i % nl], (pb[i % nb] + 1));
      denom = R_pow(pu[i % nu], (pb[i % nb] + 1)) - xmintobplus1;
      pans[i] = (R_pow(py[i % ny], pb[i % nb] + 1) - xmintobplus1) / denom;
    } else {
      logxmin = log(pl[i % nl]);
      denom = log(pu[i % nu]) - logxmin;
      pans[i] = (log(py[i % ny]) - logxmin) / denom;
    }
  }

  UNPROTECT(1);
  return(ans);
}

/*
 * adapted from
 *  https://www.int-res.com/articles/suppl/m636p019_supp1.pdf
 */
SEXP p_bpl_binned(SEXP y1, SEXP y2, SEXP b, SEXP l, SEXP u)
{
  CHECK_NUMERIC_INPUT(y1, "p_bpl_binned", "y1");
  CHECK_NUMERIC_INPUT(y2, "p_bpl_binned", "y2");
  CHECK_NUMERIC_INPUT(b, "p_bpl_binned", "b");
  CHECK_NUMERIC_INPUT(l, "p_bpl_binned", "l");
  CHECK_NUMERIC_INPUT(u, "p_bpl_binned", "u");

  int ny = LENGTH(y1), nb = LENGTH(b), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nb, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py1 = REAL(y1), *py2 = REAL(y2), *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
  double xmintobplus1, denom, logxmin;

  for(int i = 0; i < N; i++) {
    if((pl[i % nl] <= 0) || (pl[i % nl] >= pu[i % nu])) {
      pans[i] = NA_REAL;
      continue;
    }

    if(py1[i % ny] < pl[i % nl]) {
      pans[i] = 0;
      continue;
    }

    if(py2[i % ny] > pu[i % nu]) {
      pans[i] = 1;
      continue;
    }

    if(pb[i % nb] != -1) {
      pans[i] = ( R_pow(py2[i % ny], pb[i % nb] + 1) - R_pow(py1[i % ny], pb[i % nb] + 1) ) / ( R_pow(pu[i % nu], pb[i % nb] + 1) - R_pow(pl[i % nl], pb[i % nb] + 1) );
    } else {
      pans[i] = ( log(py2[i % ny]) - log(py1[i % ny]) ) / ( log(pu[i % nu]) - log(pl[i % nl]) );
    }
  }

  UNPROTECT(1);
  return(ans);
}





