#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>


// macros to find max; used to find 'max(length(a), length(b), ...)' to emulate recycling behaviour
#define MAX2(x1, x2) (((x1) > (x2)) ? (x1) : (x2))
#define MAX3(x1, x2, x3) MAX2(MAX2((x1), (x2)), (x3))
#define MAX4(x1, x2, x3, x4) MAX2(MAX3((x1), (x2), (x3)), (x4))
#define MAX5(x1, x2, x3, x4, x5) MAX2(MAX4((x1), (x2), (x3), (x4)), (x5))



// adapted from Edwards et al. 2017 (sizeSpectra::negLL.PLB.counts)
SEXP loglik_bpl(SEXP y, SEXP b, SEXP c, SEXP l, SEXP u)
{
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

// adapted from sizeSpectra::dPLB
SEXP d_bpl(SEXP y, SEXP b, SEXP l, SEXP u)
{
  int ny = LENGTH(y), nb = LENGTH(b), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nb, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pb = REAL(b), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
  double C;

  for(int i = 0; i < N; i++) {
    if((py[i % ny] < pl[i % nl]) || (py[i % ny] > pu[i % nu])) {
      pans[i] = 0.0;
      continue;
    }
    if(pb[i % nb] != -1) {
      C = (pb[i % nb] + 1) / ( R_pow(pu[i % nu], pb[i % nb] + 1) - R_pow(pl[i % nl], pb[i % nb] + 1) );pans[i] = (pb[i % nb] + 1) * R_pow(py[i % ny], pb[i % nb] + 1) / ( R_pow(pu[i % nu], pb[i % nb] + 1) - R_pow(pl[i % nl], pb[i % nb] + 1) );
    } else {
      C = 1 / (log(pu[i % nu]) - log(pl[i % nl]));
    }

    pans[i] = C * R_pow(py[i % ny], pb[i % nb]);
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
  int ny = LENGTH(y), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pc = REAL(c), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
  double sumc = 0.0, sumclogy = 0.0;

  for(int i = 0; i < N; i++) {
    sumc += pc[i % nc];
    sumclogy += pc[i % nc] * log(py[i % ny]);
  }

  for(int i = 0; i < N; i++) {
    pans[i] = 1 / (log(pl[i % nl]) - sumclogy / sumc) - 1;
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
  int ny = LENGTH(y1), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py1 = REAL(y1), *py2 = REAL(y2), *pc = REAL(c), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
  double sumc = 0.0, sumclogmids = 0.0;

  for(int i = 0; i < N; i++) {
    sumc += pc[i % nc];
    sumclogmids += pc[i % nc] * log( (py1[i % ny] + py2[i % ny]) / 2 );
  }

  for(int i = 0; i < N; i++) {
    pans[i] = 1 / ( log(pl[i % nl]) - sumclogmids / sumc ) - 1;
  }

  UNPROTECT(1);
  return(ans);
}



// adapted from 'sizeSpectra::rPLB'
// load 'rand' directly from R to avoid messing with random number generation in C
SEXP r_bpl(SEXP n, SEXP b, SEXP l, SEXP u, SEXP rand)
{
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





