#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>


// macros to find max; used to find 'max(length(a), length(b), ...)' to emulate recycling behaviour
#define MAX2(x1, x2) (((x1) > (x2)) ? (x1) : (x2))
#define MAX3(x1, x2, x3) MAX2(MAX2((x1), (x2)), (x3))
#define MAX4(x1, x2, x3, x4) MAX2(MAX3((x1), (x2), (x3)), (x4))
#define MAX5(x1, x2, x3, x4, x5) MAX2(MAX4((x1), (x2), (x3), (x4)), (x5))



// PDF of ISD; adapted from 'isdbayes' package
SEXP d_bpl(SEXP y, SEXP mu, SEXP l, SEXP u)
{
  int ny = LENGTH(y), nmu = LENGTH(mu), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nmu, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pmu = REAL(mu), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if((py[i % ny] < pl[i % nl]) || (py[i % ny] > pu[i % nu])) {
      pans[i] = 0.0;
      continue;
    }
    if(pmu[i % nmu] != -1) {
      pans[i] = (pmu[i % nmu] + 1) * pow(py[i % ny], (pmu[i % nmu] + 1)) / ( pow(pu[i % nu], (pmu[i % nmu] + 1)) - pow(pl[i % nl], (pmu[i % nmu] + 1)) );
    } else {
      pans[i] = pow(py[i % ny], -2)/(pl[i % nl] * log(pu[i % nu]/pl[i % nl]));
    }
  }

  UNPROTECT(1);
  return(ans);
}



// first derivative of ISD log-likelihood
// computed using Deriv::Deriv(<log-likelihood>)
SEXP score_bpl(SEXP y, SEXP mu, SEXP c, SEXP l, SEXP u)
{
  int ny = LENGTH(y), nmu = LENGTH(mu), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX5(ny, nmu, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double e1, e2, e3;
  double *py = REAL(y), *pc = REAL(c), *pmu = REAL(mu), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if(pmu[i % nmu] != -1) {
      e1 = 1 + pmu[i % nmu];
      e2 = pow(pl[i % nl], e1);
      e3 = pow(pu[i % nu], e1);
      pans[i] = pc[i % nc] * ((1 - e1 * (e3 * log(pu[i % nu]) - e2 * log(pl[i % nl]))/(e3 - e2))/e1 + log(py[i % ny]));
    } else {
      pans[i] = pc[i % nc] * log(py[i % ny]);
    }
  }

  UNPROTECT(1);
  return(ans);
}


// second derivative of ISD log-likelihood
// computed using Deriv::Deriv(Deriv::Deriv(<log-likelihood>))
SEXP hess_bpl(SEXP y, SEXP mu, SEXP c, SEXP l, SEXP u)
{
  int ny = LENGTH(y), nmu = LENGTH(mu), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX5(ny, nmu, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double e1, e2, e3, e4, e5, e6, e7, e8, e9;
  double *py = REAL(y), *pc = REAL(c), *pmu = REAL(mu), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if(pmu[i % nmu] != -1) {
      e1 = 1 + pmu[i % nmu];
      e2 = pow(pl[i % nl], e1);
      e3 = pow(pu[i % nu], e1);
      e4 = log(pl[i % nl]);
      e5 = log(pu[i % nu]);
      e6 = e2 * e4;
      e7 = e3 - e2;
      e8 = e3 * e5;
      e9 = e8 - e6;
      pans[i] = -(pc[i % nc] * ((e1 * (e3 * pow(e5, 2) - (pow(e9, 2)/e7 + e2 * pow(e4, 2))) + e8 - e6)/e7 + (1 - e1 * e9/e7)/e1)/e1);

      pans[i] = -pans[i];
    } else {
      pans[i] = 0.0;
    }
  }

  UNPROTECT(1);
  return(ans);
}


// log-likelihood; adapted from Edwards et al. 2017 (from 'sizeSpectra' code)
SEXP loglik_bpl(SEXP y, SEXP mu, SEXP c, SEXP l, SEXP u)
{
  int ny = LENGTH(y), nmu = LENGTH(mu), nc = LENGTH(c), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX5(ny, nmu, nc, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pc = REAL(c), *pmu = REAL(mu), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if(pmu[i % nmu] != -1) {
      pans[i] = pc[i % nc] * log((pmu[i % nmu] + 1)/(pow(pu[i % nu], pmu[i % nmu] + 1) - pow(pl[i % nl], pmu[i % nmu] + 1))) + pmu[i % nmu] * pc[i % nc] * log(py[i % ny]);
    } else {
      pans[i] =  -pc[i % nc] * log(log(pu[i % nu]) - log(pl[i % nl])) - pc[i % nc] * log(py[i % ny]);
    }
  }

  UNPROTECT(1);
  return(ans);
}


/*
 * initialize;
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


// adapted from 'sizeSpectra'
SEXP r_bpl(SEXP n, SEXP mu, SEXP l, SEXP u, SEXP rand)
{
  int N = INTEGER(n)[0];
  int nmu = LENGTH(mu), nl = LENGTH(l), nu = LENGTH(u), nrand = LENGTH(rand);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *pmu = REAL(mu), *pl = REAL(l), *pu = REAL(u), *prand = REAL(rand), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if((pl[i % nl] <= 0) || (pl[i % nl] >= pu[i % nu])) {
      pans[i] = NA_REAL;
      continue;
    }
    if(pmu[i % nmu] != -1) {
      pans[i] = pow((prand[i % nrand] * pow(pu[i % nu], (pmu[i % nmu] + 1)) + (1 - prand[i % nrand]) * pow(pl[i % nl], (pmu[i % nmu] + 1))), (1/(pmu[i % nmu] + 1)));
    } else {
      pans[i] = pow(pu[i % nu], prand[i % nrand]) * pow(pl[i % nl], (1 - prand[i % nrand]));
    }
  }

  UNPROTECT(1);
  return(ans);
}

// adapted from 'sizeSpectra'
SEXP q_bpl(SEXP p, SEXP mu, SEXP l, SEXP u)
{
  int np = LENGTH(p), nmu = LENGTH(mu), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(np, nmu, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *pp = REAL(p), *pmu = REAL(mu), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);

  for(int i = 0; i < N; i++) {
    if((pl[i % nl] <= 0) || (pl[i % nl] >= pu[i % nu]) || (pp[i % np] < 0) || (pp[i % np] > 1)) {
      pans[i] = NA_REAL;
      continue;
    }
    if(pmu[i % nmu] != -1) {
      pans[i] = pow((pp[i % np] * pow(pu[i % nu], (pmu[i % nmu] + 1)) + (1 - pp[i % np]) * pow(pl[i % nl], (pmu[i % nmu] + 1))), (1/(pmu[i % nmu] + 1)));
    } else {
      pans[i] = pow(pu[i % nu], pp[i % np]) * pow(pl[i % nl], (1 - pp[i % np]));
    }
  }

  UNPROTECT(1);
  return(ans);
}

// adapted from 'sizeSpectra'
SEXP p_bpl(SEXP y, SEXP mu, SEXP l, SEXP u)
{
  int ny = LENGTH(y), nmu = LENGTH(mu), nl = LENGTH(l), nu = LENGTH(u);
  int N = MAX4(ny, nmu, nl, nu);

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, N));

  double *py = REAL(y), *pmu = REAL(mu), *pl = REAL(l), *pu = REAL(u), *pans = REAL(ans);
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

    if(pmu[i % nmu] != -1) {
      xmintobplus1 = pow(pl[i % nl], (pmu[i % nmu] + 1));
      denom = pow(pu[i % nu], (pmu[i % nmu] + 1)) - xmintobplus1;
      pans[i] = (pow(py[i % ny], pmu[i % nmu] + 1) - xmintobplus1) / denom;
    } else {
      logxmin = log(pl[i % nl]);
      denom = log(pu[i % nu]) - logxmin;
      pans[i] = (log(py[i % ny]) - logxmin) / denom;
    }
  }

  UNPROTECT(1);
  return(ans);
}




