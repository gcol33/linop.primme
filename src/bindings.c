/*
 * .Call bindings from R to PRIMME's Hermitian eigensolver and its singular
 * value solver.
 *
 * The operator is an R function. PRIMME hands the callback a block of columns
 * with a leading dimension; the callback copies the block into an R matrix,
 * evaluates the function on it, and copies the result back at the leading
 * dimension PRIMME asked for. The copy is O(n) per column applied, against the
 * O(n * basis size) PRIMME spends orthogonalising that column.
 *
 * The R function is evaluated with R_tryEvalSilent, so an error or an
 * interrupt inside it never unwinds through PRIMME's frames. The callback sets
 * PRIMME's error flag instead, PRIMME returns, its memory is released, and the
 * condition is raised in R afterwards with the message the function produced.
 */

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "primme.h"

typedef struct {
  SEXP fn;          /* the operator: function(X) or function(X, trans) */
  SEXP env;
  int cplx;         /* 1 when the operator is complex */
  PRIMME_INT rows;  /* eigs: n. svds: m, the length of A x */
  PRIMME_INT cols;  /* eigs: n. svds: n, the length of A^H y */
  int failed;       /* 1 when the R function errored or returned the wrong shape */
  double floor;     /* c in the arithmetic floor c * eps * ||A|| of the convergence test */
  char msg[1024];
} lp_ctx;

static void lp_fail(lp_ctx *c, const char *msg)
{
  c->failed = 1;
  strncpy(c->msg, msg, sizeof(c->msg) - 1);
  c->msg[sizeof(c->msg) - 1] = '\0';
}

/* Apply the R function to nin x bs columns of x at leading dimension ldx and
   write the nout x bs result into y at leading dimension ldy. trans < 0 calls
   fn(X); trans = 0 or 1 calls fn(X, trans). Returns 0 on success. */
static int lp_apply(lp_ctx *c, const void *x, PRIMME_INT ldx, PRIMME_INT nin,
                    void *y, PRIMME_INT ldy, PRIMME_INT nout, int bs, int trans)
{
  if (c->failed) return 1;
  const SEXPTYPE type = c->cplx ? CPLXSXP : REALSXP;
  const size_t esz = c->cplx ? sizeof(Rcomplex) : sizeof(double);

  SEXP X = PROTECT(Rf_allocMatrix(type, (int)nin, bs));
  char *xd = c->cplx ? (char *)COMPLEX(X) : (char *)REAL(X);
  for (int j = 0; j < bs; j++) {
    memcpy(xd + (size_t)j * nin * esz, (const char *)x + (size_t)j * ldx * esz,
           (size_t)nin * esz);
  }

  SEXP call = PROTECT(trans < 0 ? Rf_lang2(c->fn, X)
                                : Rf_lang3(c->fn, X, Rf_ScalarLogical(trans)));
  int err = 0;
  SEXP Y = R_tryEvalSilent(call, c->env, &err);
  if (err) {
    const char *m = R_curErrorBuf();
    lp_fail(c, (m && *m) ? m : "the operator's apply was interrupted");
    UNPROTECT(2);
    return 1;
  }
  PROTECT(Y);
  if (TYPEOF(Y) != type) {
    lp_fail(c, c->cplx ? "the apply returned a non-complex block for a complex operator"
                       : "the apply returned a non-double block for a real operator");
    UNPROTECT(3);
    return 1;
  }
  if ((R_xlen_t)XLENGTH(Y) != (R_xlen_t)nout * bs) {
    lp_fail(c, "the apply returned a block of the wrong size");
    UNPROTECT(3);
    return 1;
  }
  const char *yd = c->cplx ? (const char *)COMPLEX(Y) : (const char *)REAL(Y);
  for (int j = 0; j < bs; j++) {
    memcpy((char *)y + (size_t)j * ldy * esz, yd + (size_t)j * nout * esz,
           (size_t)nout * esz);
  }
  UNPROTECT(3);
  return 0;
}

static void lp_eigs_matvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
                           int *bs, struct primme_params *pp, int *ierr)
{
  lp_ctx *c = (lp_ctx *)pp->matrix;
  *ierr = lp_apply(c, x, *ldx, c->rows, y, *ldy, c->rows, *bs, -1);
}

/* PRIMME's transpose flag is 0 for A x and 1 for A^H y. */
static void lp_svds_matvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
                           int *bs, int *transpose,
                           struct primme_svds_params *ps, int *ierr)
{
  lp_ctx *c = (lp_ctx *)ps->matrix;
  if (*transpose) {
    *ierr = lp_apply(c, x, *ldx, c->rows, y, *ldy, c->cols, *bs, 1);
  } else {
    *ierr = lp_apply(c, x, *ldx, c->cols, y, *ldy, c->rows, *bs, 0);
  }
}

/* A pair is converged when its residual meets both of the bounds linop's
   eigenpair certificate checks, each with the certificate's arithmetic floor
   f = c eps (||A|| + |value|), Higham's bound on the computed residual: the
   backward error ||r|| <= tol ||A|| + f, and the relative residual
   ||r|| <= tol |value| + f. The second is waived where |value| has itself
   fallen to c eps ||A||, below which that ratio stops measuring the solve.
   ||A|| is the caller's estimate when given, else PRIMME's largest |Ritz value|
   seen, which is a lower bound. */
static void lp_conv_test(double *eval, void *evec, double *rNorm, int *isconv,
                         struct primme_params *pp, int *ierr)
{
  (void)evec;
  const lp_ctx *c = (const lp_ctx *)pp->convtest;
  const double anorm = pp->aNorm > 0 ? pp->aNorm : pp->stats.estimateLargestSVal;
  const double tol = pp->eps;
  const double v = eval != NULL ? fabs(*eval) : 0.0;
  const double f = c->floor * DBL_EPSILON * (anorm + v);
  int ok = *rNorm <= tol * anorm + f;
  if (ok && eval != NULL && v > c->floor * DBL_EPSILON * anorm) {
    ok = *rNorm <= tol * v + f;
  }
  *isconv = ok;
  *ierr = 0;
}

static int lp_int(SEXP s, const char *what)
{
  if (!Rf_isNumeric(s) || XLENGTH(s) != 1) Rf_error("%s must be a single number", what);
  return Rf_asInteger(s);
}

static void lp_seed(PRIMME_INT *iseed, int seed)
{
  /* LAPACK's larnv requires iseed(4) odd and every entry in [0, 4095]. */
  unsigned int s = (unsigned int)seed;
  iseed[0] = (PRIMME_INT)(s % 4096u);
  iseed[1] = (PRIMME_INT)((s / 4096u) % 4096u);
  iseed[2] = (PRIMME_INT)((s / 16777216u) % 4096u);
  iseed[3] = (PRIMME_INT)(((s * 2u) + 1u) % 4096u) | 1;
}

static SEXP lp_stats(const primme_stats *st)
{
  const char *nm[] = {"numMatvecs", "numOuterIterations", "numRestarts",
                      "numOrthoInnerProds", "flopsDense", "elapsedTime",
                      "timeMatvec", "timeOrtho", "timeDense",
                      "estimateLargestSVal", "estimateMinEVal",
                      "estimateMaxEVal", "maxConvTol", "lockingIssue", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, nm));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal((double)st->numMatvecs));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal((double)st->numOuterIterations));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal((double)st->numRestarts));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(st->numOrthoInnerProds));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(st->flopsDense));
  SET_VECTOR_ELT(out, 5, Rf_ScalarReal(st->elapsedTime));
  SET_VECTOR_ELT(out, 6, Rf_ScalarReal(st->timeMatvec));
  SET_VECTOR_ELT(out, 7, Rf_ScalarReal(st->timeOrtho));
  SET_VECTOR_ELT(out, 8, Rf_ScalarReal(st->timeDense));
  SET_VECTOR_ELT(out, 9, Rf_ScalarReal(st->estimateLargestSVal));
  SET_VECTOR_ELT(out, 10, Rf_ScalarReal(st->estimateMinEVal));
  SET_VECTOR_ELT(out, 11, Rf_ScalarReal(st->estimateMaxEVal));
  SET_VECTOR_ELT(out, 12, Rf_ScalarReal(st->maxConvTol));
  SET_VECTOR_ELT(out, 13, Rf_ScalarReal((double)st->lockingIssue));
  UNPROTECT(1);
  return out;
}

static SEXP lp_params(int maxBasisSize, int minRestartSize, int maxBlockSize,
                      int locking, int orth_explicit)
{
  const char *nm[] = {"maxBasisSize", "minRestartSize", "maxBlockSize",
                      "locking", "orthExplicit", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, nm));
  SET_VECTOR_ELT(out, 0, Rf_ScalarInteger(maxBasisSize));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(minRestartSize));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(maxBlockSize));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(locking));
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(orth_explicit));
  UNPROTECT(1);
  return out;
}

/* target codes, shared with R/engine.R:
   0 largest algebraic, 1 smallest algebraic,
   2 largest magnitude, 3 smallest magnitude. */
SEXP lp_eigsh(SEXP fn, SEXP env, SEXP s_n, SEXP s_k, SEXP s_target,
              SEXP s_tol, SEXP s_maxBasis, SEXP s_minRestart, SEXP s_block,
              SEXP s_locking, SEXP s_maxMatvecs, SEXP s_v0, SEXP s_cplx, SEXP s_seed,
              SEXP s_aNorm, SEXP s_floor, SEXP s_preset)
{
  const PRIMME_INT n = (PRIMME_INT)Rf_asReal(s_n);
  const int k = lp_int(s_k, "k");
  const int target = lp_int(s_target, "target");
  const int cplx = Rf_asLogical(s_cplx);
  if (n < 1 || k < 1 || k > n) Rf_error("need 1 <= k <= n");

  lp_ctx c;
  memset(&c, 0, sizeof c);
  c.fn = fn; c.env = env; c.cplx = cplx; c.rows = n; c.cols = n;
  c.floor = Rf_asReal(s_floor);

  primme_params pp;
  primme_initialize(&pp);
  pp.n = n;
  pp.numEvals = k;
  pp.matrixMatvec = lp_eigs_matvec;
  pp.matrix = &c;
  pp.convTestFun = lp_conv_test;
  pp.convTestFun_type = primme_op_double;
  pp.convtest = &c;
  pp.eps = Rf_asReal(s_tol);
  pp.printLevel = 0;
  pp.outputFile = NULL;
  pp.maxMatvecs = (PRIMME_INT)Rf_asReal(s_maxMatvecs);
  pp.aNorm = Rf_asReal(s_aNorm);
  lp_seed(pp.iseed, Rf_asInteger(s_seed));

  double shift0 = 0.0;
  switch (target) {
  case 0: pp.target = primme_largest; break;
  case 1: pp.target = primme_smallest; break;
  case 2: pp.target = primme_largest_abs; break;
  case 3: pp.target = primme_closest_abs; break;
  default: Rf_error("unknown target code %d", target);
  }
  if (target >= 2) { pp.numTargetShifts = 1; pp.targetShifts = &shift0; }

  if (s_maxBasis != R_NilValue) pp.maxBasisSize = lp_int(s_maxBasis, "maxBasisSize");
  if (s_minRestart != R_NilValue) pp.minRestartSize = lp_int(s_minRestart, "minRestartSize");
  if (s_block != R_NilValue) pp.maxBlockSize = lp_int(s_block, "maxBlockSize");
  if (s_locking != R_NilValue) pp.locking = Rf_asInteger(s_locking);

  int preset = 0;
  if (primme_constant_info(CHAR(STRING_ELT(s_preset, 0)), &preset) != 0) {
    primme_free(&pp);
    Rf_error("unknown PRIMME preset '%s'", CHAR(STRING_ELT(s_preset, 0)));
  }
  const int pm = primme_set_method((primme_preset_method)preset, &pp);
  if (pm != 0) { primme_free(&pp); Rf_error("primme_set_method failed with %d", pm); }

  SEXP evals = PROTECT(Rf_allocVector(REALSXP, k));
  SEXP rnorms = PROTECT(Rf_allocVector(REALSXP, k));
  SEXP evecs = PROTECT(Rf_allocMatrix(cplx ? CPLXSXP : REALSXP, (int)n, k));
  if (cplx) memset(COMPLEX(evecs), 0, sizeof(Rcomplex) * (size_t)n * k);
  else memset(REAL(evecs), 0, sizeof(double) * (size_t)n * k);

  if (s_v0 != R_NilValue) {
    if (XLENGTH(s_v0) != n) Rf_error("v0 must have length n");
    if (cplx) memcpy(COMPLEX(evecs), COMPLEX(s_v0), sizeof(Rcomplex) * (size_t)n);
    else memcpy(REAL(evecs), REAL(s_v0), sizeof(double) * (size_t)n);
    pp.initSize = 1;
  }

  int ret;
  if (cplx) {
    ret = zprimme(REAL(evals), (PRIMME_COMPLEX_DOUBLE *)COMPLEX(evecs),
                  REAL(rnorms), &pp);
  } else {
    ret = dprimme(REAL(evals), REAL(evecs), REAL(rnorms), &pp);
  }

  const int nconv = pp.initSize;
  SEXP stats = PROTECT(lp_stats(&pp.stats));
  SEXP params = PROTECT(lp_params(pp.maxBasisSize, pp.minRestartSize,
                                  pp.maxBlockSize, pp.locking,
                                  pp.orth == primme_orth_explicit_I));
  primme_free(&pp);

  if (c.failed) Rf_error("%s", c.msg);

  const char *nm[] = {"values", "vectors", "resNorms", "nconv", "info",
                      "stats", "params", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, nm));
  SET_VECTOR_ELT(out, 0, evals);
  SET_VECTOR_ELT(out, 1, evecs);
  SET_VECTOR_ELT(out, 2, rnorms);
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(nconv));
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(ret));
  SET_VECTOR_ELT(out, 5, stats);
  SET_VECTOR_ELT(out, 6, params);
  UNPROTECT(6);
  return out;
}

/* target codes: 0 largest, 1 smallest. */
SEXP lp_svds(SEXP fn, SEXP env, SEXP s_m, SEXP s_n, SEXP s_k, SEXP s_target,
             SEXP s_tol, SEXP s_maxBasis, SEXP s_block, SEXP s_locking,
             SEXP s_maxMatvecs,
             SEXP s_cplx, SEXP s_seed, SEXP s_aNorm, SEXP s_preset)
{
  const PRIMME_INT m = (PRIMME_INT)Rf_asReal(s_m);
  const PRIMME_INT n = (PRIMME_INT)Rf_asReal(s_n);
  const int k = lp_int(s_k, "k");
  const int target = lp_int(s_target, "target");
  const int cplx = Rf_asLogical(s_cplx);
  const PRIMME_INT mn = m < n ? m : n;
  if (m < 1 || n < 1 || k < 1 || k > mn) Rf_error("need 1 <= k <= min(m, n)");

  lp_ctx c;
  memset(&c, 0, sizeof c);
  c.fn = fn; c.env = env; c.cplx = cplx; c.rows = m; c.cols = n;

  primme_svds_params ps;
  primme_svds_initialize(&ps);
  ps.m = m;
  ps.n = n;
  ps.numSvals = k;
  ps.matrixMatvec = lp_svds_matvec;
  ps.matrix = &c;
  ps.eps = Rf_asReal(s_tol);
  ps.printLevel = 0;
  ps.outputFile = NULL;
  ps.maxMatvecs = (PRIMME_INT)Rf_asReal(s_maxMatvecs);
  ps.aNorm = Rf_asReal(s_aNorm);
  lp_seed(ps.iseed, Rf_asInteger(s_seed));
  switch (target) {
  case 0: ps.target = primme_svds_largest; break;
  case 1: ps.target = primme_svds_smallest; break;
  default: Rf_error("unknown target code %d", target);
  }
  if (s_maxBasis != R_NilValue) ps.maxBasisSize = lp_int(s_maxBasis, "maxBasisSize");
  if (s_block != R_NilValue) ps.maxBlockSize = lp_int(s_block, "maxBlockSize");
  if (s_locking != R_NilValue) ps.locking = Rf_asInteger(s_locking);

  int preset = 0;
  if (primme_constant_info(CHAR(STRING_ELT(s_preset, 0)), &preset) != 0) {
    primme_svds_free(&ps);
    Rf_error("unknown PRIMME preset '%s'", CHAR(STRING_ELT(s_preset, 0)));
  }
  const int pm = primme_svds_set_method(primme_svds_default,
                                        (primme_preset_method)preset,
                                        (primme_preset_method)preset, &ps);
  if (pm != 0) { primme_svds_free(&ps); Rf_error("primme_svds_set_method failed with %d", pm); }

  SEXP svals = PROTECT(Rf_allocVector(REALSXP, k));
  SEXP rnorms = PROTECT(Rf_allocVector(REALSXP, k));
  const size_t len = (size_t)(m + n) * k;
  SEXP svecs = PROTECT(Rf_allocVector(cplx ? CPLXSXP : REALSXP, (R_xlen_t)len));
  if (cplx) memset(COMPLEX(svecs), 0, sizeof(Rcomplex) * len);
  else memset(REAL(svecs), 0, sizeof(double) * len);

  int ret;
  if (cplx) {
    ret = zprimme_svds(REAL(svals), (PRIMME_COMPLEX_DOUBLE *)COMPLEX(svecs),
                       REAL(rnorms), &ps);
  } else {
    ret = dprimme_svds(REAL(svals), REAL(svecs), REAL(rnorms), &ps);
  }

  const int nconv = ps.initSize;
  SEXP stats = PROTECT(lp_stats(&ps.primme.stats));
  SEXP stats2 = PROTECT(lp_stats(&ps.primmeStage2.stats));
  const char *snm[] = {"numMatvecs", "elapsedTime", ""};
  SEXP sst = PROTECT(Rf_mkNamed(VECSXP, snm));
  SET_VECTOR_ELT(sst, 0, Rf_ScalarReal((double)ps.stats.numMatvecs));
  SET_VECTOR_ELT(sst, 1, Rf_ScalarReal(ps.stats.elapsedTime));
  SEXP params = PROTECT(lp_params(ps.maxBasisSize, ps.primme.minRestartSize,
                                  ps.maxBlockSize, ps.locking,
                                  ps.primme.orth == primme_orth_explicit_I));
  primme_svds_free(&ps);

  if (c.failed) Rf_error("%s", c.msg);

  /* svecs holds the nconv left vectors (m each) followed by the nconv right
     vectors: PRIMME packs only what converged, so after an exhausted budget
     the right block starts at m * nconv, not at m * k. */
  SEXP U = PROTECT(Rf_allocMatrix(cplx ? CPLXSXP : REALSXP, (int)m, nconv));
  SEXP V = PROTECT(Rf_allocMatrix(cplx ? CPLXSXP : REALSXP, (int)n, nconv));
  const size_t esz = cplx ? sizeof(Rcomplex) : sizeof(double);
  char *sd = cplx ? (char *)COMPLEX(svecs) : (char *)REAL(svecs);
  if (nconv > 0) {
    memcpy(cplx ? (char *)COMPLEX(U) : (char *)REAL(U), sd, esz * (size_t)m * nconv);
    memcpy(cplx ? (char *)COMPLEX(V) : (char *)REAL(V), sd + esz * (size_t)m * nconv,
           esz * (size_t)n * nconv);
  }
  SEXP d = PROTECT(Rf_allocVector(REALSXP, nconv));
  SEXP rn = PROTECT(Rf_allocVector(REALSXP, nconv));
  for (int j = 0; j < nconv; j++) {
    REAL(d)[j] = REAL(svals)[j];
    REAL(rn)[j] = REAL(rnorms)[j];
  }

  const char *nm[] = {"d", "u", "v", "resNorms", "nconv", "info", "stats",
                      "stats_stage1", "stats_stage2", "params", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, nm));
  SET_VECTOR_ELT(out, 0, d);
  SET_VECTOR_ELT(out, 1, U);
  SET_VECTOR_ELT(out, 2, V);
  SET_VECTOR_ELT(out, 3, rn);
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(nconv));
  SET_VECTOR_ELT(out, 5, Rf_ScalarInteger(ret));
  SET_VECTOR_ELT(out, 6, sst);
  SET_VECTOR_ELT(out, 7, stats);
  SET_VECTOR_ELT(out, 8, stats2);
  SET_VECTOR_ELT(out, 9, params);
  UNPROTECT(12);
  return out;
}

void lp_zhegv(const int *itype, const char *jobz, const char *uplo,
              const int *n, Rcomplex *a, const int *lda, Rcomplex *b,
              const int *ldb, double *w, Rcomplex *work, const int *lwork,
              double *rwork, int *info);

/* Direct access to lp_zhegv, so its contract can be tested on its own. */
SEXP lp_zhegv_call(SEXP A, SEXP B, SEXP s_uplo)
{
  const int n = Rf_nrows(A);
  const char *uplo = CHAR(STRING_ELT(s_uplo, 0));
  SEXP a = PROTECT(Rf_duplicate(A));
  SEXP b = PROTECT(Rf_duplicate(B));
  SEXP w = PROTECT(Rf_allocVector(REALSXP, n));
  int one = 1, lwork = -1, info = 0;
  Rcomplex q;
  double *rwork = (double *)R_alloc(3 * (size_t)(n > 0 ? n : 1), sizeof(double));
  lp_zhegv(&one, "V", uplo, &n, COMPLEX(a), &n, COMPLEX(b), &n, REAL(w), &q,
           &lwork, rwork, &info);
  lwork = (int)q.r;
  if (lwork < 1) lwork = 1;
  Rcomplex *work = (Rcomplex *)R_alloc((size_t)lwork, sizeof(Rcomplex));
  lp_zhegv(&one, "V", uplo, &n, COMPLEX(a), &n, COMPLEX(b), &n, REAL(w), work,
           &lwork, rwork, &info);
  const char *nm[] = {"values", "vectors", "info", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, nm));
  SET_VECTOR_ELT(out, 0, w);
  SET_VECTOR_ELT(out, 1, a);
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(info));
  UNPROTECT(4);
  return out;
}
