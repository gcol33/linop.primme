/*
 * The complex Hermitian-definite generalized eigensolver PRIMME's projected
 * problem calls, built from routines R's bundled LAPACK exports.
 *
 * PRIMME is compiled with -DUSE_XHEEV and -DUSE_XHEGV, so its dense kernels
 * call the simple drivers zheev/dsyev and zhegv/dsygv rather than the expert
 * drivers zheevx/zhegvx. R ships zheev, dsyev and dsygv. It does not ship
 * zhegv, nor zhegst, which zhegv is built on, so this file supplies zhegv by
 * the reference LAPACK reduction, under the private name lp_zhegv that the
 * Makevars maps PRIMME's call onto:
 *
 *   B = U^H U                 (zpotrf; or B = L L^H when uplo = 'L')
 *   C = U^-H A U^-1           (two triangular solves; L^-1 A L^-H)
 *   C y = lambda y            (zheev)
 *   x = U^-1 y                (one triangular solve; L^-H y)
 *
 * Only itype = 1, the problem A x = lambda B x, is implemented, because that is
 * the only type PRIMME requests: both call sites in blaslapack.c pass the
 * constant 1. Any other itype is reported as an illegal first argument.
 *
 * The reduction forms C from the full Hermitian A with two solves, where the
 * reference zhegst updates one triangle with zhemm/zher2k. The two produce C up
 * to rounding and zheev reads only the triangle named by uplo.
 */

#ifdef zhegv_
#undef zhegv_
#endif

#include <complex.h>
#include <Rconfig.h>

/* Declared exactly as PRIMME's blaslapack_private.h declares them for double
   complex, with Fortran's hidden character lengths, so every declaration of a
   symbol in the package agrees. */
typedef double complex lp_complex;
void zpotrf_(const char *uplo, int *n, lp_complex *a, int *lda, int *info,
             FC_LEN_T);
void zheev_(const char *jobz, const char *uplo, int *n, lp_complex *a, int *lda,
            double *w, lp_complex *work, int *ldwork, double *rwork, int *info,
            FC_LEN_T, FC_LEN_T);
void ztrsm_(const char *side, const char *uplo, const char *transa,
            const char *diag, int *m, int *n, lp_complex *alpha, lp_complex *a,
            int *lda, lp_complex *b, int *ldb, FC_LEN_T, FC_LEN_T, FC_LEN_T,
            FC_LEN_T);

void lp_zhegv(int *itype, const char *jobz, const char *uplo, int *n,
              lp_complex *a, int *lda, lp_complex *b, int *ldb, double *w,
              lp_complex *work, int *lwork, double *rwork, int *info,
              FC_LEN_T jobz_len, FC_LEN_T uplo_len)
{
  const int N = *n;
  const int upper = (*uplo == 'U' || *uplo == 'u');
  const int wantz = (*jobz == 'V' || *jobz == 'v');
  lp_complex one = 1.0;
  (void)jobz_len;
  (void)uplo_len;
  int nn = *n, la = *lda, lb = *ldb;

  *info = 0;
  if (*itype != 1) *info = -1;
  else if (!wantz && !(*jobz == 'N' || *jobz == 'n')) *info = -2;
  else if (!upper && !(*uplo == 'L' || *uplo == 'l')) *info = -3;
  else if (N < 0) *info = -4;
  else if (*lda < (N > 1 ? N : 1)) *info = -6;
  else if (*ldb < (N > 1 ? N : 1)) *info = -8;
  if (*info != 0) return;

  /* Workspace query: the reduction itself needs none, so the requirement is
     zheev's. */
  if (*lwork == -1) {
    zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info, 1, 1);
    return;
  }
  if (N == 0) return;

  /* B = U^H U or L L^H. A failure at column i means B is not positive
     definite, which LAPACK reports as n + i. */
  int pinfo = 0;
  zpotrf_(uplo, n, b, ldb, &pinfo, 1);
  if (pinfo > 0) { *info = N + pinfo; return; }
  if (pinfo < 0) { *info = pinfo; return; }

  /* The full Hermitian A from the triangle that holds it. */
  const int LDA = *lda;
  for (int j = 0; j < N; j++) {
    a[j + (size_t)j * LDA] = creal(a[j + (size_t)j * LDA]);
    for (int i = j + 1; i < N; i++) {
      if (upper) a[i + (size_t)j * LDA] = conj(a[j + (size_t)i * LDA]);
      else a[j + (size_t)i * LDA] = conj(a[i + (size_t)j * LDA]);
    }
  }

  /* C = U^-H A U^-1, or L^-1 A L^-H. */
  if (upper) {
    ztrsm_("L", "U", "C", "N", &nn, &nn, &one, b, &lb, a, &la, 1, 1, 1, 1);
    ztrsm_("R", "U", "N", "N", &nn, &nn, &one, b, &lb, a, &la, 1, 1, 1, 1);
  } else {
    ztrsm_("L", "L", "N", "N", &nn, &nn, &one, b, &lb, a, &la, 1, 1, 1, 1);
    ztrsm_("R", "L", "C", "N", &nn, &nn, &one, b, &lb, a, &la, 1, 1, 1, 1);
  }

  zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info, 1, 1);
  if (*info != 0 || !wantz) return;

  /* x = U^-1 y, or L^-H y. */
  if (upper) {
    ztrsm_("L", "U", "N", "N", &nn, &nn, &one, b, &lb, a, &la, 1, 1, 1, 1);
  } else {
    ztrsm_("L", "L", "C", "N", &nn, &nn, &one, b, &lb, a, &la, 1, 1, 1, 1);
  }
}
