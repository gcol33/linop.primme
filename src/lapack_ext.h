#ifndef LP_LAPACK_EXT_H
#define LP_LAPACK_EXT_H

#include <complex.h>
#include <stddef.h>
#include <Rconfig.h>

typedef double complex lp_complex;

/* zhegv for itype = 1, with the argument types PRIMME's blaslapack_private.h
   declares for double complex and Fortran's hidden character lengths last. */
void lp_zhegv(int *itype, const char *jobz, const char *uplo, int *n,
              lp_complex *a, int *lda, lp_complex *b, int *ldb, double *w,
              lp_complex *work, int *lwork, double *rwork, int *info,
              FC_LEN_T jobz_len, FC_LEN_T uplo_len);

#endif
