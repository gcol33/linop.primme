#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP lp_eigsh(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
              SEXP, SEXP, SEXP, SEXP);
SEXP lp_svds(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
             SEXP, SEXP, SEXP);
SEXP lp_zhegv_call(SEXP, SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
  {"lp_eigsh", (DL_FUNC) &lp_eigsh, 14},
  {"lp_svds", (DL_FUNC) &lp_svds, 13},
  {"lp_zhegv_call", (DL_FUNC) &lp_zhegv_call, 3},
  {NULL, NULL, 0}
};

void R_init_linop_primme(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
