# linop.primme 0.1.0

* First release. Vendors 'PRIMME' 3.2.3 and links it against R's own 'BLAS'
  and 'LAPACK', supplying `zhegv`, which R does not ship, in `src/lapack_ext.c`.
* `primme_eigsh()` computes a few eigenpairs of a Hermitian operator and
  `primme_svds()` a few singular triplets of a rectangular one, each given as a
  function applied to a block of columns. Both default to the GD+k method,
  whose results are reproducible run to run, and accept a convergence test
  written in R.
* `primme_engine()` reports the vendored 'PRIMME' version and the version of the
  calling contract with 'linop', whose `method = "primme"` uses this package.
