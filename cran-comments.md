## New submission

linop.primme vendors the 'PRIMME' eigensolver library (BSD 3-clause,
College of William & Mary) in `src/primme/`. Its copyright holders and authors
are in `Authors@R`, and `inst/COPYRIGHTS` lists them together with the
local changes to the 'PRIMME' sources. 'PRIMME' is linked against R's own
'BLAS' and 'LAPACK'; `zhegv`, which R's 'LAPACK' does not include, is
implemented in `src/lapack_ext.c` from routines R does ship.

The package is the compiled engine for the 'linop' package (to be submitted
after this one), which will list it in Suggests.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* Possibly misspelled words in DESCRIPTION: 'PReconditioned', 'MultiMethod'
  and 'Eigensolver' are the library's name as its authors spell it
  (PReconditioned Iterative MultiMethod Eigensolver); 'Stathopoulos' and
  'McCombs' are authors of the cited papers; 'primme' is the value of the
  'linop' argument `method = "primme"`.

## Test environments

* win-builder, R-devel (2026-11-16 r90914 ucrt): 1 NOTE (above), check time 89 s
* win-builder, R 4.6.2 (2026-10-29 ucrt): 1 NOTE (above), check time 42 s
* local Windows 11, R 4.6.0: `R CMD check --as-cran`, 1 NOTE (new submission)
* GitHub Actions, R release on macos-latest (arm64), windows-latest and
  ubuntu-latest: OK
* R-hub, gcc-asan, clang-asan, clang-ubsan and valgrind: OK
