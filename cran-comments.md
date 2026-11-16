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

## Test environments

<!-- filled from the win-builder and R-hub results before submission -->
