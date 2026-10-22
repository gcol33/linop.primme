# Story: linop.primme

## Claim
The PRIMME eigensolver as the compiled engine behind `method = "primme"` in linop: a few eigenpairs or singular triplets of an operator given by its action.

## Audience
linop users who need eigenpairs or singular triplets of large operators and want an implicitly restarted engine. They meet the package through linop's `eigs()` and `svds()`, not through its own calls.

## Relation
- to: linop (consumer); PRIMME (vendored library)
- phrase: the compiled engine behind `method = "primme"` in linop
- same: linop's certificate; values are Rayleigh quotients measured on the operator
- differs: implicit restarting with a fixed small basis, against linop's reference Lanczos
- home: linop's spectral article

## Arcs

### A1 Hermitian eigenpairs and singular triplets through linop
- question: can linop's `eigs()` and `svds()` use a production engine?
- answer: installed, `method = "auto"` takes PRIMME; `method = "primme"` by name
- evidence: linop tests/testthat/test-backends.R (closed-form spectra, both engines)
- home: linop vignettes/spectral.Rmd
- also: README

### A2 The same certificate
- question: does switching engines change what I can trust?
- answer: linop's convergence test decides every pair (PRIMME's maxConvTol shortcut compiled out, src/Makevars)
- home: linop vignettes/certificates.Rmd
- also: README

## Non-goals
An interface of its own beyond the two calls linop makes.

## Documents
- DESCRIPTION: claim, PRIMME citations, vendoring
- README: claim, relation, installation with build tools, first example through linop
- help pages: `primme_eigsh()`, `primme_svds()`, `primme_engine()`

## Open
- No R-CMD-check workflow and no pkgdown site yet.
