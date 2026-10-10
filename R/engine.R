## The two calls linop makes, and the contract it checks before making them.
##
## Everything here is in PRIMME's own terms: a function applying the operator to
## a block, sizes, a target, a tolerance and the engine's knobs. Validation of
## the operator, the certificate and the choice between engines are linop's.

PRIMME_VERSION <- "3.2.3"
CONTRACT_VERSION <- "1.1"

EIGSH_TARGET <- c(LA = 0L, SA = 1L, LM = 2L, SM = 3L)
SVDS_TARGET <- c(LM = 0L, SM = 1L)

#' Engine identity
#'
#' The version of the vendored 'PRIMME' library and of the calling contract
#' between this package and 'linop'.
#'
#' @return A list with `engine_version` and `contract_version`.
#' @examples
#' primme_engine()
#' @export
primme_engine <- function() {
  list(engine_version = PRIMME_VERSION, contract_version = CONTRACT_VERSION)
}

## The apply as the C side calls it. A real run refuses a complex result rather
## than dropping its imaginary part, and a complex run accepts a real one, which
## is the same vector.
block_apply <- function(f, complex) {
  force(f)
  function(...) {
    Y <- f(...)
    if (!is.matrix(Y)) Y <- as.matrix(Y)
    if (complex) {
      if (!is.complex(Y)) storage.mode(Y) <- "complex"
    } else {
      if (is.complex(Y)) stop("the apply returned a complex block for a real operator")
      if (!is.double(Y)) storage.mode(Y) <- "double"
    }
    Y
  }
}

check_count <- function(x, what) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != round(x)) {
    stop(sprintf("%s must be a positive whole number", what), call. = FALSE)
  }
  as.integer(x)
}

check_convtest <- function(convtest) {
  if (!is.null(convtest) && !is.function(convtest)) {
    stop("convtest must be a function or NULL", call. = FALSE)
  }
}

opt_count <- function(x, what) if (is.null(x)) NULL else check_count(x, what)

opt_flag <- function(x, what) {
  if (is.null(x)) return(NULL)
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s must be TRUE, FALSE or NULL", what), call. = FALSE)
  }
  as.integer(x)
}

#' A few eigenpairs of a Hermitian operator
#'
#' Calls `dprimme()` or `zprimme()` with the default method, on an operator
#' given by its action on a block of columns.
#'
#' @param apply A function taking an `n x b` matrix and returning the operator
#'   applied to it.
#' @param n The order of the operator.
#' @param k How many eigenpairs.
#' @param which `"LA"`, `"SA"` (largest, smallest algebraic), `"LM"` or `"SM"`
#'   (largest, smallest magnitude).
#' @param tol PRIMME's `eps`. A pair is converged when its residual norm is at
#'   most `tol` times the estimate of `||A||` and at most `tol` times the
#'   magnitude of its value, each plus the arithmetic floor
#'   `floor * .Machine$double.eps * (||A|| + |value|)`, the second waived where
#'   that magnitude is below `floor * .Machine$double.eps` times the estimate of
#'   `||A||`. These are the backward-error and residual bounds of 'linop''s
#'   eigenpair certificate.
#'   PRIMME refuses a `tol` below machine epsilon.
#' @param maxit PRIMME's `maxMatvecs`, the budget of single-column applies.
#' @param ncv PRIMME's `maxBasisSize`, or `NULL` for its default.
#' @param min_restart PRIMME's `minRestartSize`, or `NULL` for its default.
#' @param block PRIMME's `maxBlockSize`, or `NULL` for its default.
#' @param locking PRIMME's `locking`, `TRUE` or `FALSE`, or `NULL` for its
#'   default.
#' @param v0 An initial guess of length `n`, or `NULL`.
#' @param complex Whether the operator is complex.
#' @param seed Seed for PRIMME's own random initial vectors.
#' @param norm An estimate of `||A||` for the convergence test, or `0` for
#'   PRIMME to estimate it.
#' @param floor `c` in the arithmetic floor of the convergence test.
#' @param convtest `NULL` for the convergence test described under `tol`, or a
#'   function `(value, x, rnorm, anorm)` returning `TRUE` when the pair is
#'   converged. It receives PRIMME's Ritz value, the unit Ritz vector or `NULL`
#'   where PRIMME passes none, PRIMME's residual norm on the operator `apply`
#'   applies, and the estimate of that operator's norm. A caller iterating on a
#'   transformed operator, such as a shift-and-invert, measures the pair on the
#'   operator it wants here.
#' @param preset The name of a PRIMME preset method, such as
#'   `"PRIMME_DEFAULT_MIN_MATVECS"` (GD+k), `"PRIMME_DEFAULT_MIN_TIME"`
#'   (JDQMR), `"PRIMME_Arnoldi"` or `"PRIMME_DYNAMIC"`. The default is GD+k
#'   rather than PRIMME's own default, `"PRIMME_DYNAMIC"`, because that one
#'   switches between the two by measured timings, so two runs of the same call
#'   need not return the same vectors.
#' @return A list with `values` (PRIMME's Ritz values), `vectors`, `resNorms`,
#'   `nconv`, `info` (PRIMME's return code, `0` on success and `-3` when the
#'   budget ran out), `stats` and the `params` PRIMME settled on.
#' @examples
#' n <- 100
#' lap <- function(X) 2 * X - rbind(X[-1, , drop = FALSE], 0) -
#'   rbind(0, X[-n, , drop = FALSE])
#' fit <- primme_eigsh(lap, n, k = 2, which = "SA")
#' fit$values - 4 * sin(1:2 * pi / (2 * (n + 1)))^2
#' @export
primme_eigsh <- function(apply, n, k, which = "LA", tol = 1e-10, maxit = 10000L,
                         ncv = NULL, min_restart = NULL, block = NULL,
                         locking = NULL, v0 = NULL, complex = FALSE, seed = 1L,
                         norm = 0, floor = 4, preset = "PRIMME_DEFAULT_MIN_MATVECS",
                         convtest = NULL) {
  if (!is.function(apply)) stop("apply must be a function", call. = FALSE)
  n <- check_count(n, "n"); k <- check_count(k, "k")
  if (k > n) stop("k cannot exceed n", call. = FALSE)
  check_convtest(convtest)
  if (!which %in% names(EIGSH_TARGET)) {
    stop("which is one of LA, SA, LM, SM", call. = FALSE)
  }
  if (!is.null(v0)) {
    v0 <- as.vector(v0)
    if (length(v0) != n) stop("v0 must have length n", call. = FALSE)
    storage.mode(v0) <- if (complex) "complex" else "double"
  }
  .Call(lp_eigsh, block_apply(apply, complex), environment(),
        as.double(n), k, EIGSH_TARGET[[which]], as.double(tol),
        opt_count(ncv, "ncv"), opt_count(min_restart, "min_restart"),
        opt_count(block, "block"), opt_flag(locking, "locking"),
        as.double(check_count(maxit, "maxit")),
        v0, isTRUE(complex), as.integer(seed), as.double(norm), as.double(floor),
        as.character(preset), convtest)
}

#' A few singular triplets of an operator
#'
#' Calls `dprimme_svds()` or `zprimme_svds()` with the default method, on an
#' operator given by its action and its adjoint.
#'
#' @param apply A function `(X, trans)`: `A %*% X` when `trans` is `FALSE`, and
#'   the conjugate transpose of `A` applied to `X` when it is `TRUE`.
#' @param m,n The operator is `m x n`.
#' @param k How many triplets.
#' @param which `"LM"` (largest) or `"SM"` (smallest).
#' @param preset The PRIMME preset method for both stages of the singular value
#'   solver; see [primme_eigsh()].
#' @param convtest `NULL` for PRIMME's own convergence test, or a function
#'   `(value, u, v, rnorm, anorm)` returning `TRUE` when the triplet is
#'   converged. It receives PRIMME's singular value, the left and right
#'   vectors or `NULL` for a side PRIMME does not pass (its first stage, on the
#'   normal equations, passes one side), PRIMME's residual norm and its
#'   estimate of `||A||`.
#' @inheritParams primme_eigsh
#' @return A list with `d`, `u`, `v`, `resNorms`, `nconv`, `info`, the overall
#'   `stats` and those of each of PRIMME's two stages, and `params`. PRIMME
#'   returns only the triplets it converged, so after an exhausted budget
#'   (`info` `-103` or `-203`, the eigensolver's `-3` from the first or second
#'   stage) `d`, `u` and `v` hold `nconv` of them, possibly none.
#' @examples
#' A <- diag(c(5, 4, 3, 2, 1))[, 1:4]
#' fit <- primme_svds(function(X, trans) if (trans) crossprod(A, X) else A %*% X,
#'                    5, 4, k = 2)
#' fit$d
#' @export
primme_svds <- function(apply, m, n, k, which = "LM", tol = 1e-10,
                        maxit = 10000L, ncv = NULL, block = NULL, locking = NULL,
                        complex = FALSE, seed = 1L, norm = 0,
                        preset = "PRIMME_DEFAULT_MIN_MATVECS", convtest = NULL) {
  if (!is.function(apply)) stop("apply must be a function", call. = FALSE)
  m <- check_count(m, "m"); n <- check_count(n, "n"); k <- check_count(k, "k")
  if (k > min(m, n)) stop("k cannot exceed min(m, n)", call. = FALSE)
  if (!which %in% names(SVDS_TARGET)) stop("which is LM or SM", call. = FALSE)
  check_convtest(convtest)
  .Call(lp_svds, block_apply(apply, complex), environment(),
        as.double(m), as.double(n), k, SVDS_TARGET[[which]], as.double(tol),
        opt_count(ncv, "ncv"), opt_count(block, "block"),
        opt_flag(locking, "locking"),
        as.double(check_count(maxit, "maxit")), isTRUE(complex),
        as.integer(seed), as.double(norm), as.character(preset), convtest)
}
