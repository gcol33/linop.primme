## The two calls linop makes, and the contract it checks before making them.
##
## Everything here is in PRIMME's own terms: a function applying the operator to
## a block, sizes, a target, a tolerance and the engine's knobs. Validation of
## the operator, the certificate and the choice between engines are linop's.

PRIMME_VERSION <- "3.2.3"
CONTRACT_VERSION <- "1.0"

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

opt_count <- function(x, what) if (is.null(x)) NULL else check_count(x, what)

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
#' @param tol PRIMME's `eps`: a pair is converged when its residual norm is at
#'   most `tol` times the estimate of `||A||`.
#' @param maxit PRIMME's `maxMatvecs`, the budget of single-column applies.
#' @param ncv PRIMME's `maxBasisSize`, or `NULL` for its default.
#' @param min_restart PRIMME's `minRestartSize`, or `NULL` for its default.
#' @param block PRIMME's `maxBlockSize`, or `NULL` for its default.
#' @param v0 An initial guess of length `n`, or `NULL`.
#' @param complex Whether the operator is complex.
#' @param seed Seed for PRIMME's own random initial vectors.
#' @param norm An estimate of `||A||` for the convergence test, or `0` for
#'   PRIMME to estimate it.
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
                         v0 = NULL, complex = FALSE, seed = 1L, norm = 0) {
  if (!is.function(apply)) stop("apply must be a function", call. = FALSE)
  n <- check_count(n, "n"); k <- check_count(k, "k")
  if (k > n) stop("k cannot exceed n", call. = FALSE)
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
        opt_count(block, "block"), as.double(check_count(maxit, "maxit")),
        v0, isTRUE(complex), as.integer(seed), as.double(norm))
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
#' @inheritParams primme_eigsh
#' @return A list with `d`, `u`, `v`, `resNorms`, `nconv`, `info`, the overall
#'   `stats` and those of each of PRIMME's two stages, and `params`.
#' @examples
#' A <- diag(c(5, 4, 3, 2, 1))[, 1:4]
#' fit <- primme_svds(function(X, trans) if (trans) crossprod(A, X) else A %*% X,
#'                    5, 4, k = 2)
#' fit$d
#' @export
primme_svds <- function(apply, m, n, k, which = "LM", tol = 1e-10,
                        maxit = 10000L, ncv = NULL, block = NULL,
                        complex = FALSE, seed = 1L, norm = 0) {
  if (!is.function(apply)) stop("apply must be a function", call. = FALSE)
  m <- check_count(m, "m"); n <- check_count(n, "n"); k <- check_count(k, "k")
  if (k > min(m, n)) stop("k cannot exceed min(m, n)", call. = FALSE)
  if (!which %in% names(SVDS_TARGET)) stop("which is LM or SM", call. = FALSE)
  .Call(lp_svds, block_apply(apply, complex), environment(),
        as.double(m), as.double(n), k, SVDS_TARGET[[which]], as.double(tol),
        opt_count(ncv, "ncv"), opt_count(block, "block"),
        as.double(check_count(maxit, "maxit")), isTRUE(complex),
        as.integer(seed), as.double(norm))
}
