lap_apply <- function(n) {
  function(X) 2 * X - rbind(X[-1, , drop = FALSE], 0) - rbind(0, X[-n, , drop = FALSE])
}
lap_values <- function(n) 4 * sin(seq_len(n) * pi / (2 * (n + 1)))^2

random_unitary <- function(m, seed) {
  set.seed(seed)
  qr.Q(qr(matrix(rnorm(m * m) + 1i * rnorm(m * m), m)))
}

test_that("primme_engine reports the vendored version and the contract", {
  e <- primme_engine()
  expect_identical(e$engine_version, "3.2.3")
  expect_identical(e$contract_version, "1.1")
})

test_that("real hermitian: every target reaches the closed-form spectrum", {
  n <- 200
  ex <- lap_values(n)
  want <- list(SA = ex[1:4], LA = rev(ex)[1:4], SM = ex[1:4], LM = rev(ex)[1:4])
  for (w in names(want)) {
    f <- primme_eigsh(lap_apply(n), n, 4, which = w, tol = 1e-12)
    expect_identical(f$info, 0L)
    expect_identical(f$nconv, 4L)
    expect_lt(max(abs(sort(f$values) - sort(want[[w]]))), 1e-12)
    R <- lap_apply(n)(f$vectors) - f$vectors %*% diag(f$values)
    expect_lt(max(sqrt(colSums(R^2))), 1e-10)
  }
})

test_that("complex hermitian, with and without the block path that calls zhegv", {
  m <- 60
  Q <- random_unitary(m, 2)
  lam <- seq(-3, 5, length.out = m)
  H <- Q %*% diag(lam) %*% Conj(t(Q))
  for (b in c(1L, 3L)) {
    f <- primme_eigsh(function(X) H %*% X, m, 3, which = "LA", complex = TRUE,
                      block = b, tol = 1e-12)
    ## maxBlockSize > 1 makes PRIMME orthogonalise explicitly and solve the
    ## projected problem against V^H V, which is the generalized driver.
    expect_identical(f$params$orthExplicit, b > 1L)
    expect_identical(f$info, 0L)
    expect_lt(max(abs(f$values - rev(lam)[1:3])), 1e-12)
    expect_true(is.complex(f$vectors))
  }
})

test_that("lp_zhegv solves A x = lambda B x for both triangles", {
  set.seed(7)
  for (n in c(1L, 5L, 23L)) {
    A <- matrix(rnorm(n * n) + 1i * rnorm(n * n), n); A <- A + Conj(t(A))
    G <- matrix(rnorm(n * n) + 1i * rnorm(n * n), n)
    B <- Conj(t(G)) %*% G + diag(n)
    ref <- sort(Re(eigen(solve(B, A), only.values = TRUE)$values))
    for (uplo in c("U", "L")) {
      ## Only the named triangle may be read: poison the other one.
      Ap <- A; Bp <- B
      if (n > 1) {
        off <- if (uplo == "U") lower.tri(A) else upper.tri(A)
        Ap[off] <- NaN; Bp[off] <- NaN
      }
      out <- .Call(linop.primme:::lp_zhegv_call, Ap, Bp, uplo)
      expect_identical(out$info, 0L)
      expect_lt(max(abs(out$values - ref)), 1e-10 * max(1, abs(ref)))
      X <- out$vectors
      expect_lt(max(Mod(A %*% X - B %*% X %*% diag(out$values, n))), 1e-9)
      ## B-orthonormal, as LAPACK's zhegv returns them for itype = 1.
      expect_lt(max(Mod(Conj(t(X)) %*% B %*% X - diag(n))), 1e-10)
    }
  }
})

test_that("lp_zhegv reports a B that is not positive definite as n + i", {
  n <- 4L
  A <- diag(as.complex(1:4))
  B <- diag(as.complex(c(1, 1, -1, 1)))
  out <- .Call(linop.primme:::lp_zhegv_call, A, B, "U")
  expect_identical(out$info, n + 3L)
})

test_that("singular triplets, real and complex, largest and smallest", {
  set.seed(3)
  A <- matrix(rnorm(50 * 30), 50)
  s <- svd(A)$d
  ap <- function(X, trans) if (trans) crossprod(A, X) else A %*% X
  f <- primme_svds(ap, 50, 30, 3, tol = 1e-12)
  expect_identical(f$info, 0L)
  expect_lt(max(abs(f$d - s[1:3])), 1e-12)
  expect_lt(max(abs(A %*% f$v - f$u %*% diag(f$d))), 1e-10)

  Ac <- matrix(rnorm(40 * 25) + 1i * rnorm(40 * 25), 40)
  sc <- svd(Ac)$d
  apc <- function(X, trans) if (trans) Conj(t(Ac)) %*% X else Ac %*% X
  f <- primme_svds(apc, 40, 25, 3, which = "SM", complex = TRUE, tol = 1e-12)
  expect_identical(f$info, 0L)
  expect_lt(max(abs(sort(f$d) - sort(sc)[1:3])), 1e-11)
})

test_that("an error inside the apply comes back as an R error", {
  expect_error(primme_eigsh(function(X) stop("boom"), 50, 2), "boom")
  expect_error(primme_svds(function(X, trans) stop("bang"), 10, 8, 2), "bang")
})

test_that("a real run refuses a complex apply rather than dropping its imaginary part", {
  expect_error(primme_eigsh(function(X) X + 1i * X, 30, 2), "complex block")
})

test_that("a wrong-sized apply is refused", {
  expect_error(primme_eigsh(function(X) X[-1, , drop = FALSE], 30, 2), "wrong size")
})

test_that("an exhausted budget is a return code, not an error", {
  n <- 200
  f <- primme_eigsh(lap_apply(n), n, 4, which = "SA", tol = 1e-14, maxit = 30)
  expect_identical(f$info, -3L)
  ## PRIMME tests the budget after the block that crosses it, so it can be
  ## exceeded by at most one block.
  expect_lte(f$stats$numMatvecs, 30 + f$params$maxBlockSize)
})

test_that("the default preset returns the same vectors on every run", {
  ## PRIMME's own default switches methods on measured timings, so it does not;
  ## GD+k, the default here, does.
  n <- 400
  set.seed(1)
  v0 <- rnorm(n)
  runs <- lapply(1:3, function(i)
    primme_eigsh(lap_apply(n), n, 3, which = "SA", tol = 1e-10, block = 3L, v0 = v0))
  expect_identical(runs[[2]]$vectors, runs[[1]]$vectors)
  expect_identical(runs[[3]]$vectors, runs[[1]]$vectors)
  A <- matrix(rnorm(60 * 40), 60)
  ap <- function(X, trans) if (trans) crossprod(A, X) else A %*% X
  s <- lapply(1:3, function(i) primme_svds(ap, 60, 40, 3, block = 3L))
  expect_identical(s[[2]]$u, s[[1]]$u)
  expect_identical(s[[3]]$v, s[[1]]$v)
})

test_that("the caller's convergence test decides, measured on another operator", {
  ## Shift-and-invert: PRIMME iterates on (L - I)^-1, and its residual there bounds
  ## the residual on L only by ||L - I|| times it. The test measures on L instead.
  n <- 60
  L <- 2 * diag(n); L[cbind(1:(n - 1), 2:n)] <- -1; L[cbind(2:n, 1:(n - 1))] <- -1
  S <- solve(L - diag(n))
  tol <- 1e-10
  seen <- 0L
  conv <- function(value, x, rnorm, anorm) {
    if (is.null(x)) return(FALSE)
    seen <<- seen + 1L
    Lx <- L %*% x
    mu <- sum(x * Lx) / sum(x^2)
    sqrt(sum((Lx - mu * x)^2)) / sqrt(sum(x^2)) <= tol * abs(mu)
  }
  f <- primme_eigsh(function(X) S %*% X, n, 3, which = "LM", tol = tol, block = 3L,
                    convtest = conv)
  expect_identical(f$info, 0L)
  expect_gt(seen, 0L)
  mu <- colSums(f$vectors * (L %*% f$vectors)) / colSums(f$vectors^2)
  rel <- sqrt(colSums((L %*% f$vectors - f$vectors %*% diag(mu))^2)) /
    (abs(mu) * sqrt(colSums(f$vectors^2)))
  expect_lte(max(rel), tol)
})

test_that("a convergence test that never accepts runs out the budget", {
  n <- 100
  f <- primme_eigsh(lap_apply(n), n, 2, which = "SA", maxit = 60,
                    convtest = function(value, x, rnorm, anorm) FALSE)
  expect_identical(f$info, -3L)
})

test_that("a convergence test that errors or answers otherwise is refused", {
  n <- 30
  expect_error(primme_eigsh(lap_apply(n), n, 2,
                            convtest = function(value, x, rnorm, anorm) stop("no verdict")),
               "no verdict")
  expect_error(primme_eigsh(lap_apply(n), n, 2,
                            convtest = function(value, x, rnorm, anorm) NA),
               "TRUE or FALSE")
  expect_error(primme_eigsh(lap_apply(n), n, 2, convtest = "yes"), "function or NULL")
})

test_that("the caller's convergence test decides for singular triplets", {
  set.seed(3)
  A <- qr.Q(qr(matrix(rnorm(40 * 30), 40))) %*% diag(seq(10, 0.1, length.out = 30)) %*%
    t(qr.Q(qr(matrix(rnorm(30 * 30), 30))))
  ap <- function(X, trans) if (trans) crossprod(A, X) else A %*% X
  tol <- 1e-10
  sides <- character(0)
  conv <- function(value, u, v, rnorm, anorm) {
    if (is.null(u) && is.null(v)) return(FALSE)
    sides <<- union(sides, c(if (!is.null(u)) "u", if (!is.null(v)) "v"))
    if (is.null(u)) u <- A %*% v
    if (is.null(v)) v <- crossprod(A, u)
    u <- u / sqrt(sum(u^2)); v <- v / sqrt(sum(v^2))
    d <- sum(u * (A %*% v))
    r <- sqrt(sum((A %*% v - d * u)^2) + sum((crossprod(A, u) - d * v)^2)) / sqrt(2)
    r <= tol * abs(d)
  }
  f <- primme_svds(ap, 40, 30, 3, which = "SM", tol = tol, block = 3L, convtest = conv)
  expect_identical(f$info, 0L)
  expect_true(all(c("u", "v") %in% sides))
  U <- f$u; V <- f$v
  d <- colSums(U * (A %*% V))
  r <- sqrt(colSums((A %*% V - U %*% diag(d))^2) + colSums((crossprod(A, U) - V %*% diag(d))^2)) /
    sqrt(2)
  expect_lte(max(r / abs(d)), tol)

  never <- primme_svds(ap, 40, 30, 2, maxit = 80, convtest = function(...) FALSE)
  expect_true(never$info %in% c(-103L, -203L))
})
