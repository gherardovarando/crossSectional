context("Utility functions")

test_that("rHurwitz produce stable matrices", {
  expect_true(all(replicate(100, {
    B <- rHurwitz(p = sample(1:50, size = 1))
    all(Re(eigen(B)$values) < 0)
  })))
})

test_that("rStableMetzler produce stable matrices", {
  expect_true(all(replicate(100, {
    B <- rStableMetzler(p = sample(1:50, size = 1), d = runif(1))
    all(Re(eigen(B)$values) < 0)
  })))
})

test_that("minus log likelihood is minimized by empirical covariance matrix", {
  B <- rHurwitz(p = 10)
  M <- gmat::chol_mh(N = 1000, p = 10)
  S <- rOU(n = 1000, B = B)
  mlls <- apply(M, MARGIN = 3, FUN = function(A){
    mll(P = A, S = var(S$dat) )
  })
  expect_true(all(mlls > mll(P = solve(var(S$data)), S = var(S$data))))
})

