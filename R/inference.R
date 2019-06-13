#' Estimate Ornstein-Uhlenbeck process
#'
#' Estimate the coefficients and the noise term of an Ornstein-Uhlenbeck
#' process (\eqn{dX(t) = BX(t) + \sqrt{C} dW(t)}) using penalized maximum-likelihood.
#'
#' @param Sigma the empirical covariance matrix
#' @param B an initial guess for the coefficient matrix B
#' @param C the initial guess for the noise matrix C
#' @param eps stopping criteria
#' @param alpha step size
#' @param maxIter maximum number of iterations
#' @param trace if >0 print info
#' @param lambda penalization coefficient
#' @param r logical
#'
#' @return the estimated B matrix (\code{estimateBLL}) or
#' the estiamted C matrix (\code{estiamteCLL}).
#' @importFrom lyapunov clyap clyap2
#' @export
estimateBLL <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                        alpha = 0.1,
                        maxIter = 1000, trace = 0,
                        lambda = 0, r = FALSE){
  p <- ncol(Sigma)

  a <- Inf
  n <- 0
  E <- matrix(nrow = p, ncol = p, 0)
  #ix <- (1: p^2 )[B != 0]
  ix <- 1:(p^2)
  while (abs(a) > eps && n < maxIter  || n < 2){
    if (r) ix <- (1: p^2 )[B != 0]
    u <- rep(0, length(ix))
    n <- n + 1

    allres <- clyap(A = B, Q = C, all = TRUE)
    S <- matrix(nrow = p, data = allres[[6]])
    AA <- matrix(nrow = p, data = allres[[4]])
    EE <- matrix(nrow = p, data = allres[[5]])
    WKV <- allres[[7]]
    P <- solve(S)

    tmp <- P %*% Sigma %*% P - P

    u <- vapply(1:length(ix), function(i){
      E[ix[i]] <- 1
      Cp <- E %*% S + S %*% t(E)
      D <- clyap2(A = AA, Q = Cp, E = EE,  WKV = WKV)
      E[ix[i]] <- 0
      sum(tmp * D )
    }, FUN.VALUE = 1)

    #Bold <- B

    B[ix] <- B[ix] + alpha * u

    ### soft thres
    B[ix] <- sign(B[ix]) * (abs(B[ix]) - alpha * lambda)
    B[abs(B) < (alpha * lambda)] <- 0

    #a <- sqrt(sum( (B - Bold)^2 ))
    a <- sqrt(sum(u ^ 2))
    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||diff||=", signif(a))
  }

  return(B)
}

#' @rdname estimateBLL
#' @param C0 penalization matrix
#' @importFrom lyapunov clyap clyap2
#' @export
estimateCLL <- function(Sigma, B, C = diag(ncol(Sigma)), C0 = diag(ncol(Sigma)),
                        eps =  1e-2,
                        alpha = 0.1,
                        maxIter = 1000, trace = 0,
                        lambda = 0, r = FALSE){
  p <- ncol(Sigma)

  a <- Inf
  n <- 0
  E <- matrix(nrow = p, ncol = p, 0)
  ix <- (1: p^2 )[B != 0]
  ixc <- (1: p^2 )[C != 0]
  allres <- clyap(A = B, Q = C, all = TRUE)

  AA <- matrix(nrow = p, data = allres[[4]])
  EE <- matrix(nrow = p, data = allres[[5]])
  WKV <- allres[[7]]
  while (abs(a) > eps && n < maxIter  || n < 2){
    if (r) ix <- (1: p^2 )[B != 0]
    u <- rep(0, length(ix))
    n <- n + 1

    S <- clyap2(A = AA, Q = C, E = EE,  WKV = WKV)
    P <- solve(S)

    tmp <- P %*% Sigma %*% P - P

    v <- vapply(1:length(ixc), function(i){
      E[ixc[i]] <- 1
      Cp <- E + t(E)
      D <- clyap2(A = AA, Q = Cp, E = EE,  WKV = WKV)
      E[ixc[i]] <- 0
      sum(tmp * D )
    }, FUN.VALUE = 1) - 2 * lambda * ( (C - C0)[ixc])

    C[ixc] <- C[ixc] + alpha * v

    a <- sqrt(sum( v ^ 2))

    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||diff||=", signif(a))
  }

  return(C)
}


