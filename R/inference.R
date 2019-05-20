#' Estimate B minimizing Frobenius norm
#'
#' @param Sigma the empirical covariance matrix
#' @param B an initial guess for matrix B
#' @param C the value of the noise matrix C
#' @param eps stopping criteria
#' @param alpha step length
#' @param maxIter
#' @param trace
#' @param lambda
#'
#' @return the estimated B matrix
#' @importFrom lyapunov clyap
#' @export
estimateBS <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                      alpha = 0.1, maxIter = 1000, trace = 0,
                      lambda = 0.1, r = TRUE){
  p <- ncol(Sigma)
  a <- Inf
  n <- 0
  Ds <- list()
  D <- matrix(nrow = p, ncol = p, 0)
  ix <- (1: p^2 )[B != 0]
  E <- matrix(nrow = p, ncol = p, 0)
  while (a > eps && n < maxIter){
    if (r) ix <- (1: p^2 )[B != 0]
    n <- n + 1
    allres <- clyap(A = B, Q = C, all = TRUE)
    S <- matrix(nrow = p, data = allres[[6]])
    AA <- matrix(nrow = p, data = allres[[4]])
    EE <- matrix(nrow = p, data = allres[[5]])
    WKV <- allres[[7]]
    Delta <- Sigma - S

    # for (i in 1:length(ix)){
    #   E[ix[i]] <- 1
    #   Cp <- E %*% S + S %*% t(E)
    #   Ds[[i]] <- clyap2(A = AA, Q = Cp, E = EE, WKV = WKV)
    #   E[ix[i]] <- 0
    # }

    u <- vapply(1:length(ix), function(i){
        E[ix[i]] <- 1
        Cp <- E %*% S + S %*% t(E)
        D <- clyap2(A = AA, Q = Cp, E = EE, WKV = WKV)
        E[ix[i]] <- 0
        sum(D * Delta )
      }, FUN.VALUE = 1)

    #U <- (W - 0.5 * (B %*% Delta + Delta %*% t(B)) ) %*% solve(S)


    Bold <- B[ix]

    B[ix] <- B[ix] + alpha * u

    B[ix] <- sign(B[ix]) * (abs(B[ix]) - alpha * lambda)
    B[abs(B) < (alpha * lambda)] <- 0

    a <- sqrt(sum( (B[ix] - Bold)^2 ) / alpha)

    if (trace > 1){
      message("Iteration: ", n, " ||u||:", signif(a), " alpha:", alpha)
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||u||=", signif(a))
  }

  return(B)
}


#' Estimate B minimizing Frobenius distance from inverse covariance
#'
#' @param Sigma the empirical covariance matrix
#' @param B an initial guess for matrix B
#' @param C the value of the noise matrix C
#' @param eps stopping criteria
#' @param alpha step length
#' @param beta increasing factor for the step size
#' @param maxIter
#' @param trace
#' @param lambda
#'
#' @return the estimated B matrix
#' @importFrom lyapunov clyap clyap2
#' @export
estimateBP <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                       alpha = 1e-2, beta = 2, maxIter = Inf, trace = 0,
                       lambda = 0.1, r = TRUE){
  Pt <- solve(Sigma)
  p <- ncol(Sigma)
  a <- Inf
  n <- 0
  aold <- 0
  Ds <- list()
  D <- matrix(nrow = p, ncol = p, 0)
  E <- matrix(nrow = p, ncol = p, 0)
  ix <- (1: p^2 )[B != 0]
  while (a > eps && n < maxIter){
    if (r) ix <- (1: p^2 )[B != 0]

    n <- n + 1


    allres <- clyap(A = B, Q = C, all = TRUE)
    S <- matrix(nrow = p, data = allres[[6]])
    AA <- matrix(nrow = p, data = allres[[4]])
    EE <- matrix(nrow = p, data = allres[[5]])
    WKV <- allres[[7]]

    P <- solve(S)
    Delta <- Pt - P

    # for (i in 1:length(ix)){
    #   E <- matrix(nrow = p, ncol = p, 0)
    #   E[ix[i]] <- 1
    #   Cp <- E %*% S + S %*% t(E)
    #   Ds[[i]] <- - P %*% clyap2(A = AA, Q = Cp, E = EE,  WKV = WKV) %*% P
    # }
    #
    # u <- vapply(1:length(ix), function(i) sum(Ds[[i]] * Delta ),
    #             FUN.VALUE = 1)

    u <- vapply(1:length(ix), function(i){
      E[ix[i]] <- 1
      Cp <- E %*% S + S %*% t(E)
      D <- - P %*% clyap2(A = AA, Q = Cp, E = EE,  WKV = WKV) %*% P
      E[ix[i]] <- 0
      sum(D * Delta )
    }, FUN.VALUE = 1)

    Bold <- B[ix]

    B[ix] <- B[ix] +  alpha * u
    B[ix] <- sign(B[ix]) * (abs(B[ix]) - alpha * lambda)
    B[abs(B) < (alpha * lambda)] <- 0

    a <- sqrt(sum((B[ix] - Bold)^2) / alpha)

    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||u||=", signif(a))
  }

  return(B)
}


#' Estimate B minimizing penalized logLik
#'
#' @param Sigma the empirical covariance matrix
#' @param B an initial guess for matrix B
#' @param C the value of the noise matrix C
#' @param eps stopping criteria
#' @param alpha step size
#' @param maxIter maximum number of iterations
#' @param trace if >0 print info
#' @param lambda sparsity
#' @param r logical
#'
#' @return the estimated B matrix
#' @importFrom lyapunov clyap
#' @export
estimateBLL <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                      alpha = 0.1, maxIter = 1000, trace = 0,
                      lambda = 0, r = FALSE){
  p <- ncol(Sigma)

  a <- Inf
  n <- 0
  E <- matrix(nrow = p, ncol = p, 0)
  ix <- (1: p^2 )[B != 0]
  while (a > eps && n < maxIter  || n < 2){
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

    ### computing gradient
    Ds <- list()

    ### obtain DS
    # for (i in 1:length(ix)){
    #   E[ix[i]] <- 1
    #   Cp <- E %*% S + S %*% t(E)
    #   Ds[[i]] <- clyap2(A = AA, Q = Cp, E = EE, WKV = WKV)
    #   E[ix[i]] <- 0
    # }

    ###compute gradient
    # for (i in 1:length(ix)){
    #   u[i] <- sum(tmp * Ds[[i]] )
    # }

    u <- vapply(1:length(ix), function(i){
      E[ix[i]] <- 1
      Cp <- E %*% S + S %*% t(E)
      D <- clyap2(A = AA, Q = Cp, E = EE,  WKV = WKV)
      E[ix[i]] <- 0
      sum(tmp * D )
    }, FUN.VALUE = 1)

    Bold <- B

    ### gradient step
    B[ix] <- B[ix] + alpha * u


    ### soft thres
    B[ix] <- sign(B[ix]) * (abs(B[ix]) - alpha * lambda)
    B[abs(B) < (alpha * lambda)] <- 0

    a <- sqrt(sum( (B - Bold)^2 ) / alpha)
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
#' @importFrom lyapunov clyap clyap2
#' @export
estimateBLL2 <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                        alpha = 0.1, maxIter = 1000, trace = 0,
                        lambda = 0, r = FALSE){
  p <- ncol(Sigma)
  Pt <- solve(Sigma)
  a <- Inf
  n <- 0
  ix <- (1: p^2 )[B != 0]
  while (a > eps && n < maxIter  || n < 2){
    if (r) ix <- (1: p^2 )[B != 0]
    u <- rep(0, length(ix))
    n <- n + 1

    S <- clyap(A = B, Q = C)
    P <- solve(S)

    #delta <- S - Sigma
    delta <- P - Pt

    Cd <- - (B %*% delta + delta %*% t(B))


    D <- -  Cd %*% P
    Bold <- B

    ### gradient step
    B <- B + alpha * D


    ### soft thres
    B[ix] <- sign(B[ix]) * (abs(B[ix]) - alpha * lambda)
    B[abs(B) < (alpha * lambda)] <- 0

    a <- sqrt(sum( (B - Bold)^2 ) / alpha)
    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||diff||=", signif(a))
  }

  return(B)
}

