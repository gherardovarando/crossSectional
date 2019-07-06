#' Estimate Ornstein-Uhlenbeck process
#'
#' Estimate the coefficients and the noise term of an Ornstein-Uhlenbeck
#' process (\eqn{dX(t) = BX(t) + \sqrt{C} dW(t)}) using penalized maximum-likelihood.
#'
#' @param Sigma the empirical covariance matrix
#' @param B an initial guess for the coefficient matrix B
#' @param C the initial guess for the noise matrix C
#' @param eps stopping criteria
#' @param alpha parameter backtracking
#' @param beta parameter backtracking
#' @param delt stepsize for internal stable approx
#' @param maxIter maximum number of iterations
#' @param trace if >0 print info
#' @param lambda penalization coefficient
#' @param r logical the set of parameter will be updated at every iteration
#' @param h logical if TRUE only the non-zero entries of B will be updated
#'
#' @return the estimated B matrix (\code{estimateBLL}) or
#' the estiamted C matrix (\code{estiamteCLL}).
#' @importFrom lyapunov clyap clyap2
#' @export
estimateBLL <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                        alpha = 0.2, beta = 0.5,
                        maxIter = 1000, trace = 0,
                        lambda = 0, r = FALSE, h = FALSE, pert = FALSE){
  p <- ncol(Sigma)
  ixd <- 0:(p - 1) * (p) + 1:p ##index of diagonal elements
  a <- Inf
  n <- 0
  E <- matrix(nrow = p, ncol = p, 0)
  if (h) ix <- (1: p^2 )[B != 0]
  else ix <- 1:(p^2)
  allres <- clyap(A = B, Q = C, all = TRUE)
  S <- matrix(nrow = p, data = allres[[6]])
  AA <- matrix(nrow = p, data = allres[[4]])
  EE <- matrix(nrow = p, data = allres[[5]])
  WKV <- allres[[7]]
  P <- solve(S)
  tk <- 1
  u <- rep(0, length(ix))
  f <- mll(P = P, S = Sigma)
  while (a > eps && n < maxIter) {
    if (r) ix <- (1: p^2 )[B != 0]
    ixnd <- ix[! ix %in% ixd ] ## not diagonal elements

    n <- n + 1

    tmp <- P %*% Sigma %*% P - P

    delta <- Sigma - S
    Cp <- delta %*% (B + C %*% P)
    Cp <- 0.5 * (Cp + t(Cp))
    AA <- S %x% diag(p)
    ixl <- matrix(nrow = p, ncol = p, 1:p^2)[lower.tri(Sigma)]
    ixu <- matrix(nrow = p, ncol = p, 1:p^2)[upper.tri(Sigma)]
    AA[, ixl] <-   AA[ , ixl] - AA[, ixu]
    AA <- AA[, -c(ixu, ixd)]
    u <- lm.fit(AA, c(Cp %*% S), offset = FALSE)$residual

    u <- u[ix]
    u <- u / sqrt(sum(u ^ 2))
    #u <- vapply(1:length(ix), function(i){
    #  E[ix[i]] <- 1
    #  Cp <- E %*% S + S %*% t(E)
    #  D <- clyap2(A = AA, Q = Cp, E = EE,  WKV = WKV)
    #  E[ix[i]] <- 0
    #  sum(tmp * D )
    #}, FUN.VALUE = 1)

    if (pert){
      ru <- rnorm(length(u))
      u <- u + ru / sqrt(sum(ru ^ 2))
    }

    Bold <- B

    #### Beck and Teboulle line search
    f <- mll(P, Sigma) + lambda * sum(abs(Bold[ixnd]))
    fnew <- Inf

    alph <- alpha
    while ( fnew   > f - sum(u * (B[ix] - Bold[ix]) ) +
            sum((B[ix] - Bold[ix]) ^ 2) / (2* alph)  || fnew > f) {

      B[ix] <- Bold[ix] + alph * u

      ### soft thres
      B[ixnd] <- sign(B[ixnd]) * (abs(B[ixnd]) - alph * lambda)
      B[ixnd][abs(B[ixnd]) < (alph * lambda)] <- 0

      ### Lyapunv solution
      allres <- clyap(A = B, Q = C, all = TRUE)
      S <- matrix(nrow = p, data = allres[[6]])
      AA <- matrix(nrow = p, data = allres[[4]])
      EE <- matrix(nrow = p, data = allres[[5]])
      WKV <- allres[[7]]
      if (all( (diag(AA) * diag(EE))  < 0 )){
        P <- solve(S)
        fnew <- mll(P, Sigma) + lambda * sum(abs(B[ixnd]))
      }else{
        fnew <- Inf
      }
      alph <- alph * beta
    }

    a <- (f - fnew ) / (abs(f))
    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a),
              " alpha:", alph / beta)
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
                        alpha = 0.2,
                        beta = 0.5,
                        maxIter = 1000, trace = 0,
                        lambda = 0, r = FALSE,
                        t0 = 1){
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
  S <- clyap2(A = AA, Q = C, E = EE,  WKV = WKV)
  P <- solve(S)
  while (a > eps && n < maxIter){
    if (r) ix <- (1: p^2 )[B != 0]
    u <- rep(0, length(ix))
    n <- n + 1

    tmp <- P %*% Sigma %*% P - P

    v <- vapply(1:length(ixc), function(i){
      E[ixc[i]] <- 1
      Cp <- E + t(E)
      D <- clyap2(A = AA, Q = Cp, E = EE,  WKV = WKV)
      E[ixc[i]] <- 0
      sum(tmp * D )
    }, FUN.VALUE = 1) - 2 * lambda * ( (C - C0)[ixc])

    Cold <- C

    ### backtracking
    f <- mll(P, Sigma) + lambda * sum((C -C0)^2)
    fnew <- Inf
    t <- t0
    while (fnew > f - t * alpha * sum(v^2)){
      C[ixc] <- Cold[ixc] + t * v
      if (all(diag(C) > 0)){
        S <- clyap2(A = AA, Q = C, E = EE,  WKV = WKV)
        P <- solve(S)
        fnew <- mll(P, Sigma) + lambda * sum((C -C0)^2)
      }else{
        fnew <- Inf
      }
      t <- t * beta
    }

    #a <- beta * sqrt(sum( (C - Cold) ^ 2)) / t
    a <- (f - fnew) / abs(f)
    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||diff||=", signif(a))
  }

  return(C)
}


#' @rdname estimateBLL
#' @importFrom lyapunov clyap clyap2
#' @export
#' @importFrom glmnet glmnet
estimateBF <- function(Sigma,
                       B,
                       C = diag(ncol(Sigma)),
                       C0 = diag(ncol(Sigma)),
                       eps =  1e-2,
                       alpha = 1,
                       maxIter = 1000,
                       trace = 0,
                       lambda = 0.1,
                       beta = 0.5,
                       r = FALSE,
                       h = FALSE) {
  p <- ncol(Sigma)
  a <- Inf
  n <- 0
  E <- matrix(nrow = p, ncol = p, 0)

  if (h) ix <- (1: p^2 )[B != 0]
  else ix <- 1:(p^2)

  allres <- clyap(A = B, Q = C, all = TRUE)
  S <- matrix(nrow = p, data = allres[[6]])
  AA <- matrix(nrow = p, data = allres[[4]])
  EE <- matrix(nrow = p, data = allres[[5]])
  WKV <- allres[[7]]
  delta <- Sigma - S
  ixd <- 0:(p - 1) * (p) + 1:p ##index of diagonal elements
  while (abs(a) > eps && n < maxIter  || n < 2){
    if (r) ix <- (1: p^2)[B != 0]
    ixnd <- ix[! ix %in% ixd ] ## not diagonal elements

    u <- rep(0, length(ix))
    n <- n + 1
    Ds <- list()


    b <- rep(0, length(ix))
    M <- matrix(nrow = length(ix), ncol = length(ix), 0)
    for (i in 1:length(ix)){
      E[ix[i]] <- 1
      Cp <- E %*% S + S %*% t(E)
      Ds[[i]] <- clyap2(A = AA, Q = Cp, E = EE,  WKV = WKV)
      E[ix[i]] <- 0
      b[i] <- sum(Ds[[i]] * delta)
      for (j in 1:i){
        M[i, j] <- M[j, i] <- sum(Ds[[i]] * Ds[[j]])
      }
    }

    #u <- b
    u <- solve(M, b)
    #u <- lm.fit(M, b)$coefficients
    #u[is.na(u)] <- 0
    #u <- glmnet(M, y = b, alpha = 1, lambda = lambda)$beta

    Bold <- B

    #### Beck and Teboulle line search
    f <- sum(delta^2)
    alph <- alpha
    fnew <- Inf
    while ( (fnew > f - sum(u * (B[ix] - Bold[ix]) ) +
             sum((B[ix] - Bold[ix]) ^ 2) / (2* alph) )){

      B[ix] <- Bold[ix] + alph * u

      ### soft thres (prox)
      B[ixnd] <- sign(B[ixnd]) * (abs(B[ixnd]) - alph * lambda)
      B[ixnd][abs(B[ixnd]) < (alph * lambda)] <- 0
      ###

      allres <- clyap(A = B, Q = C, all = TRUE)
      S <- matrix(nrow = p, data = allres[[6]])
      AA <- matrix(nrow = p, data = allres[[4]])
      EE <- matrix(nrow = p, data = allres[[5]])
      WKV <- allres[[7]]

      delta <- Sigma - S
      fnew <- sum(delta ^ 2)
      alph <- beta * alph
    }
    a <- beta * sum( (Bold - B) ^ 2) / (alph)

    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||diff||=", signif(a))
  }

  return(B)
}


#' generalized graphical lasso
#'
#' Solve the generalized graphical lasso problem with proximal gradient
#' @param Sigma the covariance matriz
#' @param P initial precision matrix
#' @param G Matrix in the generlaized penalization
#' @param lambda penalization coefficient
#' @param maxIter maximum number of iterations
#' @param eps threshold for stopping criteria
#' @param alpha param for line search
#' @param beta param for line search
#' @return The estimated precision matrix
#' @export
#' @importFrom genlasso genlasso
genGlasso <- function(Sigma, P = diag(nrow(Sigma)),
                      G = -diag(nrow(Sigma)),
                      lambda = 0.1, maxIter = 1000,
                      eps = 1e-5, alpha = 0.5, beta = 0.5, trace = 0){

  p <- nrow(Sigma)
  a <- Inf
  ixd <- 0:(p - 1) * (p) + 1:p ##index of diagonal elements
  ixl <- matrix(nrow = p, ncol = p, 1:p^2)[lower.tri(Sigma)]
  ixu <- matrix(nrow = p, ncol = p, 1:p^2)[upper.tri(Sigma)]
  GG <- diag(p) %x% G ## kronecker product
  GG <- GG[-ixd, ] ## not penalty on diagonal of B
  GG[, ixl] <- GG[, ixu] + GG[ , ixl]
  GG <- GG[, -ixu] ## eliminate strict upper part
  n <- 0
  while (abs(a) > eps && n < maxIter){
    n <- n + 1
    Pold <- P
    S <- solve(P)
    u <- S - Sigma
    ### line search
    f <- mll(P, Sigma) + sum(abs(G %*% P))
    fnew <- Inf
    alph <- alpha
    while ( (fnew   > f - sum(u * (P - Pold) ) +
             sum((P - Pold) ^ 2) / (2* alph) ) ){

      P <- Pold + alph * u


      ### proximal step (genLasso)
      b <- genlasso(P[-ixu], D = GG, minlam = alph * lambda)$beta
      P[-ixu] <- b[, ncol(b)]
      P[ixu] <- P[ixl]
      if (all(eigen(P, symmetric = TRUE, only.values = TRUE)$values > 0)){
        fnew <- mll(P, Sigma) + sum(abs(G %*% P))
      }else{
        fnew <- Inf
      }

      alph <- alph * beta
      #break ### no line search
    }



    #a <- f + sum(abs(Pold)) - fnew - sum(abs(P))
    #a <- sqrt(sum(u ^ 2))
    a <- sqrt(sum((Pold - P)^2))
    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||diff||=", signif(a))
  }

  return(P)
}


