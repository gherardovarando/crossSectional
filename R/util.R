

#'
#' @param B1 matrix
#' @param B2 matrix
#'
#' @return a graph with colored edges where
#'         in red are plotted the edges in \code{B1} and not in
#'         \code{B2}
#' @import igraph
#' @export
deltaGraph <- function(B1, B2, legend = FALSE,  ...){
  gr1 <- graph_from_adjacency_matrix(abs(sign(t(B1))))
  gr2 <- graph_from_adjacency_matrix(abs(sign(t(B2))))
  d1 <-  set_edge_attr(gr1 - gr2, name = "color", index = , value = "red" )
  d1 <- add_edges(d1, edges = c(t(get.edgelist(gr2 - gr1))),
                  attr = list(color = "blue"))
  ed <- c(t(get.edgelist(intersection(gr1, gr2))))
  d1 <- add_edges(d1, edges = ed,
                  attr = list(color = "black"))
  plot(d1, ...)
  if (legend){
    legend("left", legend = c("correct", "missing", "wrong"),
           col = c("black", "red", "blue"), lty = 1, lwd = 2 ,cex = 0.5)
  }
  return(invisible(d1))
}


#' hamming distance
#'
#' @param A
#' @param B
#'
#' @export
hamming <- function(A, B){
  sum( (A != 0) != (B != 0) )
}


#' Minus Log-Likelihood for Gaussian model
#'
#' @param P the inverse of the covariance matrix
#' @param S the empirical covariance
#'
#' @importFrom stats var
#' @export
mll <- function(P, S){
 -determinant(P, logarithm = TRUE)$modulus + sum(S * P)
}

#' Minus log-likelihood for invariant distribution of OU
#'
#' @param B
#' @param S
#' @param C
#'
#' @export
mllB <- function(B, S, C = diag(nrow(B))){
  P <- solve(lyapunov::clyap(A = B, Q = C))
  mll(P, S)
}


#' Generate data from invariant distribution of O-U process
#'
#' Sample observations from the invariant distribution of an
#' Ornstein-Uhlenbeck process:  \deqn{dX_t = B(X_t - \mu)dt + DdW_t}
#'
#' @param n number of sample
#' @param B square matrix, coefficients of the O-U equation
#' @param D noise matrix
#' @param mean invariant mean
#'
#' @return  A list with components:
#'   -  \code{sample} the produced sample
#'   -  \code{Sigma} the covariance matrix of the invariant distribution
#'   -  \code{C} the  C matrix
#'   -  \code{B} the coefficient square matrix
#'   -  \code{mean} the mean vector of the invariant distribution
#'
#'
#' @importFrom lyapunov clyap
#' @importFrom MASS mvrnorm
#'
#' @export
rOU <- function(n = 1, B,
                D = diag(nrow = nrow(B), ncol = ncol(B)),
                mean = rep(0, nrow(B)) ){
  C <- D %*% t(D)
  Sigma <- lyapunov::clyap(A = B, Q = C)
  S <- MASS::mvrnorm(n = n, Sigma = Sigma, mu = mean)
  return(list(data = S, Sigma = Sigma, C = C, B = B, mean = mean))
}


#' Random Stable (Hurwitz) matrix
#'
#' Generate a random stable matrix, that is with all
#' eigenvalues with strict negative real part
#'
#' @param p integer, the dimension of the desired matrix
#' @param d fraction of non-zero entries
#'
#' @importFrom gmat chol_mh
#' @export
rHurwitz <- function(p = 1, d = 1){
  ## each stable matrix is the product of a SPD and a generalized negative definite
  S <- gmat::chol_mh(N = 1, p = p, d = d)[,,1]
  An <- rSkewSymm(p = p)
  Sy <- - gmat::chol_mh(N = 1, p = p, d = d)[,,1]
  return( (An + Sy) %*% S)
}


#' Generate a random skew-symmetric matrix
#'
#' @param p integer the dimension of the desired matrix
#' @param rFun a function with signature \code{function(n)} that
#' generates \code{n} numbers (e.g. \code{rnorm})
#'
#' @return a skew-symmetric matrix
#'
#' @importFrom stats rnorm
#' @export
rSkewSymm <- function(p, rFun = rnorm){
  M <- matrix(nrow = p, ncol = p, data = rFun(p * p))
  return(M - t(M) )
}

#' Generate stable Metzler matrix
#'
#' @param p integere the dimension of the matrix
#' @param d proportion of non zero entries
#' @param lower logical, should the matrix be lower triangular
#'
#' @return a square matrix with positive off-diagonal entries (Metzler)
#' and eigenvalues with strict negative real part (stable)
#'
#' @importFrom stats rnorm
#' @export
rStableMetzler <- function(p = 1, d = 1, lower = FALSE, rfun = rnorm, rdiag = rnorm){
  D <- abs(rfun(p * p)) * sample(c(0,1), replace = TRUE, size = p * p,
                             prob = c(1 - d, d))
  A <- matrix(nrow = p, ncol = p, data = D)
  diag(A) <- 0
  if (lower){
    A[upper.tri(A)] <- 0
  }
  diag(A) <- -rowSums(A) - abs(rdiag(p))
  return(A)
}


#' Recover lower triangular B
#'
#' Recover the only lower triangular stable matrix B such that
#'
#'
#'@param Sigma PSD matrix, the covariance matrix
#'@param P PSD matrix, the inverse of the covariance
#'@param C PSD matrix
#'
#'@return A stable lower triangular matrix
#'
#'@importFrom gmat anti_t
#'@export
lowertriangB <- function(Sigma,
                         P = solve(Sigma),
                         C = diag(nrow = nrow(Sigma))){

  p <- nrow(Sigma)
  l <- listInverseBlocks(Sigma)
  blong <- gmat::anti_t(0.5 * C %*% P)[ upper.tri(Sigma)]
  blong <- blong[length(blong):1]
  b <- blong[1:(p-1)]
  #b <- 0.5 * C[1,1] * P[1, 2:p]
  w <- l[[1]] %*% b
  t <- p
  for (i in 2:length(l)){
    #b <- 0.5 * C[i,i] * P[i, (i+1):p ]
    b <- blong[t:(t + p - i - 1)]
    t <- t + p - i
    AA <- matrix(nrow = length(b), ncol = length(w), 0)
    for (j in (i+1):p){
      for (k in 1:(i-1)){
        ix <- (k - 1) * p - (k - 1) * k / 2 + (i - k)
        AA[j - i, ix] <-  P[k, j]
      }
    }
    b <- b + tcrossprod(AA, t(w) )
    w <- c(w, l[[i]] %*% b)
  }
  W <- matrix(nrow = p, ncol = p, 0)
  W[upper.tri(W)] <- w[length(w) : 1]
  W <- gmat::anti_t(W)
  W <- W - t(W)
  B <- (W - 0.5 * C) %*% P
  B[upper.tri(B)] <- 0
  return(B)
}


#' List all the inverse matrices of the leading sub-matrix of P=Sigma^{-1}
#' (INTERNAL)
#'
#'@param Sigma an invertible (SYMMETRIC??) matrix
#'
#'@return a list with the inverse of the leading sub-matrices of Sigma^(-1)
listInverseBlocks <- function(Sigma){
   l <- list()
   l[[1]] <- Sigma
   for (i in 2:(nrow(Sigma) )){
     l[[i]] <-l[[i - 1]][ - 1, - 1] - (l[[i - 1]][ - 1, 1] %*%
                                              t(l[[i - 1]][1, - 1])) /
       l[[i - 1]][1, 1]
   }
   return(l[-1])
}



#' Maximum likelihood optimization
#'
#' @param B initial coefficient matrix
#' @param S empirical covariance matrix
#' @param C positive definite term of the Lyapunov equation
#' @param ... additional parameters passed to \code{optim}
#'
#' @detail MLE for coefficient matrix with the same 0s pattern of \code{B}
#'
#' @return The MLE estimation of the coefficient matrix
#'
#' @export
optimizeB <- function(B, S, C = diag(nrow(B)),
                      ...){
   ix <- B != 0
   p <- nrow(B)
   f_temp <- function(pars){
     BB <- matrix(nrow = p, ncol = p, 0)
     BB[ix] <- pars
     if (any(Re(eigen(BB)$values) > 0)){
       return(Inf)
     }
     P <- solve(lyapunov::clyap(A = BB, Q = C))
     return(mll(P = P, S = S))
   }
   res  <- optim(par = B[ix], fn = f_temp, ...)
   B[ix] <- res$par
   attr(B, "objective.value") <- res$value
   attr(B, "convergence") <- res$convergence
   return(B)
}


#' Maximum likelihood optimization for triangular B
#'
#' @param B initial coefficient matrix
#' @param S empirical covariance matrix
#' @param C positive definite term of the Lyapunov equation
#' @param ... additional parameters passed to \code{optim}
#'
#' @detail MLE for coefficient matrix with the same 0s pattern of \code{B}
#'
#' @return The MLE estimation of the coefficient matrix
#'
#' @export
optimizeBtriang <- function(B, S, C = diag(nrow(B)),
                      ...){
  ix <- lower.tri(B, diag = TRUE) & (B!=0)
  p <- nrow(B)
  f_temp <- function(pars){
    BB <- matrix(nrow = p, ncol = p, 0)
    BB[ix] <- pars
    if (any(diag(BB) > 0)){
      return(Inf)
    }
    P <- solve(lyapunov::clyap(A = BB, Q = C))
    return(mll(P = P, S = S))
  }
  res  <- optim(par = B[ix], fn = f_temp, ...)
  B[ix] <- res$par
  attr(B, "objective.value") <- res$value
  attr(B, "convergence") <- res$convergence
  return(B)
}


#'
#'@param B
#'@param Sigma
#'@export
#'@importFrom lyapunov clyap
maxllB <- function(B, Sigma, C = diag(ncol(B)) ){
  P <- solve(Sigma)
  ix <- (1: p^2 )[B != 0 ]
  p <- nrow(B)
  done <- FALSE
  delta <- 1
  while (delta > 0.0001){
    obj <- mllB(B, Sigma, C)
    for (i in ix){
      S <- clyap(B, C)
      U <- matrix(nrow = p, ncol = p , 0)
      U[i] <- 1
      C1 <- U %*% S + S %*% t(U)
      D <- clyap(B, C1)
      f <- function(h){
       mll(solve(S + h * D), Sigma)
      }

      h <- optimise(f, c(-1,1))$minimum

      Btry <- B + h * U
      objn <- mllB(Btry, Sigma, C)
      if (objn < obj){
         delta <- objn - obj
         B <- Btry
         obj <- objn
      }
    }
  }

  return(B)
}




