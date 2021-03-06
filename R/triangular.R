#' Estimate triangular OU from cross sectional data
#'
#' Use l1
#'
#' @param S sample covariance
#' @param D noise matrix
#' @param alpha penalization coefficient
#' @param B initial estimate for the coefficient matrix
#' @param ... additional parameters passed to \code{optim}
#'
#' @export
l1B <-
  function(S,
           D = diag(ncol(S)),
           alpha = 0.1,
           B = NULL,
           ...) {
    p <- ncol(S)
    P <- solve(S)
    C <- D %*% t(D)
    ix <- lower.tri(C)
    if (is.null(B)) {
      B <- lowertriangB(S, P, C)
      m <- mll(P, S)
      B[upper.tri(B)] <- 0
    } else{
      B[upper.tri(B)] <- 0
      B <- optimizeBtriang(B, S, C, ...)
      m <- attr(B, "objective.value")
    }
    sorted <- sort(abs(B), decreasing = F, index.return = TRUE)$ix
    for (s in sorted) {
      j <- ceiling(s / p)
      i <- s - (j - 1) * p
      if (i > j) {
        Bnew <- B
        Bnew[i, j] <- 0
        Bnew <- optimizeBtriang(Bnew, S, C, ...)
        mnew <- attr(Bnew, "objective.value")
        if (mnew + alpha * sum(abs(Bnew[ix])) < m + (alpha) * sum(abs(B[ix]))) {
          m <- mnew
          B <- Bnew
        } else{
          #break
        }
      }
    }
    return(B)
  }


#' Estimate triangular OU from cross sectional data
#'
#' Use l1
#'
#' @param S the sample covariance matrix
#' @param D noise matrix
#' @param alpha penalization coefficient
#' @param B initial estimate for the coefficient matrix
#' @param ... additional parameters passed to \code{optim}
#'
#' @export
l1Bfast <-
  function(S,
           D = diag(ncol(S)),
           alpha = 0.1,
           B = NULL,
           ...) {
    p <- ncol(S)
    P <- solve(S)
    C <- D %*% t(D)
    ix <- lower.tri(C)
    if (is.null(B)) {
      B <- lowertriangB(S, P, C)
      m <- mll(P, S)
      B[upper.tri(B)] <- 0
    } else{
      B[upper.tri(B)] <- 0
      B <- optimizeBtriang(B, S, C, ...)
      m <- attr(B, "objective.value")
    }
    sorted <- sort(abs(B), decreasing = F, index.return = TRUE)$ix
    for (s in sorted) {
      j <- ceiling(s / p)
      i <- s - (j - 1) * p
      if (i > j) {
        Bnew <- B
        Bnew[i, j] <- 0
        #Bnew <- optimizeBtriang(Bnew, S, C, ...)
        mnew <- mllB(B = Bnew, S = S, C = C)
        if (mnew + (alpha ) * sum(abs(Bnew[ix])) < m + (alpha ) * sum(abs(B[ix]))) {
          m <- mnew
          B <- Bnew
        } else{
          #break
        }
      }
    }
    return(B)
  }

