#' @export
score_psi <- function(Z, ZtZXe, e, H, Psi0, psi0, finf = TRUE)
{
  # Define dimensions
  n <- length(e)
  q <- ncol(Psi0)
  r <- ncol(H) / q # Assumes H = [H_1, ... , H_r], where H_j is q by q
  p <- ncol(ZtZXe) - q - 1
  # Score vector to return
  s_psi <- rep(0, r + 1)

  # Fisher information to return
  I_psi <- matrix(NA, r + 1, r + 1)

  # Pre-compute Psi0ZtZ (columns 1:q), Pzi0ZtX (columns (q + 1):(q + p)),
  # and Psi0Zte (column q + p + 1)
  A <- Matrix::crossprod(Psi0, ZtZXe)

  # Matrix denoted M in manuscript is A[, 1:q]
  A <- Matrix::solve(A[, 1:q] + Matrix::Diagonal(q), A, sparse = TRUE)

  # Score for error variance psi_0
  # NB: REPLACE e by Sigma^{-1}e
  e <- (1 / psi0) * (e - Z %*% A[, q + p + 1]) # = Sigma^{-1}e
  trace_M <- sum(Matrix::diag(A[, 1:q]))
  s_psi[1] <- 0.5 * sum(e^2) - (0.5 / psi0) * (n - trace_M)


  # Use recycling to compute v'H_i v for all Hi
  v <- as(Matrix::crossprod(Z, e), "sparseVector") # sparse matrix does not recycle
  s_psi[-1] <- 0.5 * colSums(matrix(as.vector(Matrix::crossprod(v, H) * v), nrow = q))

  # B = Z'Z (M - I_q) in paper notation
  B <- A[, 1:q] - Matrix::Diagonal(q)
  B <- Matrix::crossprod(ZtZXe[, 1:q], B)

  if(!finf){
    # Compute -[ZtZ (I_q - M) * H_1, ..., ZtZ (I_q - M) * H_r] using
    # recycling, where * denotes elementwise multiplication
    # The "if" is because the calculation is a byproduct of a more expensive one
    # (B %*% H) done to get Fisher information
    H <- H * as(B, "sparseVector")

    s_psi[-1] <- s_psi[-1] + (0.5 / psi0) * colSums(matrix(Matrix::colSums(H), nrow = q))
  } else{
    H <- B %*% H
    I_psi[1, 1] <- (0.5 / psi0^2) * (n - 2 * trace_M +
    sum(Matrix::t(A[, 1:q]) * A[, 1:q]))

    D <- as(as(A[, 1:q] - Matrix::Diagonal(q), "sparseVector") * as(H, "sparseVector"), "sparseMatrix")
    dim(D) <- c(nrow(D) / r, r)
    I_psi[1, -1] <- (0.5 / psi0^2) *Matrix::colSums(D)

    for(ii in 1:r){
      first_idx <- ((ii - 1) * q + 1):(ii * q)
      s_psi[1 + ii] <- s_psi[1 + ii] + (0.5 / psi0) * sum(Matrix::diag(H[, first_idx]))
      for(jj in ii:r){
        second_idx <- ((jj - 1) * q + 1):(jj * q)
        I_psi[ii + 1, jj + 1] <- (0.5 / psi0^2) * sum(Matrix::t(H[, second_idx]) * H[, first_idx])
      }
    }
  }

  I_psi <- Matrix::forceSymmetric(I_psi, uplo = "U")
  return(list("score" = s_psi, "finf" = I_psi))

}
#' @export
loglik <- function(ZtZ, Zte, e, Psi0, psi0){
  n <- length(e)
  B <- Matrix::crossprod(Psi0, ZtZ) + Matrix::Diagonal(ncol(ZtZ))
  two_neg_ll <- n * log(psi0) + Matrix::determinant(B, logarithm = TRUE)$modulus
  two_neg_ll <- two_neg_ll + sum(e^2)/psi0

  # Profiling suggests this is faster than using
  # R <- Matrix::Cholesky(B) and then solve(R, ...) here
  two_neg_ll <- two_neg_ll - Matrix::crossprod(Zte, Matrix::solve(B,
                                            Matrix::crossprod(Psi0, Zte))) / psi0

  return(-0.5 * as.vector(two_neg_ll))

}

#' @export
get_Psi <- function(psi, H){
  q <- nrow(H)
  r <- ncol(H) / q
  Psi <- matrix(0, q, q)
  for(ii in 1:r){
    Psi <- Psi + psi[ii] * H[, ((ii - 1) * q + 1):(ii * q)]
  }
  Psi
}

#' @export
uni_test_stat <- function(test_seq, test_idx, psi, psi0, Z, ZtZXe, e, H)
{
  # first element of test_seq has to agree with psi[test_idx]
  m <- length(test_seq)
  test_stat <- rep(0, m)
  for(ii in 1:m){
    if(ii > 1){
      # Update tested parameter
      psi[test_idx] <- test_seq[ii]
      # Do one-step Fisher scoring update (with previous Information)
      # for non-tested parameters
      Psi <- getPsi(psi, H)
      s <- score_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                     Psi0 = Psi / psi0,
                     psi0 = psi0, finf = FALSE)
      psi[-test_idx] <- psi[-test_idx] +
        solve(score_inf$finf[-test_idx. -test_idx], s)
    }

    Psi <- get_Psi(psi, H)
    score_inf <- score_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                           Psi0 = Psi / psi0,
                           psi0 = psi0, finf = TRUE)

    eff_inf <- score_inf$finf[test_idx, test_idx] -
      sum(solve(score_inf$finf[-test_idx, -test_idx],
                score_inf$finf[-test_idx, test_idx]) *
            score_inf$finf[-test_idx, test_idx])

    test_stat[ii] <- score_inf$score[test_idx]^2 / eff_inf
  }
  test_stat
}
