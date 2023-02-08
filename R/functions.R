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
  v <- as.vector(Matrix::crossprod(Z, e)) # sparse matrix does not recycle
  s_psi[-1] <- 0.5 * colSums(matrix(as.vector(Matrix::crossprod(v, H)), nrow = q) * v)

  B <- A[, 1:q] - Matrix::Diagonal(q)
  B <- Matrix::crossprod(ZtZXe[, 1:q], B)

  # This wastes memory; done to use recycling of b in next step, which does
  # not currently work with B
  dim(B) <- c(prod(dim(B)), 1)
  b <- Matrix::sparseVector(x = B@x, i = B@i + 1, length = dim(B)[1])

  ## Below should equal -[ZtZ (I_q - M) * H_1, ..., ZtZ (I_q - M) * H_r] by
  ## recycling, where * denotes elementwise multiplication
  H <- H * b
  s_psi[-1] <- s_psi[-1] + (0.5 / psi0) * colSums(matrix(Matrix::colSums(H), nrow = q))
  return(s_psi)

  # May be avoided by never changing dim(B) above
  dim(B) <- c(q, q)
  I_psi[1, 1] <- (1 / psi0^2) / (n - 2 * trace_M + sum(t(A[, 1:q]) *
                                Matrix::crossprod(ZtZXe[, 1:q], A[, 1:q])))


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


