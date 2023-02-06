score_psi <- function(ZtZ, ZtX, Zte, e, H, Psi0, psi0, finf = TRUE)
{
  # Define dimensions
  p <- ncol(XtX)
  q <- ncol(ZtZ)
  r <- ncol(H) / q # Assumes H = [H_1, ... , H_r], where H_j is q by q

  # Score vector to return
  s_psi <- rep(0, r + 1)

  # Fisher information to return
  I_psi <- matrix(NA, r + 1, r + 1)

  # Pre-compute Psi0ZtZ (columns 1:r), Pzi0ZtX (columns (r + 1):(r + p)),
  # and Psi0Zte (column r + p + 1)
  A <- crossprod(Psi0, cbind(ZtZ, ZtX, Zte))

  # MAIN COMPUTATION? PROFILE
  # Sparse solve needed for all quantities
  B <- Psi0ZtZXe[, 1:r] # = Psi0ZtZ
  diag(B) <- diag(B) + 1
  A <- solve(B, A)   # Matrix denoted M in manuscript is A[, 1:r]

  # Score for error variance psi_0
  e <- (e - Z %*% A[, p + r + 1]) / psi0 # Replace residual vector
  trace_M <- sum(diag(A[, 1:r]))
  s_psi[1] <- (sum(e^2) - (n - trace_M) / psi0) / 2

  # Score for parameter vector psi
  v = crossprod(Z, e)

  ## the right-most elementwise multiplication uses recycling; make sure correct
  ## assumes H = [H_1, ..., H_r]
  s_psi[-1] <- 0.5 * colSums(matrix(as.vector(crossprod(v, H)), nrow = q) * v)

  B <- A[, 1:r] # This copying can probably be avoided but may not affect speed
  diag(B) <- diag(B) - 1
  B <- crossprod(ZtZ, B)
  ## H * B should now equal -[ZtZ (I_q - M) * H_1, ..., ZtZ (I_q - M) * H_r] by
  ## recycling, where * denotes elementwise multiplication
  s_psi[-1] <- s_psi[-1] + (0.5 / psi0) * colSums(matrix(colSums(H * B),
                                                         nrow = q))
  return(s_psi)
}

loglik <- function(ZtZ, Zte, e, Psi0, psi0){
  n <- length(e)
  R <- chol(crossprod(Psi0, ZtZ), Imult = 1)
  two_neg_ll <- n * log(psi0) + sum(log(diag(R)))
  two_neg_ll <- two_neg_ll + sum(e^2)/psi0
  # Profiling suggests this is faster than using the cholesky factor
  two_neg_ll <- two_neg_ll - crossprod(Zte, solve(R, solve(t(R),
                                            crossprod(Psi0, Zte))))/psi0

  return(-0.5 * two_neg_ll)

}
