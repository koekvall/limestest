score_psi <- function(ZtZ, ZtX, Zte, e, H, Psi0, psi0, finf = TRUE)
{
  # Define dimensions
  p <- ncol(XtX)
  q <- ncol(ZtZ)

  # Score vector to return
  s_psi <- rep(0, q + 1)

  # Fisher information to return
  I_psi <- matrix(NA, q + 1, q + 1)

  # Pre-compute Psi0ZtZ (columns 1:q), Pzi0ZtX (columns (q + 1):(q + p)),
  # and Psi0Zte (column q + p + 1)
  A <- crossprod(Psi0, cbind(ZtZ, ZtX, Zte))

  # MAIN COMPUTATION? PROFILE
  # Sparse solve needed for all quantities
  B <- Psi0ZtZXe[, 1:q] # = Psi0ZtZ
  diag(B) <- diag(B) + 1
  A <- solve(B, A)   # Matrix denoted M in manuscript is A[, 1:q]

  # Score for error variance psi_0
  e <- e - Z %*% A[, p + q + 1] # Replace residual vector
  trace_M <- sum(diag(A[, 1:q]))
  s_psi[1] <- 2 * (sum(e^2)/psi0^2 - (n - trace_M) / psi0)

  if(finf){
    B <- crossprod(ZtZ, A[, 1:q]) # Re-use storage B
    I_psi[1, 1] <- (n - 2 * trace_M + sum(diag(crossprod(A[, 1:q], B)))) /
      psi0^2
  }

}
