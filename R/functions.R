score_psi <- function(ZtZ, ZtX, Zte, H, Psi0, psi0, finf = TRUE)
{
  # Define dimensions
  p <- ncol(XtX)
  q <- ncol(ZtZ)

  # Pre-compute
  Psi0ZtZ <- tcrossprod(Psi0, ZtZ)

  # Construct M
  M <- Psi0ZtZ
  diag(M) <- diag(M) + 1
  M <- solve(M, Psi0ZtZ) # Fast because sparse?



}
