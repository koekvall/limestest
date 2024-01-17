# change of original loglik: ll -> ll[1], observed I "-"->"+"
# function with psi1 as interest, (psi0, psi2) optimize by trust
est_psi <- function(Z, ZtZXe, e, H, psi1) #loglik = TRUE,score = TRUE, finf = TRUE, expected = FALSE
{
  # Define dimensions
  n <- length(e)
  q <- ncol(Z)
  r <- ncol(H) / q # Assumes H = [H_1, ... , H_r], where H_j is q by q
  p <- ncol(ZtZXe) - q - 1
  m <- q / 2  # check

  objfun <- function(pars) { # pars = c(psi0, psi2)

    Psi1 <- matrix(c(psi1,0,0,pars[2]), nrow = 2)
    psi0 <- pars[1]
    Psi0 <- Matrix::kronecker(Matrix::Diagonal(m), Psi1 / psi0)

    # loglikelihood to return
    ll <- NA
    # Score vector for trust, assuming for psi0, par of interest, others..
    s_psi <- rep(NA, r + 1)
    # Hessian to return
    I_psi <- matrix(NA, r + 1, r + 1)

    # Pre-compute Psi0ZtZ (columns 1:q), Pzi0ZtX (columns (q + 1):(q + p)),
    # and Psi0Zte (column q + p + 1)
    A <- Matrix::crossprod(Psi0, ZtZXe)

    eig <- eigen(A[, 1:q] + Matrix::Diagonal(q))
    # check if A + I is p.d.
    if (eig$values[q] <= 0) {
      ll <- -Inf
    } else {
      # Add loglik term before overwriting
      ll <- -0.5 * Matrix::determinant(A[, 1:q] + Matrix::Diagonal(q))$modulus[1] -
        0.5 * n * log(psi0)

      # Matrix denoted M in manuscript is A[, 1:q]
      A <- Matrix::solve(A[, 1:q] + Matrix::Diagonal(q), A, sparse = TRUE)
    }





    # Score for error variance psi_0
    # NB: REPLACE e by Sigma^{-1}e
    e_save <- e
    e <- (1 / psi0) * (e - Z %*% A[, q + p + 1]) # = Sigma^{-1}e
    ll <- ll  - 0.5 * sum(e * e_save)

    trace_M <- sum(Matrix::diag(A[, 1:q]))
    s_psi[1] <- 0.5 * sum(e^2) - (0.5 / psi0) * (n - trace_M)

    # Use recycling to compute v'H_i v for all Hi
    v <- as(Matrix::crossprod(Z, e), "sparseVector") # sparse matrix does not recycle
    w <- Matrix::crossprod(v, H) # Used later if !expected
    s_psi[-1] <- 0.5 * colSums(matrix(as.vector(w * v), nrow = q))

    # B = Z'Z (M - I_q) in paper notation
    B <- A[, 1:q]
    Matrix::diag(B) <- Matrix::diag(B) - 1
    B <- Matrix::crossprod(ZtZXe[, 1:q], B)

    # if finf
    H <- B %*% H
    I_psi[1, 1] <- (0.5 / psi0^2) * (n - 2 * trace_M +
                                       sum(Matrix::t(A[, 1:q]) * A[, 1:q]))
    # Subtract identity matrix from M
    Matrix::diag(A[, 1:q]) <- Matrix::diag(A[, 1:q]) - 1
    D <- matrix(Matrix::colSums(as.vector(A[, 1:q])  * as.matrix(H)), nrow = q)
    I_psi[1, -1] <- (0.5 / psi0^2) * Matrix::colSums(D)

    for(ii in 1:r){
      first_idx <- ((ii - 1) * q + 1):(ii * q)
      s_psi[1 + ii] <- s_psi[1 + ii] + (0.5 / psi0) * sum(Matrix::diag(H[, first_idx]))
      for(jj in ii:r){
        second_idx <- ((jj - 1) * q + 1):(jj * q)
        I_psi[ii + 1, jj + 1] <- (0.5 / psi0^2) * sum(Matrix::t(H[, second_idx]) * H[, first_idx])
      }
    }
    #  !expect
    I_psi <- -I_psi
    # u = Sigma^{-2}e. Some calculations could be saved from before
    u <- (1 / psi0^2) * (e_save + Z %*% (-2 * A[, q + p + 1] +
                                           (A[, 1:q] + Matrix::Diagonal(q, 1)) %*% A[, q + p + 1]))
    I_psi[1, 1] <- I_psi[1, 1] + sum(e * u)
    v <- as.vector(Matrix::crossprod(Z, u)) # = Z' Sigma^{-2}e
    I_psi[1, -1] <- I_psi[1, -1] + colSums(matrix(v * w, ncol = r))

    for(ii in 1:r){
      first_idx <- ((ii - 1) * q + 1):(ii * q)
      for(jj in ii:r){
        second_idx <- ((jj - 1) * q + 1):(jj * q)
        I_psi[ii + 1, jj + 1] <- I_psi[ii + 1, jj + 1] + (1 / psi0) * sum(crossprod(w[first_idx], B) * w[second_idx])
      }
    }
    I_psi <- as.matrix(-Matrix::forceSymmetric(I_psi, uplo = "U"))
    # psi1 (par of interest) not need to calculate derivatives
    s_psi <- s_psi[-2]
    I_psi <- I_psi[-2,-2]
    result <- list(value = ll, gradient = s_psi, hessian = I_psi)
    result
  }
  obj <- trust(objfun, parinit = c(1,0.01), rinit = 1, rmax = 10, minimize = F)
  #print(obj)
  Psi1 <- matrix(c(psi1,0,0,obj$argument[2]), nrow = 2)
  psi0 <- obj$argument[1]
  return(list(psi0 = psi0, Psi1 = Psi1))
}




