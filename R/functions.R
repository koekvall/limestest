#' @export
loglik_psi <- function(Z, ZtZXe, e, H, Psi0, psi0, loglik = TRUE,
                       score = TRUE, finf = TRUE)
{
  # Define dimensions
  n <- length(e)
  q <- ncol(Psi0)
  r <- ncol(H) / q # Assumes H = [H_1, ... , H_r], where H_j is q by q
  p <- ncol(ZtZXe) - q - 1

  # loglikelihood to return
  ll <- NA

  # Score vector to return
  s_psi <- rep(NA, r + 1)

  # Fisher information to return
  I_psi <- matrix(NA, r + 1, r + 1)

  # Pre-compute Psi0ZtZ (columns 1:q), Pzi0ZtX (columns (q + 1):(q + p)),
  # and Psi0Zte (column q + p + 1)
  A <- Matrix::crossprod(Psi0, ZtZXe)

  # Add loglik term before overwriting
  if(loglik){
    ll <- Matrix::determinant(A[, 1:q] + Matrix::diag(q))$modulus +
      n * log(psi0)
  }

  # Matrix denoted M in manuscript is A[, 1:q]
  A <- Matrix::solve(A[, 1:q] + Matrix::Diagonal(q), A, sparse = TRUE)

  # Score for error variance psi_0
  # NB: REPLACE e by Sigma^{-1}e
  if(loglik){
   e_save <- e
  }
  e <- (1 / psi0) * (e - Z %*% A[, q + p + 1]) # = Sigma^{-1}e
  if(loglik){
    ll <- -0.5 * (ll + sum(e * e_save))
  }

  trace_M <- sum(Matrix::diag(A[, 1:q]))
  s_psi[1] <- 0.5 * sum(e^2) - (0.5 / psi0) * (n - trace_M)


  # Use recycling to compute v'H_i v for all Hi
  v <- as(Matrix::crossprod(Z, e), "sparseVector") # sparse matrix does not recycle
  s_psi[-1] <- 0.5 * colSums(matrix(as.vector(Matrix::crossprod(v, H) * v), nrow = q))

  # B = Z'Z (M - I_q) in paper notation
  B <- A[, 1:q]
  Matrix::diag(B) <- Matrix::diag(B) - 1
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

    # Subtract identity matrix from M
    Matrix::diag(A[, 1:q]) <- Matrix::diag(A[, 1:q]) - 1

    D <- matrix(Matrix::colSums(as(A[, 1:q], "sparseVector")  * H), nrow = q)

    I_psi[1, -1] <- (0.5 / psi0^2) * Matrix::colSums(D)

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
  return(list("ll" = ll,  "score" = s_psi, "finf" = I_psi))
}

# loglik <- function(ZtZ, Zte, e, Psi0, psi0){
#   n <- length(e)
#   B <- Matrix::crossprod(Psi0, ZtZ) + Matrix::Diagonal(ncol(ZtZ))
#   two_neg_ll <- n * log(psi0) + Matrix::determinant(B, logarithm = TRUE)$modulus
#   two_neg_ll <- two_neg_ll + sum(e^2)/psi0
#
#   two_neg_ll <- two_neg_ll - Matrix::crossprod(Zte, Matrix::solve(B,
#                                             Matrix::crossprod(Psi0, Zte))) / psi0
#
#   return(-0.5 * as.vector(two_neg_ll))
#
# }

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


chol_solve <- function(U, b)
{
  Matrix::solve(U, Matrix::solve(Matrix::t(U), b))
}


#' @export
res_ll <- function(XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H, Psi0, psi0, score = FALSE,
                   finf = FALSE, lik = TRUE)
{
  # Define dimensions
  n <- length(Y)
  q <- ncol(Psi0)
  r <- ncol(H) / q # Assumes H = [H_1, ... , H_r], where H_j is q by q
  p <- ncol(XtX)

  # Loglikelihood to return
  ll <- NA

  # Score vector to return
  s_psi <- rep(NA, r + 1)

  # Fisher information to return
  I_psi <- matrix(NA, r + 1, r + 1)

  # Pre-compute A = (I_q + Psi0 Z'Z)^{-1} Psi0
  A <- Matrix::crossprod(Psi0, ZtZ) + Matrix::Diagonal(q) # q x q storage

  # Add likelihood term before overwriting
  if(lik) ll <- Matrix::determinant(A, logarithm = TRUE)$modulus

  A <- Matrix::solve(A, Psi0)
  B <- XtZ %*% A # q x q

  # Create XtSiX
  U <- (1 / psi0) * (XtX - Matrix::tcrossprod(B, XtZ)) # p x p, now XtSiX
  U <- Matrix::chol(U) # replace XtSiZ by its Cholesky root

  # Create XtSiY for use in beta_tilde
  beta_tilde <- (1/ psi0) * (XtY - XtZ %*% Matrix::tcrossprod(A, YtZ)) # p x 1
  beta_tilde <- chol_solve(U, beta_tilde)

  # Replace Y by residuals
  Y <- Y - X %*% beta_tilde

  # n x 1 vector for storing \Sigma^{-1}e
  a <- (1 / psi0) * (Y - Z %*% (A %*% Matrix::crossprod(Z, Y))) # n x 1

  if(lik){
    ll <- ll + 2 * sum(log(Matrix::diag(U)))
    ll <- ll + sum(Y * a) + n * log(psi0)
    ll <- -0.5 * ll
  }

  if(score){
    # Stochastic part of restricted score for psi
    s_psi[1] <- 0.5 * sum(a^2)
    v <- as.vector(Matrix::crossprod(Z, a)) # q x 1 vector storage
    s_psi[-1] <- 0.5 * colSums(matrix(as.vector(Matrix::crossprod(v, H)) * v,
                                      nrow = q))
  }
  #############################################################################
  ## NOTHING BELOW SHOULD DEPEND ON Y.
  #############################################################################


  if(finf & score){
    A <- Matrix::tcrossprod(A, ZtZ) # q x q, called M in manuscript
    s_psi[1] <- s_psi[1] - (0.5 / psi0) * n + (0.5 / psi0) * sum(Matrix::diag(A))
    I_psi[1, 1] <- (0.5 / psi0^2) * (n - 2 * sum(Matrix::diag(A)) +
                                       sum(Matrix::t(A) * A))

    E <- Matrix::crossprod(ZtZ, A) # q x q storage

    D <- XtZ %*% A # p x q storage

    XtSiZ <- (1 / psi0) * (XtZ - D) # p x q
    XtSi2Z <- (1 / psi0)^2 * (XtZ - 2 * D + D %*% A) # p x q


    C <- Matrix::tcrossprod(B, XtZ) # p x p storage
    G <- B %*% Matrix::tcrossprod(ZtZ, B) # p x q
    XtSi2X <- (1 / psi0)^2 * (XtX - 2 * C + G) # p x p
    C <- chol_solve(U, XtSi2X) # p x p
    I_psi[1, 1] <- I_psi[1, 1] + 0.5 * sum(C * Matrix::t(C))
    s_psi[1] <- s_psi[1] + 0.5 * sum(Matrix::diag(C))

    XtSi3X <- (1 / psi0^3) * (XtX - 3 * C + 3 * G -
                                D %*% Matrix::tcrossprod(A, D)) # p x p
    # A (q x q), G (p x q) ARE FREE
    A <- (1 / psi0)^2 * (ZtZ - 2 * E + E %*% A) # ZtSi2Z right now
    E <-  (1/ psi0) * (ZtZ - E) # Now holds ZtSiZ
    D <- chol_solve(U, XtSiZ) # p x q
    A <- A - 2 * Matrix::crossprod(D, XtSi2Z) + Matrix::crossprod(XtSiZ, C %*% D)
    I_psi[-1, 1] <- 0.5 * colSums(matrix(Matrix::colSums(as.vector(A) * H), nrow = q))
    I_psi[1, 1] <- I_psi[1, 1] - sum(Matrix::diag(chol_solve(U, XtSi3X)))


    v <- as.vector(E - Matrix::crossprod(XtSiZ, D)) # pq
    s_psi[-1] <- s_psi[-1] - 0.5 * colSums(matrix(Matrix::colSums(v * H), nrow = q))

    H <- Matrix::crossprod(E, H)
    H2 <- Matrix::crossprod(XtSiZ, D %*% H) # Storage can be avoided by
                                            # muliply in loop
    for(ii in 1:r){
      idx1 <- ((ii - 1) * q + 1):(ii * q)
      for(jj in ii:r){
        idx2 <-  ((jj - 1) * q + 1):(jj * q)
        I_psi[jj + 1, ii + 1] <- 0.5 * sum(H[, idx1] * Matrix::t(H[, idx2])) -
          sum(H[, idx1] * Matrix::t(H2[, idx2])) + 0.5 *
          sum(H2[, idx1] * Matrix::t(H2[, idx2]))
      }
    }
  } else if (score){
    A <- Matrix::tcrossprod(A, ZtZ) # q x q, called M in manuscript
    s_psi[1] <- s_psi[1] - (0.5 / psi0) * n + (0.5 / psi0) * sum(Matrix::diag(A))

    D <- XtZ %*% A # p x q storage

    C <- Matrix::tcrossprod(B, XtZ) # p x p storage
    C <- chol_solve(U, (1 / psi0)^2 * (XtX - 2 * C + B %*% Matrix::tcrossprod(ZtZ, B)))
    s_psi[1] <- s_psi[1] + 0.5 * sum(Matrix::diag(C))

    A <- Matrix::crossprod(ZtZ, A)
    A <- (1/ psi0) * (ZtZ - A)

    XtSiZ <- (1 / psi0) * (XtZ - D) # p x q
    D <- chol_solve(U, XtSiZ) # p x q
    v <- as.vector(A - Matrix::crossprod(XtSiZ, D)) # pq
    s_psi[-1] <- s_psi[-1] - 0.5 * colSums(matrix(Matrix::colSums(v * H), nrow = q))
  }

  I_psi <- Matrix::forceSymmetric(I_psi, "L")
  return(list("ll" = ll[1], "score" = s_psi, "finf" = I_psi, "beta" = beta_tilde,
              "I_b_inv_chol" = U))
}
