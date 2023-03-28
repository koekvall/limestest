library(Matrix)
test_mc_sigma_crossed <- function(psi_test){
  m <- 20 # Number of units
  Z <- as(cbind(Matrix::kronecker(Matrix::diag(1, nrow = m), matrix(1, m, 1)),
                Matrix::kronecker(matrix(1, m, 1), Matrix::diag(1, nrow = m))),
          "sparseMatrix")
  X <- cbind(1, matrix(rnorm(m^2), m^2, 1))
  b <- rep(0, ncol(X))
  Psi1 <- diag(psi_test[-1], 2)
  psi0 <- psi_test[1]
  Psi <- Matrix::kronecker(Psi1, Matrix::diag(1, m))
  Psi0 <- Matrix::kronecker(Psi1 / psi0, Matrix::diag(1, m))
  H1 <- Matrix::kronecker(diag(c(1, 0), 2), Matrix::diag(1, m))
  H2 <- Matrix::kronecker(diag(c(0, 1), 2), Matrix::diag(1, m))
  H <- as(cbind(H1, H2), "sparseMatrix")

  #Psi1_eig <- eigen(Psi1)
  #R1 <- tcrossprod(diag(sqrt(Psi1_eig$values), ncol(Psi1)), Psi1_eig$vectors)
  R <-  as(diag(rep(sqrt(psi_test[-1]), each = m)), "sparseMatrix")

  num_reps <- 1e3
  one_sim_score <- function(seed){
    set.seed(seed)
    y <- rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))

    limestest::res_ll(XtX = crossprod(X),
                      XtY = crossprod(X, y),
                      XtZ = crossprod(X, Z),
                      ZtZ = crossprod(Z),
                      YtZ = crossprod(y, Z),
                      Y = y,
                      X = X,
                      Z = Z,
                      H = H,
                      Psi0 = Psi0,
                      psi0 = psi0,
                      score = TRUE,
                      finf = TRUE,
                      lik = TRUE)$score

  }

}



out <- sapply(1:num_reps, one_sim_score)
MC_inf<- cov(t(out))
AN_inf <- limestest::res_ll(XtX = crossprod(X),
                           XtY = crossprod(X, y_outside),
                           XtZ = crossprod(X, Z),
                           ZtZ = crossprod(Z),
                           YtZ = crossprod(y_outside, Z),
                           Y = y_outside,
                           X = X,
                           Z = Z,
                           H = H,
                           Psi0 = Psi0,
                           psi0 = psi0,
                           score = TRUE,
                           finf = TRUE,
                           lik = TRUE)$finf

MC_inf
AN_inf


