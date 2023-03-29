library(Matrix)
library(doParallel)
test_mc_sigma_crossed <- function(psi_test = c(1, 0, 1), num_reps = 1e3){
  m <- 20 # Number of units
  Z <- as(cbind(Matrix::kronecker(Matrix::diag(1, nrow = m), matrix(1, m, 1)),
                Matrix::kronecker(matrix(1, m, 1), Matrix::diag(1, nrow = m))),
          "sparseMatrix")
  X <- cbind(1, matrix(rnorm(m^2), m^2, 1))
  b <- runif(ncol(X), -1, 1)
  Psi1 <- diag(psi_test[-1], 2)
  psi0 <- psi_test[1]
  Psi <- Matrix::kronecker(Psi1, Matrix::diag(1, m))
  Psi0 <- Matrix::kronecker(Psi1 / psi0, Matrix::diag(1, m))
  H1 <- Matrix::kronecker(diag(c(1, 0), 2), Matrix::diag(1, m))
  H2 <- Matrix::kronecker(diag(c(0, 1), 2), Matrix::diag(1, m))
  H <- as(cbind(H1, H2), "sparseMatrix")
  R <-  as(diag(rep(sqrt(psi_test[-1]), each = m)), "sparseMatrix")
  one_sim_score <- function(seed){
    set.seed(seed)
    y <- rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))

   stuff_REML <- limestest::res_ll(XtX = crossprod(X),
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
                      finf = FALSE,
                      lik = TRUE)
    e <- y - X %*% stuff_REML$beta
    score_LL <- limestest::loglik_psi(Z = Z,
                                   ZtZXe = cbind(crossprod(Z), crossprod(Z, X), crossprod(Z, e)),
                                   e = e,
                                   H = H,
                                   Psi0 = Psi0,
                                   psi0 = psi0,
                                   loglik = TRUE,
                                   score = TRUE,
                                   finf = TRUE)$score
  c(stuff_REML$score, score_LL)
  }
  cl <- makeCluster(8)
  registerDoParallel(cl)
  out <- foreach(ii=1:num_reps, .combine = rbind, .errorhandling = "remove") %dopar% one_sim_score(ii)
  stopCluster(cl)
  MC_inf_REML<- cov(out[, 1:3])
  MC_inf_LL <- cov(out[, 4:6])
  y_outside <- rep(0, nrow(Z))
  AN_inf_REML <- limestest::res_ll(XtX = crossprod(X),
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

  AN_inf_LL <- limestest::loglik_psi(Z = Z,
                        ZtZXe = cbind(crossprod(Z), crossprod(Z, X), crossprod(Z, y_outside)),
                        e = y_outside,
                        H = H,
                        Psi0 = Psi0,
                        psi0 = psi0,
                        loglik = TRUE,
                        score = TRUE,
                        finf = TRUE)$finf
  cat("FOR REML: \n")
  cat("MC information: \n")
  print(MC_inf_REML)
  cat("Analytical information: \n")
  print(as.matrix(AN_inf_REML))

  cat("FOR LL: \n")
  cat("MC information: \n")
  print(MC_inf_LL)
  cat("Analytical information: \n")
  print(as.matrix(AN_inf_LL))
}

test_mc_sigma_crossed(psi_test = c(1, 1, 0))






