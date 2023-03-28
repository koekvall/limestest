library(lme4)
library(Matrix)
data(fev1)

test_real_data <- function(psi_test = c(1, 1, 0, 0)){
  fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
              data = fev1, REML = TRUE)
  X <- getME(fit, "X")
  Z <- getME(fit, "Z")
  Y <- getME(fit, "y")
  H1 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(1, 0, 0, 0), 2, 2))
  H2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 1, 1, 0), 2, 2))
  H3 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 0, 0, 1), 2, 2))
  H <- as(cbind(H1, H2, H3), "sparseMatrix")
  Psi_hat <- Matrix::kronecker(Matrix::Diagonal(300), matrix(VarCorr(fit)$id, 2, 2))
  psi0_hat <- attr(VarCorr(fit), "sc")^2
  psi_hat <- c(psi0_hat, Psi_hat[1, 1], Psi_hat[2, 1], Psi_hat[2, 2])
  beta_hat <- fixef(fit)
  e <- Y - X %*% beta_hat
  XtX <- Matrix::crossprod(X)
  XtY <- Matrix::crossprod(X, Y)
  XtZ <- Matrix::crossprod(X, Z)
  ZtZ <- Matrix::crossprod(Z)
  YtZ <- Matrix::crossprod(Y, Z)

  loglik_test <- function(x){
    psi0_arg <- x[1]
    Psi1_arg <- matrix(c(x[2], x[3], x[3], x[4]), 2, 2)
    Psi0_arg <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_arg / psi0_arg)
    res_ll(XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H, Psi0_arg, psi0_arg, lik = TRUE, score = FALSE,
           finf = FALSE)$ll
  }

  score_test <- function(x){
    psi0_arg <- x[1]
    Psi1_arg <- matrix(c(x[2], x[3], x[3], x[4]), 2, 2)
    Psi0_arg <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_arg / psi0_arg)
    res_ll(XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H, Psi0_arg, psi0_arg, lik = FALSE,
           score = TRUE, finf = FALSE)$score
  }

  numerical_score <- numDeriv::grad(loglik_test, psi_hat)
  analytical_score <- score_test(psi_hat)


  cat("Max abs. diff. numerical and analytical restricted score at psi_hat: ",
      max(abs(numerical_score - analytical_score)), "\n")
  cat("Max rel. diff. numerical and analytical restricted score at psi_hat: ",
      max(abs(numerical_score - analytical_score) / numerical_score), "\n")

  numerical_score <- numDeriv::grad(loglik_test, psi_test)
  analytical_score <- score_test(psi_test)

  cat("Max abs. diff. numerical and analytical restricted score at psi_test: ",
      max(abs(numerical_score - analytical_score)), "\n")
  cat("Max rel. diff. numerical and analytical restricted score at psi_test: ",
      max(abs(numerical_score - analytical_score) / numerical_score), "\n")
}

cat("Test on real data: \n")
test_real_data()


test_crossed_data <- function(psi_test = c(1, 1, 0)){
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

  XtX = crossprod(X)
  XtZ = crossprod(X, Z)
  ZtZ = crossprod(Z)
  y <-  X %*% b + rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))
  YtZ <- crossprod(y, Z)
  XtY = crossprod(X, y)
  loglik_test <- function(x){
    psi0_arg <- x[1]
    Psi1_arg <- matrix(c(x[2],0, 0, x[3]), 2, 2)
    Psi0_arg <-  Matrix::kronecker(Psi1_arg / psi0_arg, Matrix::diag(1, m))
    res_ll(XtX = XtX,
           XtY = XtY,
           XtZ = XtZ,
           ZtZ = ZtZ,
           YtZ = YtZ,
           Y = y,
           X = X,
           Z = Z,
           H = H,
           Psi0 = Psi0_arg,
           psi0 = psi0_arg,
           lik = TRUE, score = FALSE,
           finf = FALSE)$ll
  }

  score_test <- function(x){
    psi0_arg <- x[1]
    Psi1_arg <- matrix(c(x[2],0, 0, x[3]), 2, 2)
    Psi0_arg <-  Matrix::kronecker(Psi1_arg / psi0_arg, Matrix::diag(1, m))
    res_ll(XtX = XtX,
           XtY = XtY,
           XtZ = XtZ,
           ZtZ = ZtZ,
           YtZ = YtZ,
           Y = y,
           X = X,
           Z = Z,
           H = H,
           Psi0 = Psi0_arg,
           psi0 = psi0_arg,
           lik = FALSE, score = TRUE,
           finf = FALSE)$score
  }


  numerical_score <- numDeriv::grad(loglik_test, psi_test)
  analytical_score <- score_test(psi_test)

  cat("Max abs. diff. numerical and analytical restricted score at psi_test: ",
      max(abs(numerical_score - analytical_score)), "\n")
  cat("Max rel. diff. numerical and analytical restricted score at psi_test: ",
      max(abs(numerical_score - analytical_score) / numerical_score), "\n")
}
cat("Test on simulated crossed random effects data: \n")
test_crossed_data()
