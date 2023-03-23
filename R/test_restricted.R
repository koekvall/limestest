library(lme4)
library(Matrix)
data(fev1)
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
            data = fev1, REML = TRUE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
Y <- getME(fit, "y")


random_effects <- ranef(fit, condVar = TRUE)
cov_matrices <- attr(random_effects[[1]], "postVar")
H1 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 1, 1, 0), 2, 2))

H3 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 0, 0, 1), 2, 2))
H <- cbind(H1, H2, H3)
Psi_hat <- Matrix::kronecker(Matrix::Diagonal(300), matrix(VarCorr(fit)$id, 2, 2))
psi0_hat <- attr(VarCorr(fit), "sc")^2
beta_hat <- fixef(fit)
e <- Y - X %*% beta_hat
XtX <- Matrix::crossprod(X)
XtY <- Matrix::crossprod(X, Y)
XtZ <- Matrix::crossprod(X, Z)
ZtZ <- Matrix::crossprod(Z)
YtZ <- Matrix::crossprod(Y, Z)


loglik_psi <- function(psi){
  psi0 <- psi[1]
  Psi1 <- matrix(c(psi[2], psi[3], psi[3], psi[4]), 2, 2)
  Psi <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1)
  Psi0 <- Psi / psi0
  res_ll(XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H, Psi0, psi0, score = FALSE,
    finf = FALSE, lik = TRUE)$ll
}

numerical_score <- numDeriv::grad(loglik_psi, c(psi0_hat, Psi_hat[1, 1],
                                                Psi_hat[2, 1], Psi_hat[2, 2]))
analytical_score <-  res_ll(XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H,
                            Psi0 = Psi_hat / psi0_hat, psi0 = psi0_hat,
                            score = TRUE, finf = FALSE, lik = FALSE)$score


cat("The max absolute difference between numerical and analytical score is: ",
    max(abs(numerical_score - analytical_score)), "\n")
cat("The max relative difference between numerical and analytical score is: ",
    max(abs(numerical_score - analytical_score) / numerical_score), "\n")


