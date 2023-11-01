library(lme4)
library(Matrix)
library(limestest)
library(rstiefel)
data(fev1)
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
            data = fev1, REML = FALSE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
y <- getME(fit, "y")

random_effects <- ranef(fit, condVar = TRUE)
cov_matrices <- attr(random_effects[[1]], "postVar")
H1 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 0, 0, 1), 2, 2))
H <- cbind(H1, H2, H3)
Psi1_hat <-  matrix(VarCorr(fit)$id, 2, 2)
Psi_hat <- Matrix::kronecker(Matrix::Diagonal(300), Psi1_hat)
psi0_hat <- attr(VarCorr(fit), "sc")^2
beta_hat <- fixef(fit)
e <- y - X %*% beta_hat
ZtZXe <- Matrix::crossprod(Z, cbind(Z, X, e))

# Loglik for psi^0 at beta = beta_hat
loglik_test <- function(x){
  Psi1_x <- matrix(c(x[2], x[3], x[3], x[4]), 2, 2)
  Psi_x <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_x)
  loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
         Psi0 = Psi_x / x[1], psi0 = x[1], loglik = TRUE,
         score = TRUE, finf = TRUE, expected = FALSE)$ll
}

# Test derivatives at MLE
psi0_test <- psi0_hat
Psi1_test <- Psi_hat[1:2, 1:2]
Psi_test <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_test)
numerical_score <- numDeriv::grad(loglik_test, c(psi0_test, Psi_test[1, 1],
                                                Psi_test[2, 1], Psi_test[2, 2]))
analytical_score <- loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                              Psi0 = Psi_test / psi0_test, psi0 = psi0_test, loglik = F,
                              score = T, finf = F)$score

numerical_hess <- numDeriv::hessian(loglik_test, c(psi0_hat, Psi_test[1, 1],
                                                Psi_test[2, 1], Psi_test[2, 2]))
analytical_hess <- -loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                                 Psi0 = Psi_test/psi0_test, psi0 = psi0_test, loglik = F,
                                 score = T, finf = T, expected = F)$finf

cat("The max absolute difference between numerical and analytical score at MLE is: ",
    max(abs(numerical_score - analytical_score)), "\n")
cat("The max relative difference between numerical and analytical score at MLE is: ",
    max(abs(numerical_score - analytical_score) / numerical_score), "\n")



cat("The max absolute difference between numerical and analytical Hessian at MLE is: ",
    max(abs(numerical_hess - analytical_hess)), "\n")
cat("The max relative difference between numerical and analytical Hessian at MLE is: ",
    max(abs(numerical_hess - analytical_hess) / numerical_hess), "\n")


# Test derivatives at random point
psi0_test <- runif(1, min = psi0_hat / 10, psi0_hat * 10)
ed <- eigen(Psi1_hat)
U <- rustiefel(2, 2)
Psi1_test <- crossprod(U, diag(runif(2, min = min(ed$values) / 10, max = max(ed$values) * 10))) %*% U
Psi_test <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_test)
numerical_score <- numDeriv::grad(loglik_test, c(psi0_test, Psi_test[1, 1],
                                                 Psi_test[2, 1], Psi_test[2, 2]))
analytical_score <- loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                               Psi0 = Psi_test / psi0_test, psi0 = psi0_test, loglik = F,
                               score = T, finf = F)$score

numerical_hess <- numDeriv::hessian(loglik_test, c(psi0_test, Psi_test[1, 1],
                                                   Psi_test[2, 1], Psi_test[2, 2]))
analytical_hess <- -loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                               Psi0 = Psi_test / psi0_test, psi0 = psi0_test, loglik = F,
                               score = T, finf = T, expected = F)$finf

cat("The max absolute difference between numerical and analytical score at random point is: ",
    max(abs(numerical_score - analytical_score)), "\n")
cat("The max relative difference between numerical and analytical score at random point is: ",
    max(abs(numerical_score - analytical_score) / numerical_score), "\n")



cat("The max absolute difference between numerical and analytical Hessian at random point is: ",
    max(abs(numerical_hess - analytical_hess)), "\n")
cat("The max relative difference between numerical and analytical Hessian at random point is: ",
    max(abs(numerical_hess - analytical_hess) / numerical_hess), "\n")


# Simulation for Hessian
# n_reps <- 1e4
# n1 <- 20
# n2 <- 3
# n <- n1 * n2
# p <- 1
# X <- matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p)
# Z <- Matrix::bdiag(replicate(n1, cbind(1, rnorm(n2)),
#                                   simplify = FALSE))
# ZtZ <- crossprod(Z)
# ZtX <- crossprod(Z, X)
# psi0_star <- 0.1
# Psi1_star <- matrix(c(0.5, 0.2, 0.2, 0.5), 2, 2)
# Psi_star <- Matrix::kronecker(Matrix::Diagonal(n1), Psi1_star)
# Psi0_star <- Psi_star / psi0_star
# R_star <- Matrix::kronecker(Matrix::Diagonal(n1), chol(Psi1_star))
# H1 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(1, 0, 0, 0), ncol = 2,nrow =  2))
# H2 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(0, 1, 1, 0), ncol = 2,nrow =  2))
# H3 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(0, 0, 0, 1), ncol = 2,nrow =  2))
# H <- cbind(H1, H2, H3)
#
# hess_list <- list()
# y_list <- list()
# for(ii in 1:n_reps){
#   y <- rnorm(nrow(Z), sd = sqrt(psi0_star)) +
#     Z %*% crossprod(R_star,  rnorm(ncol(R_star))) # beta_star = 0
#   Zte <- crossprod(Z, y)
#   hess_list[[ii]]<- limestest::loglik_psi(Z = Z,
#                         ZtZXe = cbind(ZtZ, ZtX, Zte),
#                         e = y,
#                         H = H,
#                         Psi0 = Psi0_star,
#                         psi0 = psi0_star,
#                         loglik = FALSE,
#                         score = FALSE,
#                         finf = TRUE,
#                         expected = FALSE)$finf
#   y_list[[ii]] <- y
# }
#
# MC_inf <- Reduce("+", hess_list) / length(hess_list)
# analytical_inf <- limestest::loglik_psi(Z = Z,
#                                         ZtZXe = cbind(ZtZ, ZtX, Zte),
#                                         e = y,
#                                         H = H,
#                                         Psi0 = Psi0_star,
#                                         psi0 = psi0_star,
#                                         loglik = FALSE,
#                                         score = FALSE,
#                                         finf = TRUE,
#                                         expected = TRUE)$finf
#
#








