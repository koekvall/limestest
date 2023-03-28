set.seed(335)
num_reps <- 1e3
################################################################################
## INDEPENDENT CLUSTERS ########################################################
################################################################################
library(Matrix)
num_cluster <- 10
num_ind <- 12
Z <- Matrix::bdiag(replicate(num_cluster, cbind(1, rnorm(num_ind)),
                             simplify = FALSE))
X <- matrix(rnorm(nrow(Z) * 2), nrow = nrow(Z), ncol = 2)
b <- runif(ncol(X))
Psi1 <- matrix(c(1, 0, 0, 1), 2, 2)
psi0 <- 0.01
Psi <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1)
Psi0 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1 / psi0)
H1 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 0, 0, 1), 2, 2))
H <- cbind(H1, H2, H3)

Psi1_eig <- eigen(Psi1)
R1 <- tcrossprod(diag(sqrt(Psi1_eig$values), ncol(Psi1)), Psi1_eig$vectors)
R <-  Matrix::kronecker(Matrix::Diagonal(num_cluster), R1)

y_outside <- rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))

XtX = crossprod(X)
XtZ = crossprod(X, Z)
ZtZ = crossprod(Z)
Xb <- X %*% b

one_sim_test <- function(seed){
  set.seed(seed)
  y <- Xb + rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))
  XtY <- crossprod(X, y)
  stuff_REML <- limestest::res_ll(XtX =XtX,
                                  XtY = XtY,
                                  XtZ = XtZ,
                                  ZtZ = ZtZ,
                                  YtZ = crossprod(y, Z),
                                  Y = y,
                                  X = X,
                                  Z = Z,
                                  H = H,
                                  Psi0 = Psi0,
                                  psi0 = psi0,
                                  score = TRUE,
                                  finf = TRUE,
                                  lik = TRUE)
  e <- y - X %*% stuff_REML$beta
  stuff <- limestest::loglik_psi(Z = Z,
                                 ZtZXe = cbind(ZtZ, t(XtZ), crossprod(Z, e)),
                                 e = e,
                                 H = H,
                                 Psi0 = Psi0,
                                 psi0 = psi0,
                                 loglik = TRUE,
                                 score = TRUE,
                                 finf = TRUE)

  test_stat <- as.vector(crossprod(stuff$score, solve(stuff$finf, stuff$score)))
  test_stat_REML <- as.vector(crossprod(stuff_REML$score,
                                      solve(stuff_REML$finf, stuff_REML$score)))

  c(pchisq(test_stat, df = 4, lower = F),
    pchisq(test_stat_REML, df = 4, lower = F))

}
library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)
p_vals <- foreach(ii=1:1e3, .combine = rbind, .errorhandling = "remove") %dopar% one_sim_test(ii)
print(colMeans(p_vals < 0.05))
par(mfrow = c(1, 2))
plot(seq(0.01, 0.99, 0.01), quantile(p_vals[, 1], seq(0.01, 0.99, 0.01)), xlab =
       "Uniform quantile", ylab = "Empirical quantile", main = "RLL")
abline(a = 0, b = 1, lwd = 2)
plot(seq(0.01, 0.99, 0.01), quantile(p_vals[, 2], seq(0.01, 0.99, 0.01)), xlab =
       "Uniform quantile", ylab = "Empirical quantile", main = "LL")
abline(a = 0, b = 1, lwd = 2)
stopCluster(cl)
