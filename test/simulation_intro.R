library(Matrix)
library(doParallel)
set.seed(335)
num_reps <- 1e3
################################################################################
## INDEPENDENT CLUSTERS ########################################################
################################################################################
corr_param <- -0.95
num_cluster <- 20
num_ind <- 5
X <- cbind(1, sample(c(0, 1), num_ind * num_cluster, replace = T))
Z <- Matrix::bdiag() # FIX THIS LINE AND BELOW
b <- rep(0, ncol(X))
Psi1 <- matrix(c(1, corr_param, corr_param, 1), 2, 2)
psi0 <- 1
Psi <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1)
Psi0 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1 / psi0)
H1 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 0, 0, 1), 2, 2))
H <- as(cbind(H1, H2, H3), "sparseMatrix")

Psi1_eig <- eigen(Psi1)
R1 <- tcrossprod(diag(sqrt(Psi1_eig$values), ncol(Psi1)), Psi1_eig$vectors)
R <-  as(Matrix::kronecker(Matrix::Diagonal(num_cluster), R1), "sparseMatrix")

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
cl <- makeCluster(6)
registerDoParallel(cl)
p_vals <- foreach(ii=1:1e3, .combine = rbind, .errorhandling = "remove") %dopar% one_sim_test(ii)
print(colMeans(p_vals < 0.05))
par(mfrow = c(1, 2))
plot(seq(0.01, 0.99, 0.01), quantile(p_vals[, 1], seq(0.01, 0.99, 0.01)), xlab =
       "Uniform quantile", ylab = "Empirical quantile", main = "LL")
abline(a = 0, b = 1, lwd = 2)
plot(seq(0.01, 0.99, 0.01), quantile(p_vals[, 2], seq(0.01, 0.99, 0.01)), xlab =
       "Uniform quantile", ylab = "Empirical quantile", main = "RLL")
abline(a = 0, b = 1, lwd = 2)
stopCluster(cl)


# ################################################################################
# ## CROSSED RANDOM EFFECTS ######################################################
# ################################################################################
# m <- 20 # Number of units
# Z <- as(cbind(Matrix::kronecker(Matrix::diag(1, nrow = m), matrix(1, m, 1)),
#               Matrix::kronecker(matrix(1, m, 1), Matrix::diag(1, nrow = m))), "sparseMatrix")
# X <- matrix(1, nrow = nrow(Z))
# b <- rep(0, ncol(X))
#
# Psi1 <- matrix(c(0, 0, 0, 0.1), 2, 2) #
#
# psi0 <- 1
# Psi <- Matrix::kronecker(Psi1, Matrix::diag(1, m))
# Psi0 <- Matrix::kronecker(Psi1 / psi0, Matrix::diag(1, m))
# H1 <- Matrix::kronecker(diag(c(1, 0), 2), Matrix::diag(1, m))
# H2 <- Matrix::kronecker(diag(c(0, 1), 2), Matrix::diag(1, m))
# H <- as(cbind(H1, H2), "sparseMatrix")
#
# #Psi1_eig <- eigen(Psi1)
# #R1 <- tcrossprod(diag(sqrt(Psi1_eig$values), ncol(Psi1)), Psi1_eig$vectors)
# R <-  as(Matrix::kronecker(diag(sqrt(diag(Psi1))), Matrix::diag(1, m)), "sparseMatrix")
#
# XtX = crossprod(X)
# XtZ = crossprod(X, Z)
# ZtZ = crossprod(Z)
# Xb <- X %*% b
#
# one_sim_test <- function(seed){
#   set.seed(seed)
#   y <- Xb + rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))
#   XtY <- crossprod(X, y)
#   stuff_REML <- limestest::res_ll(XtX =XtX,
#                                   XtY = XtY,
#                                   XtZ = XtZ,
#                                   ZtZ = ZtZ,
#                                   YtZ = crossprod(y, Z),
#                                   Y = y,
#                                   X = X,
#                                   Z = Z,
#                                   H = H,
#                                   Psi0 = Psi0,
#                                   psi0 = psi0,
#                                   score = TRUE,
#                                   finf = TRUE,
#                                   lik = TRUE)
#   e <- y - X %*% stuff_REML$beta
#   stuff <- limestest::loglik_psi(Z = Z,
#                                  ZtZXe = cbind(ZtZ, t(XtZ), crossprod(Z, e)),
#                                  e = e,
#                                  H = H,
#                                  Psi0 = Psi0,
#                                  psi0 = psi0,
#                                  loglik = TRUE,
#                                  score = TRUE,
#                                  finf = TRUE)
#
#   # test_stat <- stuff$score[2]^2 / stuff$finf[2, 2]
#   # test_stat_REML <- stuff_REML$score[2]^2 / stuff_REML$finf[2, 2]
#   #
#   # c(pchisq(test_stat, df = 1, lower = F),
#   #   pchisq(test_stat_REML, df = 1, lower = F))
#
# test_stat <- as.vector(crossprod(stuff$score, solve(stuff$finf, stuff$score)))
# test_stat_REML <- as.vector(crossprod(stuff_REML$score,
#                                       solve(stuff_REML$finf, stuff_REML$score)))
#
# c(pchisq(test_stat, df = 3, lower = F),
#   pchisq(test_stat_REML, df = 3, lower = F))
#
#
# }
# cl <- makeCluster(8)
# registerDoParallel(cl)
# p_vals <- foreach(ii=1:1e3, .combine = rbind, .errorhandling = "remove") %dopar% one_sim_test(ii)
# print(colMeans(p_vals < 0.05))
# par(mfrow = c(1, 2))
# plot(seq(0.01, 0.99, 0.01), quantile(p_vals[, 1], seq(0.01, 0.99, 0.01)), xlab =
#        "Uniform quantile", ylab = "Empirical quantile", main = "LL")
# abline(a = 0, b = 1, lwd = 2)
# plot(seq(0.01, 0.99, 0.01), quantile(p_vals[, 2], seq(0.01, 0.99, 0.01)), xlab =
#        "Uniform quantile", ylab = "Empirical quantile", main = "RLL")
# abline(a = 0, b = 1, lwd = 2)
# stopCluster(cl)
