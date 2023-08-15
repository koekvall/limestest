out_dir <- "~/GitHub/lmm-crit-suppl/"
PDF <- TRUE
library(Matrix)
library(doParallel)
library(tidyverse)
library(lme4)
library(mvtnorm)
set.seed(335)
num_reps <- 1e4
################################################################################
## INDEPENDENT CLUSTERS ########################################################
################################################################################
corr_param <- 0
var_params <- c(0.01, 0.01)
num_cluster <- 20
num_ind <- 5
X <- cbind(1, sample(c(0, 1), num_ind * num_cluster, replace = T))
Z <- Matrix::bdiag(lapply(1:num_cluster, function(ii)X[((ii - 1) * num_ind + 1):(ii * num_ind), ]))
b <- rep(0, ncol(X))
Psi1 <- matrix(c(var_params[1], sqrt(prod(var_params)) * corr_param,
                 sqrt(prod(var_params)) * corr_param, var_params[2]), 2, 2)
psi0 <- 1
Psi <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1)
Sigma <- diag(psi0, nrow(X)) + Z %*% tcrossprod(Psi, Z)
M <- solve(Sigma, X)
G <- solve(crossprod(X, M), t(M))

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


  mc_data_frame <- tidyr::tibble("out" = as.vector(y),
                             "trt" = c(X[, 2]),
                             "clust" = as.factor(rep(1:num_cluster,
                                                     each = num_ind)))
  fit <- lme4::lmer(out ~ trt + (1 + trt|clust), data = mc_data_frame,
                    REML = F)

  beta_null <- G %*% y
  ll_null <- mvtnorm::dmvnorm(x = as.vector(y),
                              mean = as.vector(X %*% beta_null),
                              sigma = as.matrix(Sigma), log = TRUE)

  beta_mle <- lme4::fixef(fit)
  VC <- as.data.frame(lme4::VarCorr(fit))
  Psi1_mle <- matrix(c(VC$vcov[1], VC$vcov[3], VC$vcov[3], VC$vcov[2]), 2, 2)
  psi0_mle <- VC$vcov[4]
  Psi_mle <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1_mle)
  Sigma_mle <- diag(psi0_mle, nrow(X)) + Z %*% tcrossprod(Psi_mle, Z)
  ll_mle <- mvtnorm::dmvnorm(x = as.vector(y),
                             mean = as.vector(X %*% beta_mle),
                             sigma = as.matrix(Sigma_mle), log = TRUE)
  test_stat_lrt <- 2 * (ll_mle - ll_null)

  c(test_stat,
   test_stat_REML,
    test_stat_lrt)

}
cl <- makeCluster(6)
registerDoParallel(cl)
test_stats <- foreach(ii=1:1e3, .combine = rbind, .errorhandling = "remove") %dopar% one_sim_test(ii)

if(PDF) pdf(paste(out_dir, "fig_loglik_nocorr.pdf", sep = ""), width = 6, height = 4)
par(mfrow = c(1, 2))
ppoints <- seq(0.01, 0.99, 0.01)
plot(qchisq(ppoints, df = 4), quantile(test_stats[, 1], ppoints), xlab =
       "", ylab = "", main = "Score",
     xlim = c(0, 12), ylim = c(0, 12))
abline(a = 0, b = 1, lwd = 2)
# plot(seq(0.01, 0.99, 0.01), quantile(p_vals[, 2], seq(0.01, 0.99, 0.01)), xlab =
#        "Uniform quantile", ylab = "Empirical quantile", main = "RLL")
# abline(a = 0, b = 1, lwd = 2)
plot(qchisq(ppoints, df = 4), quantile(test_stats[, 3], ppoints), xlab =
       "", ylab = "", main = "LRT", xlim = c(0, 12), ylim = c(0, 12))
abline(a = 0, b = 1, lwd = 2)
abline(a = 0, b = 1, lwd = 2)
stopCluster(cl)
if(PDF) dev.off()

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
