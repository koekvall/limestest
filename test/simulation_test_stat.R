library(Matrix)
library(doParallel)
set.seed(335)
num_reps <- 1e3
################################################################################
## INDEPENDENT CLUSTERS ########################################################
################################################################################

## SET PARAMETERS ##############################################################
num_cluster <- 20
num_ind <- 5
p <- 2
Psi1 <- matrix(c(1, -0.99, -0.99, 1), 2, 2)
psi0 <- 1
################################################################################
n <- num_cluster * num_ind
X <- cbind(1, matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p))
Z <- Matrix::bdiag(lapply(1:num_cluster, function(ii)X[((ii - 1) * num_ind + 1):(ii * num_ind), 1:2]))
b <- rnorm(ncol(X))
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

Sigma <- diag(psi0, nrow(X)) + Z %*% tcrossprod(Psi, Z)
M <- solve(Sigma, X)
G <- solve(crossprod(X, M), t(M))

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
                                 "clust" = as.factor(rep(1:num_cluster,
                                                         each = num_ind)))

  mc_data_frame <- dplyr::bind_cols(mc_data_frame, tidyr::as_tibble(X[, -1]))
  fit <- lme4::lmer(out ~ . - clust + (1 + V1|clust), data = mc_data_frame,
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


  c(test_stat, test_stat_REML, test_stat_lrt)

}
cl <- makeCluster(6)
registerDoParallel(cl)
test_stats <- foreach(ii=1:1e3, .combine = rbind, .errorhandling = "remove") %dopar% one_sim_test(ii)
print(colMeans(test_stats > qchisq(0.95, df = 4)))

ppoints <- seq(0.01, 0.99, 0.01)
par(mfrow = c(1, 3))
plot(qchisq(ppoints, df = 4), quantile(test_stats[, 1], ppoints), xlab =
       "Chi-square (4) quantiles", ylab = "Empirical quantile", main = "LL")
abline(a = 0, b = 1, lwd = 2)
plot(qchisq(ppoints, df = 4), quantile(test_stats[, 2], ppoints), xlab =
       "Chi-square (4) quantiles", ylab = "Empirical quantile", main = "RLL")
abline(a = 0, b = 1, lwd = 2)
plot(qchisq(ppoints, df = 4), quantile(test_stats[, 3], ppoints), xlab =
       "Chi-square (4) quantiles", ylab = "Empirical quantile", main = "LRT")
abline(a = 0, b = 1, lwd = 2)
stopCluster(cl)


################################################################################
## CROSSED RANDOM EFFECTS ######################################################
################################################################################

## SET PARAMETERS ##############################################################
n1 <- 20
n2 <- 10
p <- 10
psi1 <- c(0.1, 1)
psi0 <- 1
################################################################################

Z <- as(cbind(Matrix::kronecker(Matrix::diag(1, nrow = n1), matrix(1, n2, 1)),
           Matrix::kronecker(matrix(1, n1, 1), Matrix::diag(1, nrow = n2))),
        "sparseMatrix")
X <- cbind(1, matrix(runif(n1 * n2 * p, min = -1, max = 1), n1 * n2, p))
b <- rep(0, ncol(X))

Psi <- Matrix::Diagonal(n1 + n2, c(rep(psi1[1], n1), rep(psi1[2], n2)))
Psi0 <- Psi / psi0
H1 <- Matrix::Diagonal(n1 + n2, c(rep(1, n1), rep(0, n2)))
H2 <- Matrix::Diagonal(n1 + n2, c(rep(0, n1), rep(1, n2)))
H <- as(cbind(H1, H2), "sparseMatrix")
R <- sqrt(Psi)

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

  # test_stat <- stuff$score[1]^2 / stuff$finf[1, 1]
  # test_stat_REML <- stuff_REML$score[1]^2 / stuff_REML$finf[1, 1]
  #
  # c(pchisq(test_stat, df = 1, lower = F),
  #   pchisq(test_stat_REML, df = 1, lower = F))

  test_stat <- as.vector(crossprod(stuff$score, solve(stuff$finf, stuff$score)))
  test_stat_REML <- as.vector(crossprod(stuff_REML$score,
                                        solve(stuff_REML$finf, stuff_REML$score)))

  c(pchisq(test_stat, df = 3, lower = F),
    pchisq(test_stat_REML, df = 3, lower = F))


}
cl <- makeCluster(8)
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
