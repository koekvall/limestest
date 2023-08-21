do_one_clust_sim <- function(inner_seed, n1, n2, X, Z, b, Psi1, psi0, H, R, XtX, XtZ,
                             ZtZ, Xb, Sigma, make_beta_mat){
  set.seed(inner_seed)
  n <- n1 * n2
  p <- ncol(X)
  Psi0 <- Matrix::kronecker(Matrix::Diagonal(n1), Psi1 / psi0)

  y <- Xb + rnorm(n, sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))
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
                                 "clust" = as.factor(rep(1:n1,
                                                         each = n2)))

  mc_data_frame <- dplyr::bind_cols(mc_data_frame, tidyr::as_tibble(X[, -1, drop = F]))
  fit <- lme4::lmer(out ~ . - clust + (1 + V1|clust), data = mc_data_frame,
                    REML = F)

  beta_null <- make_beta_mat %*% y
  ll_null <- mvtnorm::dmvnorm(x = as.vector(y),
                              mean = as.vector(X %*% beta_null),
                              sigma = as.matrix(Sigma), log = TRUE)

  VC <- as.data.frame(lme4::VarCorr(fit))
  Psi1_mle <- matrix(c(VC$vcov[1], VC$vcov[3], VC$vcov[3], VC$vcov[2]), 2, 2)
  psi0_mle <- VC$vcov[4]
  Psi_mle <- Matrix::kronecker(Matrix::Diagonal(n1), Psi1_mle)
  Sigma_mle <- diag(psi0_mle, n) + Z %*% tcrossprod(Psi_mle, Z)
  ll_mle <- mvtnorm::dmvnorm(x = as.vector(y),
                             mean = as.vector(X %*% stuff_REML$beta),
                             sigma = as.matrix(Sigma_mle), log = TRUE)
  test_stat_lrt <- 2 * (ll_mle - ll_null)

  mles <- c(psi0_mle, Psi1_mle[1, 1], Psi1_mle[1, 2], Psi1_mle[2, 2]) -
    c(psi0, Psi1[1, 1], Psi1[1, 2], Psi1[2, 2])
  test_stat_wald <- as.vector(crossprod(mles, stuff$finf %*% mles))

  c(test_stat, test_stat_REML, test_stat_lrt, test_stat_wald, inner_seed)
}



do_one_cross_sim <- function(inner_seed, n1, n2, X, Z, b, psi1, psi0, H, XtX, XtZ,
                             ZtZ, Xb, Sigma, make_beta_mat){
  set.seed(inner_seed)
  n <- n1 * n2
  p <- ncol(X)
  Psi <- Matrix::Diagonal(n1 + n2, c(rep(psi1[1], n1), rep(psi1[2], n2)))
  Psi0 <- Psi / psi0

  y <- Xb + rnorm(n, sd = sqrt(psi0)) + Z %*% crossprod(sqrt(Psi), rnorm(ncol(Psi)))
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
                                 "clust1" = as.factor(rep(1:n1,
                                                          each = n2)),
                                 "clust2" = as.factor(rep(1:n2, n1)))

  mc_data_frame <- dplyr::bind_cols(mc_data_frame, tidyr::as_tibble(X[, -1, drop = F]))
  fit <- lme4::lmer(out ~ . - clust1 - clust2 + (1|clust1) + (1|clust2), data = mc_data_frame,
                    REML = F)
  VC <- as.data.frame(lme4::VarCorr(fit))
  idx <- c(which(VC$grp == "clust1"), which(VC$grp == "clust2"), which(VC$grp == "Residual"))

  beta_null <- as.vector(make_beta_mat %*% y)
  ll_null <- mvtnorm::dmvnorm(x = as.vector(y),
                              mean = as.vector(X %*% beta_null),
                              sigma = as.matrix(Sigma), log = TRUE)


  psi1_mle <- VC$vcov[idx[1:2]]
  psi0_mle <- VC$vcov[idx[3]]
  Psi_mle <- Matrix::Diagonal(n1 + n2, c(rep(psi1_mle[1], n1), rep(psi1_mle[2], n2)))

  Sigma_mle <- diag(psi0_mle, nrow(X)) + Z %*% tcrossprod(Psi_mle, Z)
  ll_mle <- mvtnorm::dmvnorm(x = as.vector(y),
                             mean = as.vector(X %*% stuff_REML$beta),
                             sigma = as.matrix(Sigma_mle), log = TRUE)
  test_stat_lrt <- 2 * (ll_mle - ll_null)


  mles <- c(psi0_mle - psi0, psi1_mle - psi1)
  test_stat_wald <- as.vector(crossprod(mles, stuff$finf %*% mles))


  c(test_stat, test_stat_REML, test_stat_lrt, test_stat_wald, inner_seed)

}
