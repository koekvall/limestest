library(lme4)
library(Matrix)

# ── R and C++ implementations agree ──────────────────────────────────────────

test_that("Psi_from_H_cpp and Psi_from_Hlist agree", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Hlist   <- reconf:::get_Hlist_lmer(fit)
  H       <- do.call(cbind, Hlist)
  r       <- length(psi_hat)

  Psi_cpp <- reconf:::Psi_from_H_cpp(psi_mr = psi_hat[-r], H = H)
  Psi_R   <- reconf:::Psi_from_Hlist(psi_mr = psi_hat[-r], Hlist = Hlist)

  expect_equal(as.matrix(Psi_cpp), as.matrix(Psi_R), tolerance = 1e-12)
})

test_that("C++ and R log-likelihood values agree at MLE (ML)", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  b_hat   <- getME(fit, "beta")
  Hlist   <- reconf:::get_Hlist_lmer(fit)

  ll_cpp <- reconf:::loglikelihood(
    psi = psi_hat, b = b_hat,
    Y = getME(fit, "y"), X = getME(fit, "X"), Z = getME(fit, "Z"),
    Hlist = Hlist, REML = FALSE,
    get_val = TRUE, get_score = FALSE, get_inf = FALSE
  )$value

  # R implementation
  H     <- do.call(cbind, Hlist)
  r     <- length(psi_hat)
  e     <- getME(fit, "y") - getME(fit, "X") %*% b_hat
  Z     <- getME(fit, "Z")
  Psi_r <- reconf:::Psi_from_H_cpp(psi_mr = psi_hat[-r], H = H) / psi_hat[r]

  ll_R <- reconf:::loglik_psi(
    Z = Z,
    ZtZXe = crossprod(Z, cbind(Z, getME(fit, "X"), e)),
    e = e, H = H, Psi_r = Psi_r, psi_r = psi_hat[r],
    get_val = TRUE, get_score = FALSE, get_inf = FALSE
  )$value

  expect_equal(ll_cpp, ll_R, tolerance = 1e-6)
})

test_that("C++ and R restricted log-likelihood values agree at MLE (REML)", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Hlist   <- reconf:::get_Hlist_lmer(fit)
  H       <- do.call(cbind, Hlist)
  r       <- length(psi_hat)
  Y       <- getME(fit, "y")
  X       <- getME(fit, "X")
  Z       <- getME(fit, "Z")
  Psi_r   <- reconf:::Psi_from_H_cpp(psi_mr = psi_hat[-r], H = H) / psi_hat[r]

  ll_cpp <- reconf:::loglikelihood(
    psi = psi_hat,
    Y = Y, X = X, Z = Z, Hlist = Hlist, REML = TRUE,
    get_val = TRUE, get_score = FALSE, get_inf = FALSE
  )$value

  ll_R <- reconf:::res_ll(
    XtX = crossprod(X), XtY = crossprod(X, Y),
    XtZ = crossprod(X, Z), ZtZ = crossprod(Z),
    YtZ = crossprod(Y, Z),
    Y = Y, X = X, Z = Z, H = H,
    Psi_r = Psi_r, psi_r = psi_hat[r],
    get_val = TRUE, get_score = FALSE, get_inf = FALSE
  )$value

  expect_equal(ll_cpp, ll_R, tolerance = 1e-6)
})

test_that("C++ and R score vectors agree at MLE (REML)", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Hlist   <- reconf:::get_Hlist_lmer(fit)
  H       <- do.call(cbind, Hlist)
  r       <- length(psi_hat)
  Y       <- getME(fit, "y")
  X       <- getME(fit, "X")
  Z       <- getME(fit, "Z")
  Psi_r   <- reconf:::Psi_from_H_cpp(psi_mr = psi_hat[-r], H = H) / psi_hat[r]

  score_cpp <- reconf:::loglikelihood(
    psi = psi_hat,
    Y = Y, X = X, Z = Z, Hlist = Hlist, REML = TRUE,
    get_val = FALSE, get_score = TRUE, get_inf = FALSE
  )$score

  score_R <- reconf:::res_ll(
    XtX = crossprod(X), XtY = crossprod(X, Y),
    XtZ = crossprod(X, Z), ZtZ = crossprod(Z),
    YtZ = crossprod(Y, Z),
    Y = Y, X = X, Z = Z, H = H,
    Psi_r = Psi_r, psi_r = psi_hat[r],
    get_val = FALSE, get_score = TRUE, get_inf = FALSE
  )$score

  expect_equal(score_cpp, score_R, tolerance = 1e-6)
})
