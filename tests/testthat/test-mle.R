library(lme4)

# ── Our optimizer recovers lme4's MLE ────────────────────────────────────────

test_that("maximize_loglik (REML) recovers lme4 estimates", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  r <- length(psi_hat)

  fit_our <- reconf:::maximize_loglik(
    start_val = psi_hat,
    opt_idx   = seq_len(r),
    Y         = getME(fit, "y"),
    X         = getME(fit, "X"),
    Z         = getME(fit, "Z"),
    Hlist     = reconf:::get_Hlist_lmer(fit),
    expected  = TRUE,
    REML      = TRUE
  )

  expect_equal(unname(fit_our$arg), unname(psi_hat), tolerance = 1e-4)
})

test_that("maximize_loglik (ML) recovers lme4 estimates", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)

  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  b_hat   <- fixef(fit)
  r       <- length(psi_hat)
  p       <- length(b_hat)

  fit_our <- reconf:::maximize_loglik(
    start_val = c(b_hat, psi_hat),
    opt_idx   = seq_len(p + r),
    Y         = getME(fit, "y"),
    X         = getME(fit, "X"),
    Z         = getME(fit, "Z"),
    Hlist     = reconf:::get_Hlist_lmer(fit),
    expected  = TRUE,
    REML      = FALSE
  )

  expect_equal(unname(fit_our$arg[(p + 1):(p + r)]), unname(psi_hat), tolerance = 1e-4)
  expect_equal(unname(fit_our$arg[seq_len(p)]), unname(b_hat), tolerance = 1e-4)
})

# ── get_psi_hat_lmer returns the VarCorr vcov vector ─────────────────────────

test_that("get_psi_hat_lmer matches as.data.frame(VarCorr(...))", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  vc <- as.data.frame(VarCorr(fit), order = "lower.tri")$vcov
  expect_equal(psi_hat, vc, tolerance = 1e-12)
})
