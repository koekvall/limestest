library(lme4)

fit_ri <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
fit_rs <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

# ── score_test_lmer return structure ────────────────────────────────────────

test_that("score_test_lmer returns named numeric vector of length 3", {
  res <- score_test_lmer(fit_ri)
  expect_named(res, c("stat", "p_val", "df"))
  expect_length(res, 3L)
  expect_true(is.numeric(res))
})

test_that("p-value is in [0, 1]", {
  res <- score_test_lmer(fit_ri)
  expect_gte(res[["p_val"]], 0)
  expect_lte(res[["p_val"]], 1)
})

test_that("df equals length of test_idx", {
  res <- score_test_lmer(fit_ri, test_idx = 1L)
  expect_equal(res[["df"]], 1)
})

# ── score at MLE should be near zero ────────────────────────────────────────

test_that("score stat near zero when theta_null equals MLE", {
  psi_hat <- reconf:::get_psi_hat_lmer(fit_ri)
  # Test with null = MLE: stat should be ~0, p-value ~1
  res <- score_test_lmer(fit_ri,
                         theta_null = psi_hat,
                         test_idx = 1L,
                         profile = FALSE)
  expect_lt(abs(res[["stat"]]), 0.01)
})

# ── test rejects clearly absent random effect ────────────────────────────────

test_that("random intercept is clearly significant in sleepstudy", {
  # sleepstudy has strong between-subject variability; should reject H0: var=0
  res <- score_test_lmer(fit_ri)
  expect_lt(res[["p_val"]], 0.001)
})

# ── multiple random effects ──────────────────────────────────────────────────

test_that("works with random slope model, testing all RE params", {
  res <- score_test_lmer(fit_rs)
  expect_named(res, c("stat", "p_val", "df"))
  expect_equal(res[["df"]], 3)  # intercept var, covariance, slope var
})
