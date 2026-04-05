library(lme4)

fit_ri <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
fit_rs <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

# ── ci_lmer ──────────────────────────────────────────────────────────────────

test_that("ci_lmer returns a 1-row matrix with lower and upper columns", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 1L)
  expect_equal(colnames(ci), c("lower", "upper"))
  expect_true(is.numeric(ci))
})

test_that("ci_lmer: lower < MLE < upper for random intercept variance", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  mle <- as.data.frame(VarCorr(fit_ri), order = "lower.tri")$vcov[1]
  expect_lt(ci[1], mle)
  expect_gt(ci[2], mle)
})

test_that("ci_lmer: lower < upper", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  expect_lt(ci[1], ci[2])
})

test_that("ci_lmer: lower bound is non-negative for variance parameter", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  expect_gte(ci[1], 0)
})

# ── ci_all_lmer ──────────────────────────────────────────────────────────────

test_that("ci_all_lmer returns matrix with correct dimensions", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_ri)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2L)
  expect_equal(colnames(ci), c("lower", "upper"))
  # random intercept model: 1 RE variance + error variance
  expect_equal(nrow(ci), 2L)
})

test_that("ci_all_lmer: all lower < upper", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_rs)
  expect_true(all(ci[, "lower"] < ci[, "upper"]))
})

test_that("ci_all_lmer: all MLEs inside CIs", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_rs)
  vc <- as.data.frame(VarCorr(fit_rs), order = "lower.tri")
  mles <- vc$vcov
  expect_true(all(ci[, "lower"] <= mles & mles <= ci[, "upper"]))
})

test_that("ci_all_lmer: correct number of rows for random slope model", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_rs)
  # intercept var, covariance, slope var, error var = 4 parameters
  expect_equal(nrow(ci), 4L)
})

test_that("ci_all_lmer respects test_idx argument", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_rs, test_idx = 1L)
  expect_equal(nrow(ci), 1L)
})
