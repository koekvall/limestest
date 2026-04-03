# reconf

Score-based confidence intervals and hypothesis tests for variance components in linear mixed models fitted with [lme4](https://github.com/lme4/lme4).

## Overview

Standard methods for inference on variance components — Wald intervals and likelihood ratio tests — have known deficiencies: Wald intervals can extend below zero, and both methods rely on chi-squared approximations that are unreliable in small samples or near the boundary of the parameter space.

`reconf` implements **score-based confidence intervals** that invert a one-dimensional signed score statistic. The method:

- Works directly on the natural variance/covariance parameterisation (no Cholesky reparameterisation needed)
- Handles boundary cases (near-zero variances) without ad-hoc corrections
- Supports both REML and ML estimation
- Uses warm-started outward search for computational efficiency

## Installation

```r
# Install from GitHub
remotes::install_github("koekvall/reconf")
```

## Usage

```r
library(lme4)
library(reconf)

# Fit a linear mixed model
fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)

# 95% score-based CI for all random-effect variance/covariance parameters
ci_all_lmer(fit)

# CI for a single parameter (e.g. the random intercept variance, index 1)
ci_lmer(fit, test_idx = 1)

# Score test: H0 that all random-effect variances are zero
score_test_lmer(fit)

# Score test for a specific parameter against a custom null value
score_test_lmer(fit, test_idx = 1L)
```

The parameter ordering follows `as.data.frame(VarCorr(fit), order = "lower.tri")`, with the residual variance last (excluded from `ci_all_lmer` by default).
