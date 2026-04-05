# Profiling study for ci_lmer / ci_all_lmer
# Run this script interactively and inspect the profvis output in the viewer.

library(reconf)
library(lme4)
library(profvis)

# ---------- Fit model on FEV1 data ----------
data(fev1)
fit <- lmer(
  logfev1 ~ age + ht + baseage + baseht + (age | id),
  data = fev1,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

# ---------- Profile: single CI (random slope variance) ----------
p_single <- profvis({
  ci_lmer(fit, test_idx = 3L)
})
print(p_single)

# ---------- Profile: all CIs ----------
p_all <- profvis({
  ci_all_lmer(fit)
})
print(p_all)

# ---------- Profile: single CI with REML ----------
fit_reml <- lmer(
  logfev1 ~ age + ht + baseage + baseht + (age | id),
  data = fev1,
  REML = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

p_reml <- profvis({
  ci_lmer(fit_reml, test_idx = 3L)
})
print(p_reml)
