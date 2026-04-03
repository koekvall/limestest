## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-data----------------------------------------------------------------
library(reconf)
library(lme4)

data(fev1)
dim(fev1)
head(fev1)

## ----visits-per-subject-------------------------------------------------------
range(table(fev1$id))

## ----fit-model----------------------------------------------------------------
fit <- lmer(
  logfev1 ~ age + ht + baseage + baseht + (age | id),
  data = fev1,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

## ----mle----------------------------------------------------------------------
round(fixef(fit), 3)

as.data.frame(VarCorr(fit), order = "lower.tri")[
  , c("grp", "var1", "var2", "vcov")]

## ----ci-all, cache = TRUE-----------------------------------------------------
ci <- ci_all_lmer(fit)
ci

## ----ci-single, cache = TRUE--------------------------------------------------
ci_lmer(fit, test_idx = 3L)   # random slope variance lambda_22

## ----lr-compare, cache = TRUE-------------------------------------------------
confint(fit, method = "profile", parm = "theta_")

