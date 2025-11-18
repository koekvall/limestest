library(tictoc)
library(lme4)
data("Pixel", package="nlme")
p <- 3
r <- 6
###############################################################################
# Unconstrained optimizations
###############################################################################
tic()
fit <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side) + (1|Side:Dog), data = Pixel, REML = FALSE)
toc()

tic()
fit_REML <- update(fit, REML = TRUE)
toc()

precomp <- limestest:::get_precomp_lmer(fit)
Hlist <- limestest:::get_Hlist_lmer(fit)
psi_hat <- limestest:::get_psi_hat_lmer(fit)
b_hat <- fixef(fit)

tic()
fit_our <- limestest:::maximize_loglik(start_val = c(b_hat, psi_hat),
                            opt_idx = seq(p + r),
                            Y = getME(fit, "y"),
                            X = getME(fit, "X"),
                            Z = getME(fit, "Z"),
                            Hlist = Hlist,
                            expected = TRUE,
                            REML = FALSE)
toc()

psi_hat_REML <- limestest:::get_psi_hat_lmer(fit_REML)

tic()
fit_our_REML <- limestest:::maximize_loglik(start_val = c(rep(0, r - 1), 1),
                                       opt_idx = seq_along(psi_hat),
                                       Y = getME(fit, "y"),
                                       X = getME(fit, "X"),
                                       Z = getME(fit, "Z"),
                                       Hlist = Hlist,
                                       expected = TRUE,
                                       REML = TRUE)
toc()

###############################################################################
# Constrained optimizations -- No random effects
###############################################################################

tic()
# Testing no REs
fit_our_null <- limestest:::maximize_loglik(start_val = c(b_hat, rep(0, r - 1), psi_hat[r]),
                                       opt_idx = c(1:p, p + r),
                                       Y = getME(fit, "y"),
                                       X = getME(fit, "X"),
                                       Z = getME(fit, "Z"),
                                       Hlist = Hlist,
                                       expected = TRUE,
                                       REML = FALSE)
toc()
# Should have same coefs as the following
coef(lm(pixel ~ day + I(day^2), data = Pixel))
fit_our_null$arg

tic()
# Testing no REs with REML
fit_our_null_REML <- limestest:::maximize_loglik(start_val = c(rep(0, r - 1), psi_hat[r]),
                                            opt_idx = c(r),
                                            Y = getME(fit, "y"),
                                            X = getME(fit, "X"),
                                            Z = getME(fit, "Z"),
                                            Hlist = Hlist,
                                            expected = TRUE,
                                            REML = TRUE)
toc()
# Should have same estimate as the following
sigma(lm(pixel ~ day + I(day^2), data = Pixel))^2
fit_our_null_REML$arg

###############################################################################
# Constrained optimizations -- Test a random effect variance (Side)
###############################################################################

# Null of no side random effect
fit_our_null <- limestest:::maximize_loglik(start_val = c(b_hat, rep(0, r - 1), psi_hat[r]),
                                            opt_idx = seq(p + r)[-(p + r - 1)],
                                            Y = getME(fit, "y"),
                                            X = getME(fit, "X"),
                                            Z = getME(fit, "Z"),
                                            Hlist = Hlist,
                                            expected = TRUE,
                                            REML = FALSE)
fit_null_lme4 <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1|Side), data = Pixel, REML = FALSE)

# Compare estimates
c(fixef(fit_null_lme4), limestest:::get_psi_hat_lmer(fit_null_lme4))
fit_our_null$arg
