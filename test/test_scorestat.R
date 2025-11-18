library(tictoc)
library(lme4)
data("Pixel", package="nlme")
p <- 3
r <- 6
###############################################################################
# Test covariance between slope and intercept zero
###############################################################################
fit_lme4 <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side) + (1|Side:Dog), data = Pixel, REML = TRUE)
psi_hat <- limestest:::get_psi_hat_lmer(fit_lme4)
fit_ours <- limestest:::maximize_loglik(start_val = psi_hat,
                                        opt_idx = seq_along(psi_start),
                                        Y = getME(fit_lme4, "y"),
                                        X = getME(fit_lme4, "X"),
                                        Z = getME(fit_lme4, "Z"),
                                        Hlist = limestest:::get_Hlist_lmer(fit_lme4),
                                        expected = TRUE,
                                        REML = TRUE)
# Get constrained estimates of nuisance parameters
psi_start <- psi_hat
psi_start[3] <- 0
fit_null_ours <- limestest:::maximize_loglik(start_val = psi_start,
                                        opt_idx = seq_along(psi_start)[-3],
                                        Y = getME(fit_lme4, "y"),
                                        X = getME(fit_lme4, "X"),
                                        Z = getME(fit_lme4, "Z"),
                                        Hlist = limestest:::get_Hlist_lmer(fit_lme4),
                                        expected = TRUE,
                                        REML = TRUE)

# Score test
stat <- limestest:::score_stat(theta = fit_null_ours$arg,
                       test_idx = 3,
                       Y = getME(fit_lme4, "y"),
                       X = getME(fit_lme4, "X"),
                       Z = getME(fit_lme4, "Z"),
                       Hlist = limestest:::get_Hlist_lmer(fit_lme4),
                       expected = TRUE,
                       REML = TRUE,
                       signed = FALSE)
pchisq(stat, df = 1, lower = FALSE)

# RLRT
stat_LRT <- 2 * (fit_ours$value - fit_null_ours$value)
pchisq(stat_LRT, df = 1, lower = FALSE)

###############################################################################
# Get confidence interval for covariance between slope and intercept
###############################################################################
