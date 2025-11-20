# devtools::install_github("koekvall/limestest@develop")
library(tictoc)
library(lme4)
data("Pixel", package="nlme")
p <- 3
r <- 6
###############################################################################
# Test covariance between slope and intercept zero
###############################################################################
fit_lme4 <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side) +
                   (1|Side:Dog), data = Pixel, REML = TRUE)
psi_hat <- limestest:::get_psi_hat_lmer(fit_lme4)
psi_start <- psi_hat
fit_ours <- limestest::maximize_loglik(start_val = psi_hat,
                                        opt_idx = seq_along(psi_start),
                                        Y = getME(fit_lme4, "y"),
                                        X = getME(fit_lme4, "X"),
                                        Z = getME(fit_lme4, "Z"),
                                        Hlist = limestest:::get_Hlist_lmer(fit_lme4),
                                        expected = TRUE,
                                        REML = TRUE)
# Get constrained estimates of nuisance parameters
psi_start[3] <- 0
fit_null_ours <- limestest::maximize_loglik(start_val = psi_start,
                                        opt_idx = seq_along(psi_start)[-3],
                                        Y = getME(fit_lme4, "y"),
                                        X = getME(fit_lme4, "X"),
                                        Z = getME(fit_lme4, "Z"),
                                        Hlist = limestest:::get_Hlist_lmer(fit_lme4),
                                        expected = TRUE,
                                        REML = TRUE)

# Score test
stat <- limestest::score_stat(theta = fit_null_ours$arg,
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
# Check confidence interval fixed intercept (sanity check)
###############################################################################
fit_lme4_ml <- update(fit_lme4, REML = F)
theta_hat <- c(fixef(fit_lme4_ml), limestest:::get_psi_hat_lmer(fit_lme4_ml))
theta_test <- theta_hat
theta_test[1] <- 1073
stat_values <- limestest::score_nuisance(theta_start = theta_test,
                                         test_idx = 1,
                                         max_radius = 100,
                                         num_points = 100,
                                         Y = getME(fit_lme4_ml, "y"),
                                         X = getME(fit_lme4_ml, "X"),
                                         Z = getME(fit_lme4_ml, "Z"),
                                         Hlist = limestest:::get_Hlist_lmer(fit_lme4_ml),
                                         REML = FALSE, iterlim = 1000, signed = FALSE,
                                         # Changing known changes behavior.
                                         # Try known_idx = NULL (all treated unknown)
                                         # and known_idx = 2:9 (all treated known)
                                         known_idx = 2:9)
plot(x = as.numeric(attr(stat_values, "null_values")), abs(stat_values), xlab = "Parameter",
     ylab = "Absolute value of test statistic")
abline(h = 1.96^2)

## TROUBESHOOTING BASED ON ABOVE PLOT
theta_test <- theta_hat
theta_test[1] <- 1050
fit_null_ours <- limestest:::maximize_loglik(start_val = theta_test,
                                        opt_idx = seq_along(theta_test)[-1],
                                        Y = getME(fit_lme4_ml, "y"),
                                        X = getME(fit_lme4_ml, "X"),
                                        Z = getME(fit_lme4_ml, "Z"),
                                        Hlist = limestest:::get_Hlist_lmer(fit_lme4_ml),
                                        expected = TRUE,
                                        REML = FALSE, iterlim = 1000)

# This stat is greater than the value in the plot, indicating issue in
# score_nuisance
stat <- limestest:::score_stat(theta = fit_null_ours$arg,
                               test_idx = 1,
                               Y = getME(fit_lme4_ml, "y"),
                               X = getME(fit_lme4_ml, "X"),
                               Z = getME(fit_lme4_ml, "Z"),
                               Hlist = limestest:::get_Hlist_lmer(fit_lme4_ml),
                               expected = TRUE,
                               REML = FALSE,
                               signed = FALSE)



###############################################################################
# Check confidence interval for covariance between slope and intercept
###############################################################################
psi_null <- psi_hat
psi_null[3] <- 0
stat_values <- limestest::score_nuisance(theta_start = psi_null,
                          test_idx = 3,
                          max_radius = 200,
                          num_points = 1e2,
                          Y = getME(fit_lme4, "y"),
                          X = getME(fit_lme4, "X"),
                          Z = getME(fit_lme4, "Z"),
                          Hlist = limestest:::get_Hlist_lmer(fit_lme4),
                          iterlim = 1e3)
plot(x = as.numeric(names(stat_values)), abs(stat_values), xlab = "Parameter",
     ylab = "Absolute value of test statistic")
abline(h = 1.96)

###############################################################################
# Check confidence interval for variance of (1|Side)
###############################################################################
psi_null <- psi_hat
psi_null[5] <- 0
stat_values <- limestest::score_nuisance(theta_start = psi_null,
                                         test_idx = 5,
                                         max_radius = c(0, 1000),
                                         num_points = 1e3,
                                         Y = getME(fit_lme4, "y"),
                                         X = getME(fit_lme4, "X"),
                                         Z = getME(fit_lme4, "Z"),
                                         Hlist = limestest:::get_Hlist_lmer(fit_lme4),
                                         efficient = FALSE)
plot(x = as.numeric(names(stat_values)), abs(stat_values), xlab = "Parameter",
     ylab = "Absolute value of test statistic")
abline(h = 1.96)

