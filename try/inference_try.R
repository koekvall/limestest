library(lme4)
data("Pixel", package="nlme")
mform <- pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog)
fit <- lmer(mform, data = Pixel, REML = FALSE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
Y <- getME(fit, "y")

precomp_lmer <- limestest:::get_precomp_lmer(fit)
psi_hat <- limestest:::get_psi_hat_lmer(fit)
psi_start <- psi_hat

fix_idx <- c(1, 3)
fix_vals <- c(200, 0)
psi_start <- psi_hat
psi_start[fix_idx] <- fix_vals
Hlist <- limestest:::get_Hlist_lmer(fit)

fit_null <- limestest:::partial_min_psi(psi_start = psi_start,
                                        opt_idx = seq_len(6)[-fix_idx],
                                        b = NULL,
                                        Y = Y,
                                        X = X,
                                        Z = Z,
                                        Hlist = Hlist,
                                        precomp = precomp_lmer,
                                        REML = FALSE,
                                        expected = TRUE)
psi_tilde <- fit_null$psi_hat

fit_our <- limestest:::partial_min_psi(psi_start = psi_hat,
                                       opt_idx = seq_len(6),
                                       b = NULL,
                                       Y = Y,
                                       X = X,
                                       Z = Z,
                                       Hlist = Hlist,
                                       precomp = precomp_lmer,
                                       REML = TRUE,
                                       expected = TRUE)
psi_hat_our <- fit_our$psi_hat

# These should be in decreasing order
limestest:::loglikelihood(psi = psi_hat_our, b = NULL, Y = Y, X = X, Z = Z, Hlist = Hlist, REML = FALSE)$value
limestest:::loglikelihood(psi = psi_hat, b = NULL, Y = Y, X = X, Z = Z, Hlist = Hlist, REML = FALSE)$value
limestest:::loglikelihood(psi = psi_tilde, b = NULL, Y = Y, X = X, Z = Z, Hlist = Hlist, REML = FALSE)$value
limestest:::loglikelihood(psi = psi_start, b = NULL, Y = Y, X = X, Z = Z, Hlist = Hlist, REML = FALSE)$value

# Try tests
limestest:::score_test_lmer(fit, test_idx = 3)


###############################################################################
# Test with example that failed for trust
###############################################################################
data(fev1)
fit <- lme4::lmer(exp(logfev1) ~ age + ht + baseage + baseht + (1|id)  + (0 + age|id),
            data = fev1, REML = FALSE)
psi_hat <- limestest:::get_psi_hat_lmer(fit)

# This fails
limestest:::score_test_lmer(fit, test_idx = 1)

# Search for error; is starting point valid?
limestest:::test_all_lmer(fit)


