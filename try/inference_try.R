library(lme4)
data("Pixel", package="nlme")
mform <- pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog)
fit <- lmer(mform, data = Pixel)
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

fit_null <- limestest:::partial_min(opt_idx = seq_len(6)[-fix_idx],
                                    precomp = precomp_lmer,
                        psi_start = psi_start, REML = TRUE,
                        expected = TRUE)
psi_tilde <- fit_null$psi_hat

fit_our <- limestest:::partial_min(opt_idx = seq_len(6), precomp = precomp_lmer,
                                    psi_start = psi_hat, REML = TRUE,
                                   expected = TRUE)
psi_hat_our <- fit_our$psi_hat

# These should be in decreasing order
limestest:::loglikelihood(psi = psi_hat_our, b = NULL, precomp = precomp_lmer, REML = TRUE)$value
limestest:::loglikelihood(psi = psi_hat, b = NULL, precomp = precomp_lmer, REML = TRUE)$value
limestest:::loglikelihood(psi = psi_tilde, b = NULL, precomp = precomp_lmer, REML = TRUE)$value
limestest:::loglikelihood(psi = psi_start, b = NULL, precomp = precomp_lmer, REML = TRUE)$value

# Try tests
limestest:::score_test_lmer(fit, joint = FALSE)
