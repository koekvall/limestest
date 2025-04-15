library(lme4)
data("Pixel", package="nlme")
mform <- pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog)
fit <- lmer(mform, data = Pixel)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
Y <- getME(fit, "y")

precomp_lmer <- list(XtX = crossprod(X),
                XtY = crossprod(X, Y),
                XtZ = crossprod(X, Z),
                ZtZ = crossprod(Z),
                YtZ = crossprod(Y, Z),
                Y = Y,
                X = X,
                Z = Z,
                Hlist = limestest:::get_Hlist(fit))
psi_hat <- as.data.frame(VarCorr(fit), order = "lower.tri")$vcov
psi_start <- psi_hat

fix_idx <- c(1, 3)
fix_vals <- c(200, 0)
psi_start <- psi_hat
psi_start[fix_idx] <- fix_vals

profvis::profvis(fit_null <- limestest:::partial_min(opt_idx = seq_len(6)[-fix_idx], precomp = precomp_lmer,
                        psi_start = psi_start, REML = TRUE,
                        expected = TRUE))
psi_tilde <- fit_null$psihat

fit_our <- limestest:::partial_min(opt_idx = seq_len(6), precomp = precomp_lmer,
                                    psi_start = psi_hat, REML = TRUE,
                                   expected = TRUE)
psi_hat_our <- fit_our$psihat

# These should be in decreasing order
limestest:::loglikelihood(psi = psi_hat_our, b = NULL, precomp = precomp_lmer, REML = TRUE)$value
limestest:::loglikelihood(psi = psi_hat, b = NULL, precomp = precomp_lmer, REML = TRUE)$value
limestest:::loglikelihood(psi = psi_tilde, b = NULL, precomp = precomp_lmer, REML = TRUE)$value
limestest:::loglikelihood(psi = psi_start, b = NULL, precomp = precomp_lmer, REML = TRUE)$value

