library(lme4)
data("Pixel", package="nlme")
mform <- pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog)
fit <- lmer(mform, data = Pixel, REML = FALSE)
precomp <- limestest:::get_precomp_lmer(fit)
Hlist <- limestest:::get_Hlist_lmer(fit)
psi_hat <- limestest:::get_psi_hat_lmer(fit)
b_hat <- fixef(fit)

limestest:::maximize_loglik(start_val = c(b_hat, psi_hat),
                            opt_idx = seq(length(psi_hat) + length(b_hat)),
                            Y = getME(fit, "y"),
                            X = getME(fit, "X"),
                            Z = getME(fit, "Z"),
                            Hlist = Hlist,
                            expected = TRUE,
                            REML = FALSE)
