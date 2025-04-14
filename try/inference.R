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
psi_null <- as.data.frame(VarCorr(fit), order = "lower.tri")$vcov
psi_null[3] <- 0
limestest:::partial_min(opt_idx = seq_len(6)[-3], precomp = precomp_lmer,
                        psi_start = psi_null, REML = TRUE, expected = TRUE)
