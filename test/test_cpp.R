library(lme4)
data("Pixel", package="nlme")
mform <- pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog)
fit <- lmer(mform, data = Pixel)
precomp <- limestest:::get_precomp_lmer(fit)
H <- do.call(cbind, precomp$Hlist)
psi_hat <- limestest:::get_psi_hat_lmer(fit)

# Helper functions
Psi_cpp <- limestest:::Psi_from_H_cpp(psi_mr = psi_hat[-length(psi_hat)], H = H)
Psi_R <- limestest:::Psi_from_Hlist(psi_mr = psi_hat[-length(psi_hat)], Hlist = precomp$Hlist)
cat(max(abs(Psi_cpp - Psi_R)), "\n")

# Test likelihoods
Z <- getME(fit, "Z")
X <- getME(fit, "X")
Y <- getME(fit, "y")
b <- getME(fit, "beta")
Y <- Y - X %*% b
r <- length(psi_hat)

# Regular likelihood
loglik_R <- limestest::loglik_psi(Z = Z,
                                   ZtZXe = crossprod(Z, cbind(Z, X, Y)),
                                   e = Y,
                                   H = H,
                                   Psi_r = Psi_cpp / psi_hat[r],
                                   psi_r = psi_hat[r],
                                   get_val = TRUE,
                                   get_score = TRUE,
                                   get_inf = TRUE,
                                   expected = TRUE)
loglik_cpp <- limestest:::loglik_psi_cpp(ZtZ = as(crossprod(Z), "generalMatrix"),
                                         XtZ = as.matrix(crossprod(X, Z)),
                                         Zte = as.vector(crossprod(Z, Y)),
                                         Z = Z,
                                         e = Y,
                                         H = H,
                                         Psi_r = Psi_cpp / psi_hat[r],
                                         psi_r = psi_hat[r],
                                         get_val = TRUE,
                                         get_score = TRUE,
                                         get_inf = TRUE,
                                         expected = TRUE)

abs(loglik_R$value - loglik_cpp$value)

max(abs(loglik_R$score - loglik_cpp$score))

round(abs((as.matrix(loglik_R$inf_mat) - loglik_cpp$inf_mat) / as.matrix(loglik_R$inf_mat)), 1)

# Restricted likelihood
loglik_R <- limestest:::res_ll(XtX = crossprod(X),
                               XtY = crossprod(X, Y),
                               XtZ = crossprod(X, Z),
                               ZtZ = crossprod(Z),
                               crossprod(Y, Z),
                               Y = Y,
                               X = X,
                               Z = Z,
                               H = H,
                               Psi_r = Psi_cpp/psi_hat[r],
                               psi_r = psi_hat[r], get_val = TRUE, get_score = TRUE,
                               get_inf = TRUE)
loglik_cpp <- limestest:::res_ll_cpp(Y = Y,
                                     X = X,
                                     Z = Z,
                                     XtY = as.vector(crossprod(X, Y)),
                                     ZtY = as.vector(crossprod(Z, Y)),
                                     XtX = as.matrix(crossprod(X)),
                                     XtZ = as.matrix(crossprod(X, Z)),
                                     ZtZ = as(crossprod(Z), "generalMatrix"),
                                     H = H,
                                     Psi_r = Psi_cpp/ psi_hat[r],
                                     psi_r = psi_hat[r],
                                     get_val = TRUE,
                                     get_score = TRUE,
                                     get_inf = TRUE)
abs(loglik_R$value - loglik_cpp$value)

loglik_R$score - loglik_cpp$score

round(abs((as.matrix(loglik_R$inf_mat) - loglik_cpp$inf_mat) / as.matrix(loglik_R$inf_mat)), 1)
