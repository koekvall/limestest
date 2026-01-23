library(lme4)
data("Pixel", package="nlme")
mform <- pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog)
fit <- lmer(mform, data = Pixel)
precomp <- limestest:::get_precomp_lmer(fit)
Hlist <- limestest:::get_Hlist_lmer(fit)
H <- do.call(cbind, Hlist)
psi_hat <- limestest:::get_psi_hat_lmer(fit)
Z <- getME(fit, "Z")
X <- getME(fit, "X")
p <- ncol(X)
Y <- getME(fit, "y")
b <- getME(fit, "beta")
Y <- Y - X %*% b
r <- length(psi_hat)

# Test helper functions
Psi_cpp <- limestest:::Psi_from_H_cpp(psi_mr = psi_hat[-r], H = H)
Psi_R <- limestest:::Psi_from_Hlist(psi_mr = psi_hat[-r], Hlist = Hlist)

cat("Max difference in R and Cpp Psi: ", max(abs(Psi_cpp - Psi_R)), "\n")

# Test likelihoods
# Regular likelihood
loglik_R <- limestest:::loglik_psi(Z = Z,
                                   ZtZXe = crossprod(Z, cbind(Z, X, Y)),
                                   e = Y,
                                   H = H,
                                   Psi_r = Psi_cpp / psi_hat[r],
                                   psi_r = psi_hat[r],
                                   get_val = TRUE,
                                   get_score = TRUE,
                                   get_inf = TRUE,
                                   expected = TRUE)

loglik_cpp <- limestest:::loglik(Psi_r = Psi_cpp / psi_hat[r],
                                 psi_r = psi_hat[r],
                                 H = H,
                                 e = Y,
                                 X = X,
                                 Z = Z,
                                 XtX = crossprod(X),
                                 XtZ = as.matrix(crossprod(X, Z)),
                                 ZtZ = as(crossprod(Z), "generalMatrix"),
                                 get_inf = TRUE,
                                 expected = TRUE)

(as.matrix(loglik_R$inf_mat) - loglik_cpp$inf_mat[4:9, 4:9]) / as.matrix(loglik_R$inf_mat)

cat("Difference in R and Cpp loglik: ", abs(loglik_R$value - loglik_cpp$value) , "\n")
cat("Max difference in R and Cpp score: ", max(abs(loglik_R$score - loglik_cpp$score[-c(1:p)])) , "\n")
cat("Max difference in R and Cpp information: ", max(abs(loglik_R$inf_mat - loglik_cpp$inf_mat[-c(1:p), -c(1:p)])) , "\n")

loglik_R <- limestest:::loglik_psi(Z = Z,
                                  ZtZXe = crossprod(Z, cbind(Z, X, Y)),
                                  e = Y,
                                  H = H,
                                  Psi_r = Psi_cpp / psi_hat[r],
                                  psi_r = psi_hat[r],
                                  get_val = TRUE,
                                  get_score = TRUE,
                                  get_inf = TRUE,
                                  expected = FALSE)

loglik_cpp <- loglik_cpp <- limestest:::loglik(Psi_r = Psi_cpp / psi_hat[r],
                                               psi_r = psi_hat[r],
                                               H = H,
                                               e = Y,
                                               X = X,
                                               Z = Z,
                                               XtX = crossprod(X),
                                               XtZ = as.matrix(crossprod(X, Z)),
                                               ZtZ = as(crossprod(Z), "generalMatrix"),
                                               get_inf = TRUE,
                                               expected = FALSE)

cat("Max difference in R and Cpp Hessian: ", max(abs(loglik_R$inf_mat - loglik_cpp$inf_mat[-c(1:p), -c(1:p)])) , "\n")

# Restricted likelihood
Y <- Y + X %*% b # Test with nonzero mean to avoid errors otherwise hidden

resloglik_R <- limestest:::res_ll(XtX = crossprod(X),
                               XtY = crossprod(X, Y),
                               XtZ = crossprod(X, Z),
                               ZtZ = crossprod(Z),
                               crossprod(Y, Z),
                               Y = Y,
                               X = X,
                               Z = Z,
                               H = H,
                               Psi_r = Psi_cpp / psi_hat[r],
                               psi_r = psi_hat[r],
                               get_val = TRUE,
                               get_score = TRUE,
                               get_inf = TRUE)

resloglik_cpp <- limestest:::loglik_res(Y = Y,
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

cat("Difference in R and Cpp res. loglik: ", abs(resloglik_R$value - resloglik_cpp$value) , "\n")
cat("Max difference in R and Cpp res. score: ", max(abs(resloglik_R$score - resloglik_cpp$score)) , "\n")
cat("Max difference in R and Cpp res. information: ", max(abs(resloglik_R$inf_mat - resloglik_cpp$inf_mat)) , "\n")

