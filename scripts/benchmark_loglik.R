# Benchmark loglik with different get_* flag combinations
# to identify which parts of the C++ function are most expensive.

library(reconf)
library(lme4)

# ---------- Setup ----------
data(fev1)
fit <- lmer(
  logfev1 ~ age + ht + baseage + baseht + (age | id),
  data = fev1,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

Y     <- getME(fit, "y")
X     <- getME(fit, "X")
Z     <- getME(fit, "Z")
Hlist <- reconf:::get_Hlist_lmer(fit)
psi   <- reconf:::get_psi_hat_lmer(fit)
b     <- fixef(fit)
r     <- length(psi)
p     <- ncol(X)

H     <- methods::as(do.call(cbind, Hlist), "generalMatrix")
Psi_r <- (1 / psi[r]) * reconf:::Psi_from_H_cpp(psi[-r], H)
e     <- Y - X %*% b

XtX <- as.matrix(crossprod(X))
XtZ <- as.matrix(crossprod(X, Z))
ZtZ <- methods::as(crossprod(Z), "generalMatrix")

n_rep <- 200L

# ---------- Benchmark: ML loglik ----------
cat("=== ML loglik (", n_rep, "reps) ===\n\n")

# Value only
t_val <- system.time(for (i in seq_len(n_rep)) {
  reconf:::loglik(Psi_r = Psi_r, psi_r = psi[r], H = H, e = e,
                  X = X, Z = Z, XtX = XtX, XtZ = XtZ, ZtZ = ZtZ,
                  get_val = TRUE, get_score = FALSE, get_inf = FALSE,
                  expected = TRUE)
})
cat("val only:          ", round(t_val["elapsed"] / n_rep * 1000, 2), "ms\n")

# Value + score
t_vs <- system.time(for (i in seq_len(n_rep)) {
  reconf:::loglik(Psi_r = Psi_r, psi_r = psi[r], H = H, e = e,
                  X = X, Z = Z, XtX = XtX, XtZ = XtZ, ZtZ = ZtZ,
                  get_val = TRUE, get_score = TRUE, get_inf = FALSE,
                  expected = TRUE)
})
cat("val + score:       ", round(t_vs["elapsed"] / n_rep * 1000, 2), "ms\n")

# Value + score + expected info
t_vsi <- system.time(for (i in seq_len(n_rep)) {
  reconf:::loglik(Psi_r = Psi_r, psi_r = psi[r], H = H, e = e,
                  X = X, Z = Z, XtX = XtX, XtZ = XtZ, ZtZ = ZtZ,
                  get_val = TRUE, get_score = TRUE, get_inf = TRUE,
                  expected = TRUE)
})
cat("val + score + info:", round(t_vsi["elapsed"] / n_rep * 1000, 2), "ms\n")

# Value + score + observed info
t_vso <- system.time(for (i in seq_len(n_rep)) {
  reconf:::loglik(Psi_r = Psi_r, psi_r = psi[r], H = H, e = e,
                  X = X, Z = Z, XtX = XtX, XtZ = XtZ, ZtZ = ZtZ,
                  get_val = TRUE, get_score = TRUE, get_inf = TRUE,
                  expected = FALSE)
})
cat("val + score + obs: ", round(t_vso["elapsed"] / n_rep * 1000, 2), "ms\n")

cat("\n--- Incremental costs ---\n")
cat("Factorization + value:  ", round(t_val["elapsed"] / n_rep * 1000, 2), "ms\n")
cat("Adding score:           ", round((t_vs["elapsed"] - t_val["elapsed"]) / n_rep * 1000, 2), "ms\n")
cat("Adding expected info:   ", round((t_vsi["elapsed"] - t_vs["elapsed"]) / n_rep * 1000, 2), "ms\n")
cat("Observed vs expected:   ", round((t_vso["elapsed"] - t_vsi["elapsed"]) / n_rep * 1000, 2), "ms\n")

# ---------- Benchmark: REML loglik_res ----------
cat("\n\n=== REML loglik_res (", n_rep, "reps) ===\n\n")

XtY <- as.vector(crossprod(X, Y))
ZtY <- as.vector(crossprod(Z, Y))

t_val_r <- system.time(for (i in seq_len(n_rep)) {
  reconf:::loglik_res(Psi_r = Psi_r, psi_r = psi[r], H = H, Y = Y,
                      X = X, Z = Z, XtX = XtX, XtZ = XtZ, ZtZ = ZtZ,
                      XtY = XtY, ZtY = ZtY,
                      get_val = TRUE, get_score = FALSE, get_inf = FALSE)
})
cat("val only:          ", round(t_val_r["elapsed"] / n_rep * 1000, 2), "ms\n")

t_vs_r <- system.time(for (i in seq_len(n_rep)) {
  reconf:::loglik_res(Psi_r = Psi_r, psi_r = psi[r], H = H, Y = Y,
                      X = X, Z = Z, XtX = XtX, XtZ = XtZ, ZtZ = ZtZ,
                      XtY = XtY, ZtY = ZtY,
                      get_val = TRUE, get_score = TRUE, get_inf = FALSE)
})
cat("val + score:       ", round(t_vs_r["elapsed"] / n_rep * 1000, 2), "ms\n")

t_vsi_r <- system.time(for (i in seq_len(n_rep)) {
  reconf:::loglik_res(Psi_r = Psi_r, psi_r = psi[r], H = H, Y = Y,
                      X = X, Z = Z, XtX = XtX, XtZ = XtZ, ZtZ = ZtZ,
                      XtY = XtY, ZtY = ZtY,
                      get_val = TRUE, get_score = TRUE, get_inf = TRUE)
})
cat("val + score + info:", round(t_vsi_r["elapsed"] / n_rep * 1000, 2), "ms\n")

cat("\n--- Incremental costs ---\n")
cat("Factorization + value:  ", round(t_val_r["elapsed"] / n_rep * 1000, 2), "ms\n")
cat("Adding score:           ", round((t_vs_r["elapsed"] - t_val_r["elapsed"]) / n_rep * 1000, 2), "ms\n")
cat("Adding expected info:   ", round((t_vsi_r["elapsed"] - t_vs_r["elapsed"]) / n_rep * 1000, 2), "ms\n")

# ---------- Dimensions for context ----------
cat("\n\n=== Problem dimensions ===\n")
cat("n (observations):", length(Y), "\n")
cat("p (fixed effects):", p, "\n")
cat("q (random effects):", ncol(Z), "\n")
cat("r (covariance params):", r, "\n")
