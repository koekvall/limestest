library(tictoc)
library(lme4)
data("Pixel", package="nlme")
p <- 3
r <- 6
###############################################################################
# Unconstrained optimizations
###############################################################################

# lme4
tic()
fit <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side) + (1|Side:Dog), data = Pixel, REML = FALSE)
toc()
psi_hat_lme4 <- limestest:::get_psi_hat_lmer(fit)
b_hat_lme4 <- fixef(fit)

tic()
fit_REML <- update(fit, REML = TRUE)
toc()
psi_hat_lme4_REML <- limestest:::get_psi_hat_lmer(fit_REML)


precomp <- limestest:::get_precomp_lmer(fit)
Hlist <- limestest:::get_Hlist_lmer(fit)


tic()
fit_our <- limestest:::maximize_loglik(start_val = c(b_hat_lme4, psi_hat_lme4),
                            opt_idx = seq(p + r),
                            Y = getME(fit, "y"),
                            X = getME(fit, "X"),
                            Z = getME(fit, "Z"),
                            Hlist = Hlist,
                            expected = TRUE,
                            REML = FALSE)
toc()



tic()
fit_our_REML <- limestest:::maximize_loglik(start_val = c(rep(0, r - 1), 1),
                                       opt_idx = seq_along(psi_hat_lme4),
                                       Y = getME(fit, "y"),
                                       X = getME(fit, "X"),
                                       Z = getME(fit, "Z"),
                                       Hlist = Hlist,
                                       expected = TRUE,
                                       REML = TRUE)
toc()

cat("ML estimates from lme4: ", c(b_hat_lme4, psi_hat_lme4), "\n")
cat("ML estimates from our: ", fit_our$arg, "\n")
cat("REML estimates from lme4: ", psi_hat_lme4, "\n")
cat("REML estimates from our: ", fit_our_REML$arg, "\n")

###############################################################################
# Constrained optimizations -- No random effects
###############################################################################

tic()
# Testing no REs
fit_our_null <- limestest:::maximize_loglik(start_val = c(b_hat_lme4, rep(0, r - 1), psi_hat_lme4[r]),
                                       opt_idx = c(1:p, p + r),
                                       Y = getME(fit, "y"),
                                       X = getME(fit, "X"),
                                       Z = getME(fit, "Z"),
                                       Hlist = Hlist,
                                       expected = TRUE,
                                       REML = FALSE)
toc()
# Should have same coefs as the following
fit_lm <- lm(pixel ~ day + I(day^2), data = Pixel)

tic()
# Testing no REs with REML
fit_our_null_REML <- limestest:::maximize_loglik(start_val = c(rep(0, r - 1), psi_hat_lme4[r]),
                                            opt_idx = c(r),
                                            Y = getME(fit, "y"),
                                            X = getME(fit, "X"),
                                            Z = getME(fit, "Z"),
                                            Hlist = Hlist,
                                            expected = TRUE,
                                            REML = TRUE)
toc()

cat("No RE ML estimates from our: ", fit_our_null$arg, "\n")
cat("No RE lm estimates: ", coef(lm(pixel ~ day + I(day^2), data = Pixel)), "\n")
cat("REML estimate error variance from lm: ", sigma(lm(pixel ~ day + I(day^2), data = Pixel))^2, "\n")
cat("REML estimate error variance from our: ", fit_our_null_REML$arg[r], "\n")

###############################################################################
# Constrained optimizations -- Test a random effect variance (Side)
###############################################################################

# Null of no side random effect
fit_null_lme4 <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side), data = Pixel, REML = FALSE)
b_hat_null <- fixef(fit_null_lme4)
psi_hat_null <- limestest:::get_psi_hat_lmer(fit_null_lme4)
psi_hat_null <- append(psi_hat_null, 0, length(psi_hat_null) - 1)
fit_our_null <- limestest:::maximize_loglik(start_val = c(b_hat_null, psi_hat_null),
                                                 opt_idx = seq(p + r)[-(p + r - 1)],
                                                 Y = getME(fit, "y"),
                                                 X = getME(fit, "X"),
                                                 Z = getME(fit, "Z"),
                                                 Hlist = Hlist,
                                                 expected = TRUE,
                                                 REML = FALSE)



# Compare estimates
c(b_hat_null, psi_hat_null)
fit_our_null$arg


# Null of no side random effect REML
fit_null_lme4 <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side), data = Pixel, REML = TRUE)
psi_hat_null <- limestest:::get_psi_hat_lmer(fit_null_lme4)
psi_hat_null <- append(psi_hat_null, 0, length(psi_hat_null) - 1)
fit_our_null_REML <- limestest:::maximize_loglik(start_val = psi_hat_null,
                                            opt_idx = seq(r)[-(r - 1)],
                                            Y = getME(fit, "y"),
                                            X = getME(fit, "X"),
                                            Z = getME(fit, "Z"),
                                            Hlist = Hlist,
                                            expected = TRUE,
                                            REML = TRUE)


fit_our_null_REML_simple <- limestest:::maximize_loglik(start_val = psi_hat_null[-(r - 1)],
                                                 opt_idx = seq(r - 1),
                                                 Y = getME(fit_null_lme4, "y"),
                                                 X = getME(fit_null_lme4, "X"),
                                                 Z = getME(fit_null_lme4, "Z"),
                                                 Hlist = limestest:::get_Hlist_lmer(fit_null_lme4),
                                                 expected = TRUE,
                                                 REML = TRUE)

# Compare estimates
psi_hat_null
fit_our_null_REML_simple$arg

