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
fit_lme4 <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side) + 
    (1|Side:Dog), data = Pixel, REML = FALSE)
toc()
psi_hat <- limestest:::get_psi_hat_lmer(fit_lme4)
b_hat <- fixef(fit_lme4)

tic()
fit_lme4_REML <- update(fit_lme4, REML = TRUE)
toc()
psi_hat_REML <- limestest:::get_psi_hat_lmer(fit_lme4_REML)


precomp <- limestest:::get_precomp_lmer(fit_lme4)
Hlist <- limestest:::get_Hlist_lmer(fit_lme4)


tic()
fit_our <- limestest:::maximize_loglik(start_val = c(b_hat, psi_hat),
                            opt_idx = seq(p + r),
                            Y = getME(fit_lme4, "y"),
                            X = getME(fit_lme4, "X"),
                            Z = getME(fit_lme4, "Z"),
                            Hlist = Hlist,
                            expected = TRUE,
                            REML = FALSE)
toc()



tic()
fit_our_REML <- limestest:::maximize_loglik(start_val = c(rep(0, r - 1), 1),
                                       opt_idx = seq_along(psi_hat),
                                       Y = getME(fit_lme4, "y"),
                                       X = getME(fit_lme4, "X"),
                                       Z = getME(fit_lme4, "Z"),
                                       Hlist = Hlist,
                                       expected = TRUE,
                                       REML = TRUE)
toc()

cat("ML estimates from lme4: ", c(b_hat, psi_hat), "\n")
cat("ML estimates from our: ", fit_our$arg, "\n")
cat("REML estimates from lme4: ", psi_hat_REML, "\n")
cat("REML estimates from our: ", fit_our_REML$arg, "\n")

###############################################################################
# Constrained optimizations -- No random effects
###############################################################################

tic()
# Testing no REs
fit_our_no_re <- limestest:::maximize_loglik(start_val = c(b_hat, rep(0, r - 1), psi_hat[r]),
                                       opt_idx = c(1:p, p + r),
                                       Y = getME(fit_lme4, "y"),
                                       X = getME(fit_lme4, "X"),
                                       Z = getME(fit_lme4, "Z"),
                                       Hlist = Hlist,
                                       expected = TRUE,
                                       REML = FALSE)
toc()
# Should have same coefs as the following
fit_lm <- lm(pixel ~ day + I(day^2), data = Pixel)

tic()
# Testing no REs with REML
fit_our_no_re_REML <- limestest:::maximize_loglik(start_val = c(rep(0, r - 1), psi_hat[r]),
                                            opt_idx = c(r),
                                            Y = getME(fit_lme4, "y"),
                                            X = getME(fit_lme4, "X"),
                                            Z = getME(fit_lme4, "Z"),
                                            Hlist = Hlist,
                                            expected = TRUE,
                                            REML = TRUE)
toc()

cat("No RE ML estimates from our: ", fit_our_no_re$arg, "\n")
cat("No RE lm estimates: ", coef(fit_lm), "\n")
cat("REML estimate error variance from lm: ", sigma(fit_lm)^2, "\n")
cat("REML estimate error variance from our: ", fit_our_no_re_REML$arg[r], "\n")

###############################################################################
# Constrained optimizations -- Test a random effect variance (Side)
###############################################################################

# Null of no side random effect
# In lme4, the first random effect parameter is the variance for (1|Side:Dog)
fit_lme4_null <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side),
                      data = Pixel, REML = FALSE)
b_hat_null <- fixef(fit_lme4_null)
psi_hat_null <- limestest:::get_psi_hat_lmer(fit_lme4_null)
psi_hat_null <- append(psi_hat_null, 0, 0)

fit_our_null <- limestest:::maximize_loglik(start_val = c(b_hat_null, psi_hat_null),
                                                 opt_idx = seq(p + r)[-(p + 1)],
                                                 Y = getME(fit_lme4, "y"),
                                                 X = getME(fit_lme4, "X"),
                                                 Z = getME(fit_lme4, "Z"),
                                                 Hlist = Hlist,
                                                 expected = TRUE,
                                                 REML = FALSE)



# Compare estimates
cat("ML estimates under null lme4: ", c(b_hat_null, psi_hat_null), "\n")
cat("ML estimates under null our: ", fit_our_null$arg, "\n")


# Null of no side random effect REML
fit_lme4_null_REML <- update(fit_lme4_null, REML = TRUE)
psi_hat_null_REML <- limestest:::get_psi_hat_lmer(fit_lme4_null_REML)
psi_hat_null_REML <- append(psi_hat_null_REML, 0, 0)
fit_our_null_REML <- limestest:::maximize_loglik(start_val = psi_hat_null_REML,
                                            opt_idx = seq(r)[-1],
                                            Y = getME(fit_lme4, "y"),
                                            X = getME(fit_lme4, "X"),
                                            Z = getME(fit_lme4, "Z"),
                                            Hlist = Hlist,
                                            expected = TRUE,
                                            REML = TRUE)

# Compare estimates
cat("REML estimates under null lme4: ", psi_hat_null_REML, "\n")
cat("REML estimates under null our: ", fit_our_null_REML$arg, "\n")

