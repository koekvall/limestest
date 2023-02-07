library(foreign)
library(lme4)
library(Matrix)
fev1_dat <- read.dta("http://www.hsph.harvard.edu/fitzmaur/ala2e/fev1.dta")
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
            data = fev1_dat, REML = FALSE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
y <- getME(fit, "y")


random_effects <- ranef(fit, condVar = TRUE)
cov_matrices <- attr(random_effects[[1]], "postVar")
H1 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 0, 0, 1), 2, 2))
Psi_hat <- Matrix::kronecker(Matrix::Diagonal(300), matrix(VarCorr(fit)$id, 2, 2))
psi0_hat <- attr(VarCorr(fit), "sc")^2
beta_hat <- fixef(fit)
e <- y - X %*% beta_hat
ZtZXe <- Matrix::crossprod(Z, cbind(Z, X, e))


#loglik(Z = Z, ZtZXe = ZtZXe, e = e, Psi0 = Psi_hat / psi0_hat, psi0 = psi0_hat);
score_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = cbind(H1, H2, H3),
          Psi0 = Psi_hat/psi0_hat, psi0 = psi0_hat, FALSE)
