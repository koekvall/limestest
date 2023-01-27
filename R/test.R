library(foreign)
library(lme4)
fev1_dat <- read.dta("http://www.hsph.harvard.edu/fitzmaur/ala2e/fev1.dta")
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
            data = fev1_dat, REML = FALSE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
y <- getME(fit, "y")


random_effects <- ranef(fit, condVar = TRUE)
cov_matrices <- attr(random_effects[[1]], "postVar")

