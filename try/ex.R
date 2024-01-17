library(lme4)
library(Matrix)
#library(limestest)
library(rstiefel)
data(fev1)
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
            data = fev1, REML = FALSE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
y <- getME(fit, "y")

random_effects <- ranef(fit, condVar = TRUE)
cov_matrices <- attr(random_effects[[1]], "postVar")
H1 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 0, 0, 1), 2, 2))
H <- cbind(H1, H2, H3)
Psi1_hat <-  matrix(VarCorr(fit)$id, 2, 2)
Psi_hat <- Matrix::kronecker(Matrix::Diagonal(300), Psi1_hat)
psi0_hat <- attr(VarCorr(fit), "sc")^2
beta_hat <- fixef(fit)
e <- y - X %*% beta_hat
ZtZXe <- Matrix::crossprod(Z, cbind(Z, X, e))

# Loglik for psi^0 at beta = beta_hat
loglik_test <- function(x){
  Psi1_x <- matrix(c(x[2], x[3], x[3], x[4]), 2, 2)
  Psi_x <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_x)
  loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
             Psi0 = Psi_x / x[1], psi0 = x[1], loglik = TRUE,
             score = TRUE, finf = TRUE, expected = FALSE)
}


psi0_test <- psi0_hat
Psi1_test <- Psi_hat[1:2, 1:2]
Psi_test <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_test)

library(Rcpp)
sourceCpp("./src/functions_cpp.cpp")
source("./R/functions.R")

starttime <- Sys.time()
for (i in 1:100) {
  l <- loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                  Psi0 = Psi_test / psi0_test, psi0 = psi0_test, loglik = T,
                  finf = T, expected=T)
}
Sys.time()- starttime # 1.69 s

starttime <- Sys.time()
for (i in 1:100) {
lr <- loglik_psiRcpp(Z = Z, #ZtZXe = ZtZXe,
                     e = e, H = H,
                     Psi0 = Psi_test / psi0_test, psi0 = psi0_test, loglik = T,
                     finf = T, expected=T)
}
Sys.time()- starttime # 0.56 s
#res_ll
# check I(1:r,1:r) in finf=T
starttime <- Sys.time()
for (i in 1:100) {
  l <- res_ll(XtX = crossprod(X), XtY = crossprod(X,y),XtZ = crossprod(X,Z),ZtZ = crossprod(Z),YtZ = crossprod(y,Z),
              Y=y, X=X, Z = Z, H = H, Psi0 = Psi_test / psi0_test, psi0 = psi0_test, lik = T,
              score = T, finf = F)
}
Sys.time()- starttime # 11.02915 secs
starttime <- Sys.time()
for (i in 1:100) {
  lr <- res_llRcpp(X=X, Y=y, Z = Z, H = H,
                   Psi0 = Psi_test / psi0_test, psi0 = psi0_test, lik = T,
                   score = T, finf = F)
}
Sys.time()- starttime # 1.423625 secs



