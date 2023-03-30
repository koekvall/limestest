library(trust)
library(lme4)
library(Matrix)
data(fev1)

# INDEP REs
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (1|id) +
              (0 + age|id), data = fev1, REML = TRUE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
Y <- getME(fit, "y")

H1 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(1, 0, 0, 0), 2, 2))
# H2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 0, 0, 1), 2, 2))
H <- as(cbind(H1, H3), "sparseMatrix")
VC <- as.data.frame(VarCorr(fit))
Psi_hat <- Matrix::kronecker(Matrix::Diagonal(300), diag(VC[1:2, 4], 2))
psi0_hat <- attr(VarCorr(fit), "sc")^2
psi_hat <- c(psi0_hat, Psi_hat[1, 1], Psi_hat[2, 2])
beta_hat <- fixef(fit)
e <- Y - X %*% beta_hat
XtX <- Matrix::crossprod(X)
XtY <- Matrix::crossprod(X, Y)
XtZ <- Matrix::crossprod(X, Z)
ZtZ <- Matrix::crossprod(Z)
YtZ <- Matrix::crossprod(Y, Z)


loglik_test <- function(x){
  if(any(x < 0)){
    return(-Inf)
  }
  psi0_arg <- x[1]
  Psi1_arg <- matrix(c(x[2], 0, 0, x[3]), 2, 2)
  Psi0_arg <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_arg / psi0_arg)
  limestest::res_ll(XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H, Psi0_arg, psi0_arg, lik = TRUE, score = FALSE,
         finf = FALSE)$ll
}

score_test <- function(x){
  psi0_arg <- x[1]
  Psi1_arg <- matrix(c(x[2], 0, 0, x[3]), 2, 2)
  Psi0_arg <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_arg / psi0_arg)
  limestest::res_ll(XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H, Psi0_arg, psi0_arg, lik = FALSE, score = TRUE,
         finf = FALSE)$score
}

my_obj <- function(x){

    psi0_arg <- x[1]
    Psi1_arg <- matrix(c(x[2], 0, 0, x[3]), 2, 2)
    Psi0_arg <-  as(Matrix::kronecker(Matrix::Diagonal(300), Psi1_arg / psi0_arg), "sparseMatrix")
    stuff <- limestest::res_ll(XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H, Psi0_arg, psi0_arg, lik = TRUE, score = TRUE,
           finf = TRUE)
    return(list("value" = -stuff$ll, "gradient" = -stuff$score, "hessian" = as.matrix(stuff$finf)))

}

trust_fit <- trust(my_obj, psi_hat, 1, 100)
optim_fit <- optim(par = psi_hat, fn = function(x)-loglik_test(x), gr = function(x)-score_test(x),
                   lower = c(0, 0, 0), method = "L-BFGS-B")

my_obj(trust_fit$argument)$value
-my_obj(optim_fit$par)$value
my_obj(psi_hat)$value
