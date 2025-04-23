library(lme4)
library(Matrix)
library(limestest)
library(rstiefel)
library(merDeriv)
data(fev1)
###############################################################################
# Test derivatives of regular log-likelihood
###############################################################################
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
            data = fev1, REML = FALSE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
Y <- getME(fit, "y")

H1 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 0, 0, 1), 2, 2))
H <- cbind(H1, H2, H3)

Psi1_hat <-  matrix(VarCorr(fit)$id, 2, 2)
Psi_hat <- Matrix::kronecker(Matrix::Diagonal(300), Psi1_hat)
psir_hat <- attr(VarCorr(fit), "sc")^2
beta_hat <- fixef(fit)
e <- Y - X %*% beta_hat
ZtZXe <- Matrix::crossprod(Z, cbind(Z, X, e))
# Loglik for psi^0 at beta = beta_hat
loglik_test <- function(x){
  Psi1_x <- matrix(c(x[1], x[2], x[2], x[3]), 2, 2)
  Psi_x <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_x)
  limestest:::loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
         Psi_r = Psi_x / x[4], psi_r = x[4], get_val = TRUE,
         get_score = TRUE, get_inf = TRUE, expected = FALSE)$value
}

# Test derivatives at MLE
psir_test <- psir_hat
Psi1_test <- Psi_hat[1:2, 1:2]
Psi_test <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_test)
test_point <- c(Psi_test[1, 1], Psi_test[2, 1], Psi_test[2, 2], psir_test)
numerical_score <- numDeriv::grad(loglik_test, test_point)
analytical_score <- limestest:::loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                              Psi_r = Psi_test / psir_test, psi_r = psir_test, get_val = F,
                              get_score = T, get_inf = F)$score
mer_score <- colSums(merDeriv::estfun.lmerMod(fit))[6:9]

numerical_hess <- numDeriv::hessian(loglik_test, test_point)
analytical_hess <- -limestest:::loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                               Psi_r = Psi_test / psir_test, psi_r = psir_test, get_val = F,
                               get_score = T, get_inf = TRUE, expected = F)$inf_mat
mer_hess <- -solve(merDeriv::vcov.lmerMod(object = fit, full = TRUE, expected = FALSE))[6:9, 6:9]

cat("The max absolute difference between numerical and analytical score at MLE is: ",
    max(abs(numerical_score - analytical_score)), "\n")
cat("The max relative difference between numerical and analytical score at MLE is: ",
    max(abs(numerical_score - analytical_score) / numerical_score), "\n")
cat("The max absolute difference between merDeriv and analytical score at MLE is: ",
    max(abs(mer_score - analytical_score)), "\n")
cat("The max relative difference between merDeriv and analytical score at MLE is: ",
    max(abs(mer_score - analytical_score) / mer_score), "\n")


cat("The max absolute difference between numerical and analytical Hessian at MLE is: ",
    max(abs(numerical_hess - analytical_hess)), "\n")
cat("The max relative difference between numerical and analytical Hessian at MLE is: ",
    max(abs(numerical_hess - analytical_hess) / numerical_hess), "\n")
cat("The max absolute difference between merDeriv and analytical Hessian at MLE is: ",
    max(abs(mer_hess - analytical_hess)), "\n")
cat("The max relative difference between merDeriv and analytical Hessian at MLE is: ",
    max(abs(mer_hess - analytical_hess) / mer_hess), "\n")

# Test wrapper function
value_raw <- limestest:::loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                               Psi_r = Psi_test / psir_test, psi_r = psir_test, get_val = TRUE,
                               get_score = T, get_inf = F)$value
Hlist <- limestest:::get_Hlist_lmer(fit)
psi_hat <- limestest:::get_psi_hat_lmer(fit)
ll_things <- limestest::loglikelihood(psi = psi_hat,
                                      b = lme4::getME(fit, "beta"),
                                      Y = Y,
                                      X = X,
                                      Z = Z,
                                      Hlist = Hlist,
                                      REML = FALSE,
                                      expected = FALSE,
                                      get_inf = TRUE)

cat("The difference in raw and wrapper value is: ", abs(ll_things$value - value_raw), "\n")
cat("The max difference in raw and wrapper score is: ", max(abs(ll_things$score - analytical_score)), "\n")
cat("The max difference in raw and wrapper Hessian is: ", max(abs(-ll_things$inf_mat - analytical_hess)), "\n")

analytical_inf <- limestest:::loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                                           Psi_r = Psi_test / psir_test, psi_r = psir_test, get_val = F,
                                           get_score = T, get_inf = TRUE, expected = T)$inf_mat
ll_things <- limestest::loglikelihood(psi = psi_hat,
                                      b = lme4::getME(fit, "beta"),
                                      Y = Y,
                                      X = X,
                                      Z = Z,
                                      Hlist = Hlist,
                                      REML = FALSE,
                                      expected = TRUE,
                                      get_inf = TRUE)

cat("The max difference in raw and wrapper Information is: ", max(abs(ll_things$inf_mat - analytical_inf)), "\n")

# Test derivatives at random point
psir_test <- runif(1, min = psir_hat / 10, psir_hat * 10)
ed <- eigen(Psi1_hat)
U <- rustiefel(2, 2)
Psi1_test <- crossprod(U, diag(runif(2, min = min(ed$values) / 10, max = max(ed$values) * 10))) %*% U
Psi_test <- Matrix::kronecker(Matrix::Diagonal(300), Psi1_test)
test_point <- c(Psi_test[1, 1], Psi_test[2, 1], Psi_test[2, 2], psir_test)
numerical_score <- numDeriv::grad(loglik_test, test_point)
analytical_score <- limestest:::loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                               Psi_r = Psi_test / psir_test, psi_r = psir_test, get_val = F,
                               get_score = T, get_inf = F)$score

numerical_hess <- numDeriv::hessian(loglik_test, test_point)
analytical_hess <- -limestest:::loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                               Psi_r = Psi_test / psir_test, psi_r = psir_test, get_val = F,
                               get_score = T, get_inf = T, expected = F)$inf_mat

cat("The max absolute difference between numerical and analytical score at random point is: ",
    max(abs(numerical_score - analytical_score)), "\n")
cat("The max relative difference between numerical and analytical score at random point is: ",
    max(abs(numerical_score - analytical_score) / numerical_score), "\n")



cat("The max absolute difference between numerical and analytical Hessian at random point is: ",
    max(abs(numerical_hess - analytical_hess)), "\n")
cat("The max relative difference between numerical and analytical Hessian at random point is: ",
    max(abs(numerical_hess - analytical_hess) / numerical_hess), "\n")

###############################################################################
# Test derivatives of restricted log-likelihood
###############################################################################
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
            data = fev1, REML = TRUE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
Y <- getME(fit, "y")

H1 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 0, 0, 1), 2, 2))
H <- cbind(H1, H2, H3)

Psi1_hat <-  matrix(VarCorr(fit)$id, 2, 2)
Psi_hat <- Matrix::kronecker(Matrix::Diagonal(300), Psi1_hat)
psir_hat <- attr(VarCorr(fit), "sc")^2

loglik_test <- function(x){
  Psi1_x <- matrix(c(x[1], x[2], x[2], x[3]), 2, 2)
  Psi_x <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_x)
  limestest:::res_ll(XtX = crossprod(X), XtY = crossprod(X, Y), XtZ = crossprod(X, Z),
         ZtZ = crossprod(Z), YtZ = crossprod(Y, Z), Y = Y, X = X, Z = Z, H = H,
             Psi_r = Psi_x / x[4], psi_r = x[4], get_val = TRUE,
             get_score = FALSE, get_inf = FALSE)$value
}

# Test derivatives at MLE
psir_test <- psir_hat
Psi1_test <- Psi_hat[1:2, 1:2]
Psi_test <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_test)
test_point <- c(Psi_test[1, 1], Psi_test[2, 1], Psi_test[2, 2], psir_test)
numerical_score <- numDeriv::grad(loglik_test, test_point)

analytical_score <- limestest:::res_ll(XtX = crossprod(X), XtY = crossprod(X, Y), XtZ = crossprod(X, Z),
                           ZtZ = crossprod(Z), YtZ = crossprod(Y, Z), Y = Y,
                           X = X, Z = Z, H = H, Psi_r = Psi_test / psir_test,
                           psi_r = psir_test, get_val = FALSE,
                           get_score = TRUE, get_inf = FALSE)$score
mer_score <- colSums(merDeriv::estfun.lmerMod(fit))[6:9]

numerical_hess <- numDeriv::hessian(loglik_test, test_point)


mer_inf <- solve(merDeriv::vcov.lmerMod(object = fit, full = TRUE, expected = TRUE)[6:9, 6:9])
analytical_inf <- limestest:::res_ll(XtX = crossprod(X), XtY = crossprod(X, Y), XtZ = crossprod(X, Z),
                                             ZtZ = crossprod(Z), YtZ = crossprod(Y, Z), Y = Y,
                                             X = X, Z = Z, H = H, Psi_r = Psi_test / psir_test,
                                             psi_r = psir_test, get_val = FALSE,
                                             get_score = TRUE, get_inf = TRUE)$inf_mat

cat("The max absolute difference between numerical and analytical score at MLE is: ",
    max(abs(numerical_score - analytical_score)), "\n")
cat("The max relative difference between numerical and analytical score at MLE is: ",
    max(abs(numerical_score - analytical_score) / numerical_score), "\n")
cat("The max absolute difference between merDeriv and analytical score at MLE is: ",
    max(abs(mer_score - analytical_score)), "\n")
cat("The max relative difference between merDeriv and analytical score at MLE is: ",
    max(abs(mer_score - analytical_score) / mer_score), "\n")

cat("The max absolute difference between merDeriv and analytical information at MLE is: ",
    max(abs(mer_inf - analytical_inf)), "\n")
cat("The max relative difference between merDeriv and analytical information at MLE is: ",
    max(abs(mer_inf - analytical_inf) / mer_inf), "\n")

# Test wrapper function
psi_hat <- limestest:::get_psi_hat_lmer(fit)
value_raw <- limestest:::res_ll(XtX = crossprod(X),
                                XtY = crossprod(X, Y),
                                XtZ = crossprod(X, Z),
                                ZtZ = crossprod(Z),
                                YtZ = crossprod(Y, Z),
                                Y = Y,
                                X = X,
                                Z = Z,
                                H = H,
                                Psi_r = Psi_test / psir_test,
                                psi_r = psir_test,
                                get_val = TRUE,
                                get_score = FALSE,
                                get_inf = FALSE)$value
Hlist <- limestest:::get_Hlist_lmer(fit)
ll_things <- limestest::loglikelihood(psi = psi_hat,
                                      Y = Y,
                                      X = X,
                                      Z = Z,
                                      Hlist = Hlist,
                                      REML = TRUE)

cat("The difference in raw and wrapper value is: ", abs(ll_things$value - value_raw), "\n")

cat("The max difference in raw and wrapper score is: ", max(abs(ll_things$score - analytical_score)), "\n")
cat("The max difference in raw and wrapper Hessian is: ", max(abs(-ll_things$inf_mat - analytical_hess)), "\n")

analytical_inf <- limestest:::loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                                         Psi_r = Psi_test / psir_test, psi_r = psir_test, get_val = F,
                                         get_score = T, get_inf = TRUE, expected = T)$inf_mat
ll_things <- limestest::loglikelihood(psi = psi_hat,
                                      b = lme4::getME(fit, "beta"),
                                      Y = Y,
                                      X = X,
                                      Z = Z,
                                      Hlist = Hlist,
                                      REML = FALSE,
                                      expected = TRUE,
                                      get_inf = TRUE)

cat("The max difference in raw and wrapper Information is: ", max(abs(ll_things$inf_mat - analytical_inf)), "\n")


# Test derivatives at random point
psir_test <- runif(1, min = psir_hat / 10, psir_hat * 10)
ed <- eigen(Psi1_hat)
U <- rustiefel(2, 2)
Psi1_test <- crossprod(U, diag(runif(2, min = min(ed$values) / 10, max = max(ed$values) * 10))) %*% U
Psi_test <- Matrix::kronecker(Matrix::Diagonal(300), Psi1_test)
test_point <- c(Psi_test[1, 1], Psi_test[2, 1], Psi_test[2, 2], psir_test)
numerical_score <- numDeriv::grad(loglik_test, test_point)

analytical_score <- limestest:::res_ll(XtX = crossprod(X), XtY = crossprod(X, Y), XtZ = crossprod(X, Z),
                           ZtZ = crossprod(Z), YtZ = crossprod(Y, Z), Y = Y,
                           X = X, Z = Z, H = H, Psi_r = Psi_test / psir_test,
                           psi_r = psir_test, get_val = FALSE,
                           get_score = TRUE, get_inf = FALSE)$score

numerical_hess <- numDeriv::hessian(loglik_test, test_point)

analytical_inf <- limestest:::res_ll(XtX = crossprod(X), XtY = crossprod(X, Y), XtZ = crossprod(X, Z),
                         ZtZ = crossprod(Z), YtZ = crossprod(Y, Z), Y = Y,
                         X = X, Z = Z, H = H, Psi_r = Psi_test / psir_test,
                         psi_r = psir_test, get_val = FALSE,
                         get_score = TRUE, get_inf = TRUE)$inf_mat

cat("The max absolute difference between numerical and analytical score at random point is: ",
    max(abs(numerical_score - analytical_score)), "\n")
cat("The max relative difference between numerical and analytical score at random point is: ",
    max(abs(numerical_score - analytical_score) / numerical_score), "\n")


# Test wrapper function



###############################################################################
# Simulation for Hessian
###############################################################################

# Simulation for Hessian
# n_reps <- 1e4
# n1 <- 20
# n2 <- 3
# n <- n1 * n2
# p <- 1
# X <- matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p)
# Z <- Matrix::bdiag(replicate(n1, cbind(1, rnorm(n2)),
#                                   simplify = FALSE))
# ZtZ <- crossprod(Z)
# ZtX <- crossprod(Z, X)
# psi0_star <- 0.1
# Psi1_star <- matrix(c(0.5, 0.2, 0.2, 0.5), 2, 2)
# Psi_star <- Matrix::kronecker(Matrix::Diagonal(n1), Psi1_star)
# Psi0_star <- Psi_star / psi0_star
# R_star <- Matrix::kronecker(Matrix::Diagonal(n1), chol(Psi1_star))
# H1 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(1, 0, 0, 0), ncol = 2,nrow =  2))
# H2 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(0, 1, 1, 0), ncol = 2,nrow =  2))
# H3 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(0, 0, 0, 1), ncol = 2,nrow =  2))
# H <- cbind(H1, H2, H3)
#
# hess_list <- list()
# y_list <- list()
# for(ii in 1:n_reps){
#   y <- rnorm(nrow(Z), sd = sqrt(psi0_star)) +
#     Z %*% crossprod(R_star,  rnorm(ncol(R_star))) # beta_star = 0
#   Zte <- crossprod(Z, y)
#   hess_list[[ii]]<- limestest::loglik_psi(Z = Z,
#                         ZtZXe = cbind(ZtZ, ZtX, Zte),
#                         e = y,
#                         H = H,
#                         Psi0 = Psi0_star,
#                         psi0 = psi0_star,
#                         loglik = FALSE,
#                         score = FALSE,
#                         finf = TRUE,
#                         expected = FALSE)$finf
#   y_list[[ii]] <- y
# }
#
# MC_inf <- Reduce("+", hess_list) / length(hess_list)
# analytical_inf <- limestest::loglik_psi(Z = Z,
#                                         ZtZXe = cbind(ZtZ, ZtX, Zte),
#                                         e = y,
#                                         H = H,
#                                         Psi0 = Psi0_star,
#                                         psi0 = psi0_star,
#                                         loglik = FALSE,
#                                         score = FALSE,
#                                         finf = TRUE,
#                                         expected = TRUE)$finf
#
#








