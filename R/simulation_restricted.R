# set.seed(335)
# library(Matrix)
# num_reps <- 1e3
# num_cluster <- 100
# num_ind <- 10
# Z <- Matrix::bdiag(replicate(num_cluster, cbind(1, rnorm(num_ind)),
#                              simplify = FALSE))
# X <- matrix(rnorm(nrow(Z) * 2), nrow = nrow(Z), ncol = 2)
# Psi1 <- matrix(c(1, -0.8, -0.8, 1), 2, 2)
# psi0 <- 0.5
# Psi <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1)
# Psi0 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1 / psi0)
# H1 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(1, 0, 0, 0), 2, 2))
# H2 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 1, 1, 0), 2, 2))
# H3 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 0, 0, 1), 2, 2))
# H <- cbind(H1, H2, H3)
# R <-  Matrix::kronecker(Matrix::Diagonal(num_cluster), chol(Psi1))
#
# y_outside <- rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))

# ###############################################################################
# ## CHECK IF MC ESTIMATE OF FINF IS CLOSE #######################################
# ################################################################################
# one_sim_score <- function(seed){
#   set.seed(seed)
#   y <- rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))
#
#   limestest::res_ll(XtX = crossprod(X),
#                                  XtY = crossprod(X, y),
#                                  XtZ = crossprod(X, Z),
#                                  ZtZ = crossprod(Z),
#                                  YtZ = crossprod(y, Z),
#                                  Y = y,
#                                  X = X,
#                                  Z = Z,
#                                  H = H,
#                                  Psi0 = Psi0,
#                                  psi0 = psi0,
#                                  score = TRUE,
#                                  finf = TRUE,
#                                  lik = TRUE)$score
#
# }
# out <- sapply(1:num_reps, one_sim_score)
# MC_inf<- cov(t(out))
# AN_inf <- limestest::res_ll(XtX = crossprod(X),
#                            XtY = crossprod(X, y_outside),
#                            XtZ = crossprod(X, Z),
#                            ZtZ = crossprod(Z),
#                            YtZ = crossprod(y_outside, Z),
#                            Y = y_outside,
#                            X = X,
#                            Z = Z,
#                            H = H,
#                            Psi0 = Psi0,
#                            psi0 = psi0,
#                            score = TRUE,
#                            finf = TRUE,
#                            lik = TRUE)$finf
#
# MC_inf
# AN_inf
# ###############################################################################

################################################################################
## CHECK IF P-VALUES ARE UNIFORM ###############################################
################################################################################
# one_sim_test <- function(seed){
#   set.seed(seed)
#   y <- rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))
#
#   stuff <- limestest::res_ll(XtX = crossprod(X),
#                              XtY = crossprod(X, y),
#                              XtZ = crossprod(X, Z),
#                              ZtZ = crossprod(Z),
#                              YtZ = crossprod(y, Z),
#                              Y = y,
#                              X = X,
#                              Z = Z,
#                              H = H,
#                              Psi0 = Psi0,
#                              psi0 = psi0,
#                              score = TRUE,
#                              finf = TRUE,
#                              lik = TRUE)
#   test_stat <- as.vector(crossprod(stuff$score, solve(stuff$finf, stuff$score)))
#
#   pchisq(test_stat, df = 4, lower = F)
#
# }
# p_vals <- sapply(1:num_reps, one_sim_test)
# print(mean(p_vals < 0.05))
# plot(quantile(p_vals, seq(0.01, 0.99, 0.01)))

###############################################################################



# # TEST SIGMA CORRECT
# Sigma <- diag(psi0, nrow(Z)) + Z %*% tcrossprod(Psi, Z)
#
# MC <- t(matrix(replicate(1e4,
#                           as.vector(rnorm(nrow(Z), sd = sqrt(psi0)) +
#                                       Z %*% crossprod(R, rnorm(ncol(R))))),
#                ncol = 1e4))
# Sigma_mc <- cov(MC)
# (Sigma[1:10, 1:10] - Sigma_mc[1:10, 1:10]) / Sigma[1:10, 1:10]
# colMeans(MC)



### DERIVATIVE ON SIMULATED DATA ###############################################
# y <- rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))
#
# loglik_psi <- function(x){
#   psi0_arg <- x[1]
#   Psi1_arg <- matrix(c(x[2], x[3], x[3], x[4]), 2, 2)
#   Psi0_arg <-  Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1_arg / psi0_arg)
#   res_ll(XtX = crossprod(X),
#          XtY = crossprod(X, y),
#          XtZ = crossprod(X, Z),
#          ZtZ = crossprod(Z),
#          YtZ = crossprod(y, Z),
#          Y = y,
#          X = X,
#          Z = Z,
#          H = H,
#          Psi0 = Psi0_arg,
#          psi0 = psi0_arg,
#          score = TRUE,
#          finf = TRUE,
#          lik = TRUE)$ll
# }
#
# score_psi <- function(x){
#   psi0_arg <- x[1]
#   Psi1_arg <- matrix(c(x[2], x[3], x[3], x[4]), 2, 2)
#   Psi0_arg <-  Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1_arg / psi0_arg)
#   res_ll(XtX = crossprod(X),
#          XtY = crossprod(X, y),
#          XtZ = crossprod(X, Z),
#          ZtZ = crossprod(Z),
#          YtZ = crossprod(y, Z),
#          Y = y,
#          X = X,
#          Z = Z,
#          H = H,
#          Psi0 = Psi0_arg,
#          psi0 = psi0_arg,
#          score = TRUE,
#          finf = TRUE,
#          lik = TRUE)$score
# }
# test_point <- c(psi0, c(Psi1[1, 1], Psi1[1, 2], Psi1[2, 2]))
#
# score_psi(test_point)
# numDeriv::grad(loglik_psi, test_point)
