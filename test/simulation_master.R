library(Matrix)
library(doParallel)
library(lme4)
library(limestest)
library(doRNG)

set.seed(1336)

###############################################################################
# Settings
###############################################################################

# Common settings
num_cores <- 8
num_reps <- 1e2
out_dir <- "~/GitHub/lmm-crit-suppl/" # has to end in /
fun_dir <- "~/GitHub/limestest/test/" # has to end in /

cl <- makeCluster(num_cores)
registerDoParallel(cl)

source(paste0(fun_dir, "sim_funs.R"))

# Common data generating parameters
n1 <- 20
n2 <- 5
n <- n1 * n2
p <- 2
psi0 <- 1
X <- cbind(1, matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p))
XtX <- crossprod(X)
b <- rnorm(ncol(X))
Xb <- X %*% b

###############################################################################
# Independent cluster simulation with correlation
###############################################################################
Z <- Matrix::bdiag(lapply(1:n1, function(ii)X[((ii - 1) * n2 + 1):(ii * n2), 1:2]))
ZtZ <- crossprod(Z)
XtZ <- crossprod(X, Z)
H1 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(0, 0, 0, 1), 2, 2))
H <- as(cbind(H1, H2, H3), "sparseMatrix")

correlations <- c(0.99, 0.95, seq(0.9, -0.9, by  =  -0.1), -0.95, -0.99)
num_settings <- length(correlations)
out_mat_corr <- matrix(0, nrow = num_reps * length(correlations), ncol = 10)
for(jj in 1:num_settings){
  Psi1 <- matrix(c(1, correlations[jj], correlations[jj], 1), 2, 2)
  Psi1_eig <- eigen(Psi1)
  R1 <- tcrossprod(diag(sqrt(Psi1_eig$values), ncol(Psi1)), Psi1_eig$vectors)
  R <-  as(Matrix::kronecker(Matrix::Diagonal(n1), R1), "sparseMatrix")
  Psi <- Matrix::kronecker(Matrix::Diagonal(n1), Psi1)
  Sigma <- diag(psi0, n) + Z %*% tcrossprod(Psi, Z)
  M <- solve(Sigma, X)
  make_beta_mat <- solve(crossprod(X, M), t(M))

  out_idx <- ((jj - 1) * num_reps + 1):(jj * num_reps)
  out_mat_corr[out_idx, ] <- foreach(kk = 1:num_reps, .combine = rbind,
                        .errorhandling = "remove",
                        .packages = c("limestest", "lme4", "Matrix")) %dorng%{
                          c(do_one_clust_sim(inner_seed = kk, n1, n2, X, Z, b, Psi1, psi0, H, R, XtX, XtZ,
                                           ZtZ, Xb, Sigma, make_beta_mat),
                            n1, n2, p, psi0, correlations[jj])
                        }

}
###############################################################################


###############################################################################
# Independent cluster simulation without correlation
###############################################################################
Z <- Matrix::bdiag(lapply(1:n1, function(ii)X[((ii - 1) * n2 + 1):(ii * n2), 1:2]))
ZtZ <- crossprod(Z)
XtZ <- crossprod(X, Z)
H1 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(0, 0, 0, 1), 2, 2))
H <- as(cbind(H1, H2, H3), "sparseMatrix")

variances <- c(0.01, 0.05, 0.1, 0.15, 0.2, seq(0.3, 1, by = 0.1))
num_settings <- length(variances)
out_mat_indep <- matrix(0, nrow = num_reps * length(variances), ncol = 10)
for(jj in 1:num_settings){
  Psi1 <- matrix(c(variances[jj], 0, 0, variances[jj]), 2, 2)
  R1 <- sqrt(Psi1)
  R <-  as(Matrix::kronecker(Matrix::Diagonal(n1), R1), "sparseMatrix")
  Psi <- Matrix::kronecker(Matrix::Diagonal(n1), Psi1)
  Sigma <- diag(psi0, n) + Z %*% tcrossprod(Psi, Z)
  M <- solve(Sigma, X)
  make_beta_mat <- solve(crossprod(X, M), t(M))

  out_idx <- ((jj - 1) * num_reps + 1):(jj * num_reps)
  out_mat_indep[out_idx, ] <- foreach(kk = 1:num_reps, .combine = rbind,
                                     .errorhandling = "remove",
                                     .packages = c("limestest", "lme4", "Matrix")) %dorng%{
                                       c(do_one_clust_sim(inner_seed = kk, n1, n2, X, Z, b, Psi1, psi0, H, R, XtX, XtZ,
                                                          ZtZ, Xb, Sigma, make_beta_mat),
                                         n1, n2, p, psi0, variances[jj])
                                     }
}
###############################################################################


###############################################################################
# Crossed random effects simulation
###############################################################################
Z <- as(cbind(Matrix::kronecker(Matrix::diag(1, nrow = n1), matrix(1, n2, 1)),
              Matrix::kronecker(matrix(1, n1, 1), Matrix::diag(1, nrow = n2))),
        "sparseMatrix")
ZtZ <- crossprod(Z)
XtZ <- crossprod(X, Z)
H1 <- Matrix::Diagonal(n1 + n2, c(rep(1, n1), rep(0, n2)))
H2 <- Matrix::Diagonal(n1 + n2, c(rep(0, n1), rep(1, n2)))
H <- as(cbind(H1, H2), "sparseMatrix")

variances <- c(0.01, 0.05, 0.1, 0.15, 0.2, seq(0.3, 1, by = 0.1))
num_settings <- length(variances)
out_mat_cross <- matrix(0, nrow = num_reps * length(variances), ncol = 10)
for(jj in 1:num_settings){
  psi1 <- c(variances[jj], variances[jj])
  Psi <- Matrix::Diagonal(n1 + n2, c(rep(psi1[1], n1), rep(psi1[2], n2)))
  Sigma <- diag(psi0, n) + Z %*% tcrossprod(Psi, Z)
  M <- solve(Sigma, X)
  make_beta_mat <- solve(crossprod(X, M), t(M))

  out_idx <- ((jj - 1) * num_reps + 1):(jj * num_reps)
  out_mat_cross[out_idx, ] <- foreach(kk = 1:num_reps, .combine = rbind,
                                      .errorhandling = "remove",
                                      .packages = c("limestest", "lme4", "Matrix")) %dorng%{
                                        c(do_one_cross_sim(kk, n1, n2, X, Z, b, psi1, psi0, H, XtX, XtZ,
                                                           ZtZ, Xb, Sigma, make_beta_mat),
                                          n1, n2, p, psi0, variances[jj])
                                      }
}
###############################################################################

stopCluster(cl)

###############################################################################
# Output
###############################################################################
out <- tidyr::tibble(tidyr::as_tibble(rbind(out_mat_corr, out_mat_indep, out_mat_cross)),
                     "type" = c(rep("corr", nrow(out_mat_corr)),
                                rep("inde", nrow(out_mat_indep)),
                                rep("cross", nrow(out_mat_cross))))
rm(out_mat_corr, out_mat_indep, out_mat_cross)
colnames(out) <- c("LL", "RLL", "LRT", "WLD", "seed", "n1", "n2", "p", "psi0", "param", "type")

today <- as.numeric(format(Sys.time(), "%H%d%m%y"))
saveRDS(out, paste0(out_dir, "lmm_sims_", today, ".Rds"))
