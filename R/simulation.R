set.seed(335)
library(Matrix)
num_cluster <- 10
num_ind <- 5
Z <- Matrix::bdiag(replicate(num_cluster, cbind(1, rnorm(num_ind)),
                             simplify = FALSE))
X <- matrix(1, nrow = nrow(Z), ncol = 1)
Psi1 <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
psi0 <- 1
Psi <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1)
H1 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 0, 0, 1), 2, 2))
H <- cbind(H1, H2, H3)
R <-  Matrix::kronecker(Matrix::Diagonal(num_cluster), chol(Psi1))

one_sim <- function(seed){
  set.seed(seed)
  y <- rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))
  ZtZXe <- crossprod(Z, cbind(Z, X, y))

  score_inf <- limestest::score_psi(Z = Z, ZtZXe = ZtZXe, e = y, H = H,
                                    Psi0 = Psi / psi0, psi0 = psi0, finf = TRUE)
  #test_stat <- score_inf$score[3]^2 / score_inf$finf[3, 3]

  #pchisq(test_stat, 1, lower = F)

  c(score_inf$score, diag(score_inf$finf))
}

out <- do.call(rbind, lapply(1:1000, one_sim))
out[1, 5:8]
diag(cov(out[, 1:4]))

#p_vals <- sapply(1:1000, one_sim)
#print(mean(p_vals < 0.05))

# Test log-likelihood
# y <- rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(Psi)))
# ZtZXe <- crossprod(Z, cbind(Z, X, y))
# loglik_psi <- function(psi){
#   psi0 <- psi[1]
#   Psi1 <- matrix(c(psi[2], psi[3], psi[3], psi[4]), 2, 2)
#   Psi <-  Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1)
#   loglik(ZtZ = ZtZXe[, 1:ncol(Z)], Zte = ZtZXe[, ncol(ZtZXe)], e = y,
#          Psi0 = Psi / psi0, psi0 = psi0)
# }
#
# numerical_score <- numDeriv::grad(loglik_psi, c(psi0, Psi[1, 1],
#                                                 Psi[2, 1], Psi[2, 2]))
# analytical_score <- score_psi(Z = Z, ZtZXe = ZtZXe, e = y, H = H,
#                               Psi0 = Psi / psi0, psi0 = psi0, FALSE)$score
