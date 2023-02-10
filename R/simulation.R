set.seed(335)
library(Matrix)
num_cluster <- 100
num_ind <- 5
Z <- Matrix::bdiag(replicate(num_cluster, matrix(rnorm(num_ind * 2), num_ind, 2),
                             simplify = FALSE))
X <- matrix(1, nrow = nrow(Z), ncol = 1)
Psi1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
psi0 <- 1.5
Psi <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1)
H1 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 0, 0, 1), 2, 2))
H <- cbind(H1 + H3, H2) # Same variance
R <- chol(Psi)
one_sim <- function(seed){
  set.seed(seed)
  y <- rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(Psi)))
  ZtZXe <- crossprod(Z, cbind(Z, X, y))

  score_inf <- limestest::score_psi(Z = Z, ZtZXe = ZtZXe, e = y, H = H,
                                    Psi0 = Psi / psi0, psi0 = psi0, finf = TRUE)
  test_stat <- score_inf$score[3]^2 / score_inf$finf[3, 3]

  pchisq(test_stat, 1, lower = F)
}

p_vals <- sapply(1:1000, one_sim)
print(mean(p_vals < 0.05))

