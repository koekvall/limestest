set.seed(335)
library(Matrix)
num_reps <- 1e3
num_cluster <- 100
num_ind <- 10
Z <- Matrix::bdiag(replicate(num_cluster, cbind(1, rnorm(num_ind)),
                             simplify = FALSE))
X <- matrix(1, nrow = nrow(Z), ncol = 1)
Psi1 <- matrix(c(1, -0.99, -0.99, 1), 2, 2)
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
  ZtZXe <- crossprod(Z, cbind(Z, X, e))

  score_inf <- limestest::loglik_psi(Z = Z, ZtZXe = ZtZXe, e = y, H = H,
                                    Psi0 = Psi / psi0, psi0 = psi0, finf = TRUE)
  test_stat <- as.vector(crossprod(score_inf$score, solve(score_inf$finf, score_inf$score)))

  pchisq(test_stat, 4, lower = F)

}

p_vals <- sapply(1:num_reps, one_sim)
print(mean(p_vals < 0.05))
plot(quantile(p_vals, seq(0.01, 0.99, 0.01)))
