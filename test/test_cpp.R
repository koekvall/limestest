library(Matrix)
library(doParallel)
set.seed(335)
num_reps <- 1e3
################################################################################
## INDEPENDENT CLUSTERS ########################################################
################################################################################
num_cluster <- 30
num_ind <- 5
Z <- Matrix::bdiag(replicate(num_cluster, cbind(1, rnorm(num_ind)),
                             simplify = FALSE))
X <- matrix(rnorm(nrow(Z) * 20), nrow = nrow(Z), ncol = 20)
b <- runif(ncol(X))
Psi1 <- matrix(c(1, 0, 0, 1), 2, 2)
psi0 <- 0.5
Psi <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1)
Psi0 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1 / psi0)
H1 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 0, 0, 1), 2, 2))
H <- as(cbind(H1, H2, H3), "sparseMatrix")

Psi1_eig <- eigen(Psi1)
R1 <- tcrossprod(diag(sqrt(Psi1_eig$values), ncol(Psi1)), Psi1_eig$vectors)
R <-  as(Matrix::kronecker(Matrix::Diagonal(num_cluster), R1), "sparseMatrix")

XtX = crossprod(X)
XtZ = crossprod(X, Z)
ZtZ = crossprod(Z)
Xb <- X %*% b
y <- Xb + rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))
XtY <- crossprod(X, y)
stuff_REML <- limestest::res_ll(XtX =XtX,
                                XtY = XtY,
                                XtZ = XtZ,
                                ZtZ = ZtZ,
                                YtZ = crossprod(y, Z),
                                Y = y,
                                X = X,
                                Z = Z,
                                H = H,
                                Psi0 = Psi0,
                                psi0 = psi0,
                                score = TRUE,
                                finf = TRUE,
                                lik = TRUE)
stuff_RML_cpp <- limestest:::res_ll_cpp(XtX = XtX,
                                        XtY = as.matrix(XtY),
                                        XtZ = as.matrix(XtZ),
                                        ZtZ = as(ZtZ, "dgCMatrix"))
