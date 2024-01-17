library(lme4)
library(Matrix)
#library(limestest)
library(rstiefel)
library(trust)

data(fev1)
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
            data = fev1, REML = FALSE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
y <- getME(fit, "y")

#random_effects <- ranef(fit, condVar = TRUE)
#cov_matrices <- attr(random_effects[[1]], "postVar")
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
             score = TRUE, finf = TRUE, expected = FALSE)$ll
}

# # Test derivatives at MLE
# psi0_test <- psi0_hat
# Psi1_test <- Psi_hat[1:2, 1:2]
# Psi_test <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_test)




# parameters: pars (psi0, psi1, psi2), interest (psi3), X, Y, Z, H
opt = optim(par = c(psi0_hat, Psi1_hat[1,1], Psi1_hat[1,2]), fn = neg_res, interest = Psi1_hat[2,2], X=X, Y=y, Z=Z, H=H, method = "BFGS")
# parameters: Y, X, Z, H, Psi0, psi0
result <- scoreandfinf(X=X, Y=y, Z = Z, H = H,
                       pars = opt$par, interest = Psi1_hat[2,2])
teststat <- t(result$score) %*% chol2inv(chol(result$finf)) %*% result$score
teststat > qchisq(0.95, 1)



## use trust for regular likelihood, beta=0, simulation_intro.R with only diagonal value



# simulation:

n <- 600
r <- 3
p <- 1
q1 <- 2
m <- 200 # dim(Z) = n*(mq1), length(U)=m*q1
psi0 <- 1
Psi1 <- matrix(rnorm(4), nrow = 2)
Psi1 <- Psi1 %*% Psi1
Psi <- Matrix::kronecker(Matrix::Diagonal(m), Psi1)
X <- matrix(1, nrow = n)
H1 <- Matrix::kronecker(Matrix::Diagonal(m), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(m), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(m), matrix(c(0, 0, 0, 1), 2, 2))
H <- cbind(H1, H2, H3)

Z <- list()
for (i in 1:m) {
  Z[[i]] <- matrix(rnorm(6), nrow = 3)
}
Z <- bdiag(Z)
eigen(Psi)$values
eigen(Z %*% Psi %*% t(Z))$values
eigen(Z %*% Psi %*% t(Z) + diag(psi0, nrow = n))$values
y <- mvrnorm(mu = 0, Sigma = Z %*% Psi %*% t(Z) + diag(psi0, nrow = n))



## try
testfun <- function() {
  a <- 2
  objfun <- function(x) {
    m <- matrix(c(x[1],0,0,x[2]), nrow = 2)
    m <- Matrix::kronecker(Matrix::Diagonal(10), m)
    f <- sum(diag(crossprod(m)))
    g <- c(a*10*x[1],a*10*x[2])
    h <- matrix(c(a*10,0,0,a*10), nrow = 2)
    list(value = f, gradient = g, hessian = h)
  }
  obj <- trust(objfun, parinit = c(3,1), rinit = 1, rmax = 5)
  print(obj)
  xhat <- obj$argument
  return(xhat)
}
testfun()



### simpler example
set.seed(335)
N <- 1e2
psi1 <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05)
psi2 <- 0.01
psi0 <- 1
num_cluster <- 20
num_ind <- 5
X <- cbind(1, sample(c(0, 1), num_ind * num_cluster, replace = T))
Z <- Matrix::bdiag(lapply(1:num_cluster, function(ii)X[((ii - 1) * num_ind + 1):(ii * num_ind), ]))
beta <- rep(0, 2)
ZtZXe <- Matrix::crossprod(Z, cbind(Z, X))

test <- matrix(0, nrow = N, ncol = length(psi1))
starttime <- Sys.time()
for (i in 1:length(psi1)) {
  Psi1 <- matrix(c(psi1[i], 0,0, psi2), 2, 2)
  Psi <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1)
  Sigma <- diag(psi0, nrow(X)) + Z %*% tcrossprod(Psi, Z)
  M <- solve(Sigma, X) # S^{-1}X
  G <- solve(crossprod(X, M), t(M)) # (X'S^{-1}X)^{-1}X'S^{-1}

  Psi1_eig <- eigen(Psi1)
  R1 <- tcrossprod(diag(sqrt(Psi1_eig$values), ncol(Psi1)), Psi1_eig$vectors)
  R <-  as(Matrix::kronecker(Matrix::Diagonal(num_cluster), R1), "sparseMatrix")

  H1 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(1, 0, 0, 0), 2, 2))
  H2 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), matrix(c(0, 0, 0, 1), 2, 2))
  H <- as(cbind(H1, H2), "sparseMatrix")

  for (k in 1:N) {
    y <- X %*% beta + rnorm(nrow(Z), sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))
    beta_hat <- G %*% y
    e <- y - X %*% beta_hat
    ZtZXe <- cbind(ZtZXe, crossprod(Z, e))
    # estimation of nuisance parmeter when true psi1 is psi1[i]
    Psilist <- est_psi(Z, ZtZXe, e, H, psi1[i])
    Psi0 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psilist$Psi1 / Psilist$psi0)
    loglik <- loglik_psi(Z, ZtZXe, e, H, Psi0, psi0, loglik = TRUE,
                         score = TRUE, finf = TRUE, expected = TRUE)
    finv <- Matrix::solve(loglik$finf)[2,2]
    test[k, i] <- loglik$score[2]^2 / finv
    if (k %% 10 == 0) {
      print(paste(i, ",", k,"th finished at", Sys.time()-starttime))
    }
  }
}

phat <- test < qchisq(0.95, 1)
coverage <- apply(phat, 2, mean)
mcerr <- sqrt(coverage*(1-coverage)/N) * 1.96

plot(psi1, coverage, xlab = TeX("$psi1$"), ylab = "Coverage Probability", pch = 20)
lines(psi1, coverage, lty=2)
abline(h = 0.95, col = "red")
arrows(psi1, coverage-mcerr, psi1, coverage+mcerr, angle=90, code=3, length = 0.05)

# check
eigen(Z %*% Psi0 %*% t(Z) + diag(num_cluster*num_ind))$values
