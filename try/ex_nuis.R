library(lme4)
library(Matrix)
#library(limestest)
library(rstiefel)
library(trust)
library(latex2exp)
source("R/functions.R")
library(profvis)

data(fev1)
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
            data = fev1, REML = FALSE)
X <- getME(fit, "X")
Z <- getME(fit, "Z")
y <- getME(fit, "y")
Lambda <- getME(fit, "Lambda")

psi0_hat <- getME(fit, "sigma")
Psi_hat <- psi0_hat^2 * tcrossprod(Lambda)

getME(fit, "GP") # groups pointer vector
tcrossprod(getME(fit, "Tlist")[[1]])
# getME(fit, "sigma") = attr(VarCorr(fit), "sc") = psi0_hat
# getME(fit, "sigma")^2 * tcrossprod(Lambda)[1:2,1:2]
# = as.matrix(Matrix::bdiag(VarCorr(fit)))
# = Psi1_hat <-  matrix(VarCorr(fit)$id, 2, 2)

#random_effects <- ranef(fit, condVar = TRUE)
#cov_matrices <- attr(random_effects[[1]], "postVar")[1]
H1 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 0, 0, 1),  2, 2))
H <- cbind(H1, H2, H3)
Psi1_hat <-  matrix(VarCorr(fit)$id, 2, 2)
Psi_hat <- Matrix::kronecker(Matrix::Diagonal(300), Psi1_hat)
psi0_hat <- attr(VarCorr(fit), "sc")^2
beta_hat <- fixef(fit)
e <- y - X %*% beta_hat
ZtZXe <- Matrix::crossprod(Z, cbind(Z, X, e))


S1 <- getME(fit, "sigma")^2 * tcrossprod(Lambda)
S2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(VarCorr(fit)$id, 2, 2))
all.equal(S1,S2)

# # parameters: pars (psi0, psi1, psi2), interest (psi3), X, Y, Z, H
# opt = optim(par = c(psi0_hat, Psi1_hat[1,1], Psi1_hat[1,2]), fn = neg_res, interest = Psi1_hat[2,2], X=X, Y=y, Z=Z, H=H, method = "BFGS")
# # parameters: Y, X, Z, H, Psi0, psi0
# result <- scoreandfinf(X=X, Y=y, Z = Z, H = H,
#                        pars = opt$par, interest = Psi1_hat[2,2])
# teststat <- t(result$score) %*% chol2inv(chol(result$finf)) %*% result$score
# teststat > qchisq(0.95, 1)



### data2
fit2 <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin)
getME(fit2, "Lind")
getME(fit2, "Lambda")
getME(fit, "sigma")^2 * tcrossprod(getME(fit2, "Lambda"))


### data3
data(Orthodont,package="nlme")
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
Orthodont$nsexage <- with(Orthodont, nsex*age)
fit3 <- lmer(distance ~ age + (age|Subject) + (0+nsex|Subject)+(0 + nsexage|Subject) , data=Orthodont) #
as.matrix(Matrix::bdiag(VarCorr(fit3)))

par <- 3 # par to test, psi1

score <- scoretest_par(object = fit2, par = par, order.par = 1, reml = F)
score

Psi0 <- tcrossprod(getME(fit3, "Lambda"))
gp <- getME(fit3, "Gp")   # total of q
# num of unique block length(getME(fit3, "Tlist"))
getME(fit3, "Tlist")

# as.matrix(Matrix::bdiag(VarCorr(fit3))) = tcrossprod(getME(fit3, "Lambda"))[1:2,1:2]*getME(fit3, "sigma")^2
tcrossprod(getME(fit3, "Tlist")[[1]]) # = tcrossprod(getME(fit3, "Lambda"))[1:2,1:2]
tcrossprod(getME(fit3, "Tlist")[[2]])
tcrossprod(getME(fit3, "Tlist")[[3]])
lind <- getME(fit3, "Lind")
getME(fit3, "theta")
View(getME(fit3, "Lambda"))

# construct H
gp <- getME(fit3, "Gp") # pointer to the beginning of each group, (gp[i+1]-gp[i])/dp[i] = number of blocks in each group
# dp <- NULL
# for (t in getME(fit3, "Tlist")) {
#   dp <- c(dp, dim(t)[1])
# }
dp <- sapply(getME(fit3, "Tlist"), dim)[1,]  # store the dimensions of blocks
H <- NULL
prevDim <- 0
for (i in 1:length(getME(fit3, "Tlist"))) {          # num of unique blocks is length(getME(fit3, "Tlist"))
  for (row in 1:dp[i]) {
    for (col in 1:row) {
      b <- sparseMatrix(row, col, dims = c(dp[i], dp[i]), x = 1, symmetric = T) # a single block
      block <- bdiag(rep(list(b), (gp[i+1]-gp[i])/dp[i]))                       # block diags of current group
      upcorner <- sparseMatrix(i = NULL, j = NULL, dims = c(prevDim, prevDim))  # complete dims
      downcorner <- sparseMatrix(i = NULL, j = NULL, dims = c(gp[length(gp)]-prevDim-nrow(block), gp[length(gp)]-prevDim-nrow(block)))
      H <- cbind(H, bdiag(upcorner, block, downcorner))
    }
  }
  prevDim <- prevDim + nrow(block)
}

######################### simpler example
set.seed(1)
N <- 4000
psi1 <- c(0.0001, 0.001, 0.01, 0.02, 0.03, 0.05, 0.1) #, 0.2, 0.5, 1, 10)
psi2 <- 0.01
psi0 <- 1
num_cluster <- 20
num_ind <- 5
X <- cbind(1, sample(c(0, 1), num_ind * num_cluster, replace = T))
Z <- Matrix::bdiag(lapply(1:num_cluster, function(ii)X[((ii - 1) * num_ind + 1):(ii * num_ind), ]))
beta <- rep(0, 2)
ZtZXe <- Matrix::crossprod(Z, cbind(Z, X))
Hlist <- list(matrix(c(1, 0, 0, 0), 2, 2), matrix(c(0, 0, 0, 1), 2, 2))

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
    loglik <- loglik_psi(Z, ZtZXe, e, H, Psi0, Psilist$psi0, loglik = TRUE,
                         score = TRUE, finf = TRUE, expected = TRUE)
    finv <- Matrix::solve(loglik$finf)[2,2]
    test[k, i] <- loglik$score[2]^2 * finv
    if (k %% 10 == 0) {
      print(paste(i, ",", k,"th finished at", Sys.time()-starttime))
    }
  }
}

save(test, file = "psi1Test_7.Rdata")
phat <- test < qchisq(0.95, 1)
coverage <- apply(phat, 2, mean)
mcerr <- sqrt(coverage*(1-coverage)/N) * 1.96

plot(psi1, coverage, xlab = TeX("$psi1$"), ylab = "Coverage Probability", pch = 20, ylim = c(0.85,1))
lines(psi1, coverage, lty=2)
abline(h = 0.95, col = "red")
arrows(psi1, coverage-mcerr, psi1, coverage+mcerr, angle=90, code=3, length = 0.05)






############### put together
####### example
data(fev1)
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
            data = fev1, REML = FALSE)
# Hlist <- list(matrix(c(1, 0, 0, 0), 2, 2), matrix(c(0, 1, 1, 0), 2, 2), matrix(c(0, 0, 0, 1), 2, 2))
par <- 5 # par to test, psi1

# Hlist: default first mat is par of interest, others are nuisance pars, RESTRICTION: sigmas are linear and additive
score <- scoretest_par(object = fit, par = par, order.par = 1)
score






#### check H is same
object <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (age|id),
              data = fev1, REML = FALSE)
gp <- getME(object, "Gp")        # pointer to the beginning of each group, (gp[i+1]-gp[i])/dp[i] = number of blocks in each group
dp <- sapply(getME(object, "Tlist"), dim)[1,]         # store the dimensions of blocks
H <- NULL
prevDim <- 0
for (i in 1:length(getME(object, "Tlist"))) {          # num of unique blocks is length(getME(object, "Tlist"))
  for (row in 1:dp[i]) {
    for (col in 1:row) {
      b <- sparseMatrix(row, col, dims = c(dp[i], dp[i]), x = 1, symmetric = T) # a single block
      block <- bdiag(rep(list(b), (gp[i+1]-gp[i])/dp[i]))                       # block diags of current group
      upcorner <- sparseMatrix(i = NULL, j = NULL, dims = c(prevDim, prevDim))  # complete dims
      downcorner <- sparseMatrix(i = NULL, j = NULL, dims = c(gp[length(gp)]-prevDim-nrow(block), gp[length(gp)]-prevDim-nrow(block)))
      H <- cbind(H, bdiag(upcorner, block, downcorner))
      print(block)
    }
  }
  prevDim <- prevDim + nrow(block)
}
H1 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(1, 0, 0, 0), 2, 2))
H2 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 1, 1, 0), 2, 2))
H3 <- Matrix::kronecker(Matrix::Diagonal(300), matrix(c(0, 0, 0, 1), 2, 2))
Hpre <- cbind(H1, H2, H3)
all.equal(H, Hpre)
#######


psi0_test <- attr(VarCorr(fit), "sc")^2
Psi1_test <- matrix(VarCorr(fit)$id, 2, 2)
Psi_test <-  Matrix::kronecker(Matrix::Diagonal(300), Psi1_test)
ll <- loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
           Psi0 = Psi_test / psi0_test, psi0 = psi0_test, loglik = T,
           finf = T, expected=T)
(t(ll$score) %*% solve(ll$finf) %*% ll$score)


estedvar <- sparseMatrix(i = NULL, j = NULL, dims = c(nrow(H), nrow(H)))
for (i in 1:(ncol(H)/nrow(H))) {
  estedvar <- estedvar + estimates[i+1] * H[,((i-1)*nrow(H)+1):(i*nrow(H))]
}
all.equal(tcrossprod(getME(fit3, "Lambda")) *getME(fit3, "sigma")^2, estedvar)
