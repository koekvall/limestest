
# # function with psi1 as interest, (psi0, psi2) optimize by trust
# est_psi <- function(Z, ZtZXe, e, H, psi1) #loglik = TRUE,score = TRUE, finf = TRUE, expected = FALSE
# {
#   # Define dimensions
#   n <- length(e)
#   q <- ncol(Z)
#   r <- ncol(H) / q # Assumes H = [H_1, ... , H_r], where H_j is q by q
#   p <- ncol(ZtZXe) - q - 1
#   m <- q / 2  # check
#
#   objfun <- function(pars) { # pars = c(psi0, psi2)
#
#     Psi1 <- matrix(c(psi1,0,0,pars[2]), nrow = 2)
#     psi0 <- pars[1]
#     Psi0 <- Matrix::kronecker(Matrix::Diagonal(m), Psi1 / psi0)
#
#     # loglikelihood to return
#     ll <- NA
#     # Score vector for trust, assuming for psi0, par of interest, others..
#     s_psi <- rep(NA, r + 1)
#     # fisher information
#     I_psi <- matrix(NA, r + 1, r + 1)
#
#     # Pre-compute Psi0ZtZ (columns 1:q), Pzi0ZtX (columns (q + 1):(q + p)),
#     # and Psi0Zte (column q + p + 1)
#     A <- Matrix::crossprod(Psi0, ZtZXe)
#
#     eig <- eigen(A[, 1:q] + Matrix::Diagonal(q))
#
#     ##### check if A + I is p.d.
#     if (eig$values[q] <= 0) {
#       ll <- -Inf
#       # print(paste("Force ll to -Inf"))
#     } else {
#       # Add loglik term before overwriting
#       ll <- -0.5 * Matrix::determinant(A[, 1:q] + Matrix::Diagonal(q))$modulus[1] -
#         0.5 * n * log(psi0)
#
#       ##### Matrix denoted M in manuscript is A[, 1:q],
#       ##### main step that needs to check eig$value > 0
#       A <- Matrix::solve(A[, 1:q] + Matrix::Diagonal(q), A, sparse = TRUE)
#
#       # Score for error variance psi_0
#       # NB: REPLACE e by Sigma^{-1}e
#       e_save <- e
#       e <- (1 / psi0) * (e - Z %*% A[, q + p + 1]) # = Sigma^{-1}e
#       ll <- ll  - 0.5 * sum(e * e_save)
#
#       trace_M <- sum(Matrix::diag(A[, 1:q]))
#       s_psi[1] <- 0.5 * sum(e^2) - (0.5 / psi0) * (n - trace_M)
#
#       # Use recycling to compute v'H_i v for all Hi
#       v <- as(Matrix::crossprod(Z, e), "sparseVector") # sparse matrix does not recycle
#       w <- Matrix::crossprod(v, H) # Used later if !expected
#       s_psi[-1] <- 0.5 * colSums(matrix(as.vector(w * v), nrow = q))
#
#       # B = Z'Z (M - I_q) in paper notation
#       B <- A[, 1:q]
#       Matrix::diag(B) <- Matrix::diag(B) - 1
#       B <- Matrix::crossprod(ZtZXe[, 1:q], B)
#
#       # if finf
#       H <- B %*% H
#       I_psi[1, 1] <- (0.5 / psi0^2) * (n - 2 * trace_M +
#                                          sum(Matrix::t(A[, 1:q]) * A[, 1:q]))
#       # Subtract identity matrix from M
#       Matrix::diag(A[, 1:q]) <- Matrix::diag(A[, 1:q]) - 1
#       D <- matrix(Matrix::colSums(as.vector(A[, 1:q])  * as.matrix(H)), nrow = q)
#       I_psi[1, -1] <- (0.5 / psi0^2) * Matrix::colSums(D)
#
#       for(ii in 1:r){
#         first_idx <- ((ii - 1) * q + 1):(ii * q)
#         s_psi[1 + ii] <- s_psi[1 + ii] + (0.5 / psi0) * sum(Matrix::diag(H[, first_idx]))
#         for(jj in ii:r){
#           second_idx <- ((jj - 1) * q + 1):(jj * q)
#           I_psi[ii + 1, jj + 1] <- (0.5 / psi0^2) * sum(Matrix::t(H[, second_idx]) * H[, first_idx])
#         }
#       }
#       # !expect, observed I, or Hessian of -ll
#       I_psi <- -I_psi
#       # u = Sigma^{-2}e. Some calculations could be saved from before
#       u <- (1 / psi0^2) * (e_save + Z %*% (-2 * A[, q + p + 1] +
#                                              (A[, 1:q] + Matrix::Diagonal(q, 1)) %*% A[, q + p + 1]))
#       I_psi[1, 1] <- I_psi[1, 1] + sum(e * u)
#       v <- as.vector(Matrix::crossprod(Z, u)) # = Z' Sigma^{-2}e
#       I_psi[1, -1] <- I_psi[1, -1] + colSums(matrix(v * w, ncol = r))
#
#       for(ii in 1:r){
#         first_idx <- ((ii - 1) * q + 1):(ii * q)
#         for(jj in ii:r){
#           second_idx <- ((jj - 1) * q + 1):(jj * q)
#           I_psi[ii + 1, jj + 1] <- I_psi[ii + 1, jj + 1] - (1 / psi0) * sum(crossprod(w[first_idx], B) * w[second_idx])
#         }
#       }
#       I_psi <- as.matrix(-Matrix::forceSymmetric(I_psi, uplo = "U"))
#       # psi1 (par of interest) not need to calculate derivatives
#       s_psi <- s_psi[-2]
#       I_psi <- I_psi[-2,-2]
#     }
#
#     result <- list(value = ll, gradient = s_psi, hessian = -I_psi)
#     result
#   }
#   obj <- trust(objfun, parinit = c(1,0.01), rinit = 1, rmax = 10, minimize = F)
#   # if (obj$converged == FALSE) {print(obj)}
#   return(list(psi0 = obj$argument[1], Psi1 = matrix(c(psi1,0,0,obj$argument[2]), nrow = 2)))
# }

## outer function
# scoretest_par <- function(object, Hlist = NULL, par, order.par = 1) {
#   # num_cluster <- dim(ranef(object, condVar = TRUE)[[1]])[1]
#   X <- getME(object, "X")
#   Z <- getME(object, "Z")
#   y <- getME(object, "y")
#   beta_hat <- fixef(object)
#   e <- y - X %*% beta_hat
#   ZtZXe <- Matrix::crossprod(Z, cbind(Z, X, e))
#
#   if (is.null(Hlist)) {
#     H <- covInfo(object)$H
#     inits <- covInfo(object)$estimates[-(order.par+1)]
#   } else{}      # could user specify Hlist and inits?
#
#   Psilist <- est_par(Z, ZtZXe, e, H, par, order.par = order.par, inits = inits)
#
#   #Psi0 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psilist$Psi1 / Psilist$psi0)
#   loglik <- loglik_psi(Z, ZtZXe, e, H, Psilist$Psi0, Psilist$psi0, loglik = TRUE,
#                        score = TRUE, finf = TRUE, expected = TRUE)
#   finv <- Matrix::solve(loglik$finf)[2,2]
#   test <- loglik$score[2]^2 * finv
#   return(test)
# }




library(profvis)


####### generalized estimates of nuisance par
# Z, ZtZXe, e, H, par, order.par = order.par, inits = inits
est_par <- function(Z, ZtZXe, e, H, psi1, order.par, inits = NULL) #, tolr = 1e-4,,,,loglik = TRUE,score = TRUE, finf = TRUE, expected = FALSE
{
  # Define dimensions
  n <- length(e)
  q <- ncol(Z)
  r <- ncol(H) / q  # Assumes H = [H_1, ... , H_r], where H_j is q by q, r >= 2 to get here
  p <- ncol(ZtZXe) - q - 1
  rinds <- (1:r)[-order.par]
  if (is.null(inits)) {inits = c(1,rep(0.01, r-1))}

  # optimizer for nuisance pars
  objfun <- function(pars) { # pars = c(psi0, psi2,...), number of pars = r-1+1

    psi0 <- pars[1] # default pars[1] is the error scale

    {# Psi1 <- psi1 * Hlist[[1]]                                ############ need to change, not only one block
    # if (r >= 2) { # there are nuisance pars apart from psi0
    #   for (i in 2:r) {
    #     Psi1 <- Psi1 + pars[i] * Hlist[[i]]
    #   }
    # }
    # Psi0pre <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psi1 / psi0)
    }
    {# Psi0pre <- psi1 * H[,1:q]
    # if (r > 1) {          # there are nuisance pars other than psi0 ############ change order
    #   for (i in 2:r) {
    #     Psi0pre <- Psi0pre + pars[i] * H[,((i-1)*q+1):(i*q)]
    #   }
    # }
    }
    Psi0 <- psi1 * H[,((order.par-1)*q+1):(order.par*q)]
    for (i in 1:length(rinds)) {
      Psi0 <- Psi0 + pars[i+1] * H[,((rinds[i]-1)*q+1):(rinds[i]*q)]        # pars: 2:r <-> H:(1:r)[-order.par]
    }
    print(paste("pars: ", pars))

    # begin ...
    # loglikelihood to return
    ll <- NA
    # Score vector for trust, assuming for psi0, par of interest, others..
    s_psi <- rep(NA, r + 1)
    # fisher information
    I_psi <- matrix(NA, r + 1, r + 1)

    # Pre-compute Psi0ZtZ (columns 1:q), Pzi0ZtX (columns (q + 1):(q + p)),
    # and Psi0Zte (column q + p + 1)
    A <- Matrix::crossprod(Psi0, ZtZXe)

    eig <- eigen(A[, 1:q] + Matrix::Diagonal(q))

    ##### check if A + I is p.d.
    if (Im(eig$values[q]) != 0 || Re(eig$values[q]) <= 0) {
      ll <- -Inf
      # print(paste("Force ll to -Inf"))
    } else {
      # Add loglik term before overwriting
      print(paste("psi0= ", psi0))
      ll <- -0.5 * Matrix::determinant(A[, 1:q] + Matrix::Diagonal(q))$modulus[1] -
        0.5 * n * log(psi0)

      ##### Matrix denoted M in manuscript is A[, 1:q],
      ##### main step that needs to check eig$value > 0
      A <- Matrix::solve(A[, 1:q] + Matrix::Diagonal(q), A, sparse = TRUE)

      # Score for error variance psi_0
      # NB: REPLACE e by Sigma^{-1}e
      e_save <- e
      e <- (1 / psi0) * (e - Z %*% A[, q + p + 1]) # = Sigma^{-1}e
      ll <- ll  - 0.5 * sum(e * e_save)

      trace_M <- sum(Matrix::diag(A[, 1:q]))
      s_psi[1] <- 0.5 * sum(e^2) - (0.5 / psi0) * (n - trace_M)

      # Use recycling to compute v'H_i v for all Hi
      v <- as.vector(Matrix::crossprod(Z, e)) # sparse matrix does not recycle
      w <- as.vector(Matrix::crossprod(v, H)) # Used later if !expected
      s_psi[-1] <- 0.5 * colSums(matrix(as.vector(w * v), nrow = q))

      # B = Z'Z (M - I_q) in paper notation
      B <- A[, 1:q]
      Matrix::diag(B) <- Matrix::diag(B) - 1
      B <- Matrix::crossprod(ZtZXe[, 1:q], B)

      # if finf
      H <- B %*% H
      I_psi[1, 1] <- (0.5 / psi0^2) * (n - 2 * trace_M +
                                         sum(Matrix::t(A[, 1:q]) * A[, 1:q]))
      # Subtract identity matrix from M
      Matrix::diag(A[, 1:q]) <- Matrix::diag(A[, 1:q]) - 1
      D <- matrix(Matrix::colSums(as.vector(A[, 1:q])  * as.matrix(H)), nrow = q)
      I_psi[1, -1] <- (0.5 / psi0^2) * Matrix::colSums(D)

      for(ii in 1:r){
        first_idx <- ((ii - 1) * q + 1):(ii * q)
        s_psi[1 + ii] <- s_psi[1 + ii] + (0.5 / psi0) * sum(Matrix::diag(H[, first_idx]))
        for(jj in ii:r){
          second_idx <- ((jj - 1) * q + 1):(jj * q)
          I_psi[ii + 1, jj + 1] <- (0.5 / psi0^2) * sum(Matrix::t(H[, second_idx]) * H[, first_idx])
        }
      }
      # !expect, observed I, or Hessian of -ll
      I_psi <- -I_psi
      # u = Sigma^{-2}e. Some calculations could be saved from before
      u <- (1 / psi0^2) * (e_save + Z %*% (-2 * A[, q + p + 1] +
                                             (A[, 1:q] + Matrix::Diagonal(q, 1)) %*% A[, q + p + 1]))
      I_psi[1, 1] <- I_psi[1, 1] + sum(e * u)
      v <- as.vector(Matrix::crossprod(Z, u)) # = Z' Sigma^{-2}e
      I_psi[1, -1] <- I_psi[1, -1] + colSums(matrix(v * w, ncol = r))

      for(ii in 1:r){
        first_idx <- ((ii - 1) * q + 1):(ii * q)
        for(jj in ii:r){
          second_idx <- ((jj - 1) * q + 1):(jj * q)
          I_psi[ii + 1, jj + 1] <- I_psi[ii + 1, jj + 1] - (1 / psi0) * sum(crossprod(w[first_idx], B) * w[second_idx])
        }
      }
      I_psi <- as.matrix(Matrix::forceSymmetric(I_psi, uplo = "U"))#as.matrix(-Matrix::forceSymmetric(I_psi, uplo = "U"))
      # psi1 (par of interest) not need to calculate derivatives
      s_psi <- s_psi[-2]
      I_psi <- I_psi[-2,-2]
    }
    result <- list(value = ll, gradient = s_psi, hessian = -I_psi)
    result
  }
  obj <- trust(objfun, parinit = inits, rinit = 1, rmax = 10, minimize = F)
  {############ need to change, not only one block
  # Psi1_opt <- psi1 * Hlist[[1]]
  # for (i in 2:r) {
  #   Psi1_opt <- Psi1_opt + obj$argument[i] * Hlist[[i]]
  # }
  ############ need to change, not only one block, return Psi0
  }
  {# Psi0_opt <- psi1 * H[,1:q]
  # if (r > 1) {          ############# change order
  #   for (i in 2:r) {
  #     Psi0_opt <- Psi0_opt + obj$argument[i] * H[,((i-1)*q+1):(i*q)]
  #   }
  # }
    }
  Psi0_opt <- psi1 * H[,((order.par-1)*q+1):(order.par*q)]
  for (i in 1:length(rinds)) {
    Psi0_opt <- Psi0_opt + obj$argument[i+1] * H[,((rinds[i]-1)*q+1):(rinds[i]*q)]        # pars: 2:r <-> H:(1:r)[-order.par]
  }
  Psi0_opt <- Psi0_opt / obj$argument[1]
  return(list(psi0 = obj$argument[1], Psi0 = Psi0_opt))
}


## construct H
covInfo <- function(object) {
  gp <- getME(object, "Gp")        # pointer to the beginning of each group, (gp[i+1]-gp[i])/dp[i] = number of blocks in each group
  # dp <- sapply(getME(object, "Tlist"), dim)[1,]         # store the dimensions of blocks
  H <- NULL
  prevDim <- 0
  estimates <- getME(object, "sigma")^2
  for (i in 1:length(getME(object, "Tlist"))) {          # num of unique blocks is length(getME(object, "Tlist"))
    vcor <- tcrossprod(getME(fit3, "Tlist")[[i]]) #* getME(object, "sigma")^2
    d <- dim(vcor)
    for (row in 1:d[1]) {
      for (col in 1:row) {
        estimates <- c(estimates, vcor[row, col])

        b <- list(sparseMatrix(row, col, dims = d, x = 1, symmetric = T)) # a single block)
        block <- bdiag(rep(b, (gp[i+1]-gp[i])/d[1]))                       # block diags of current group
        upcorner <- sparseMatrix(i = NULL, j = NULL, dims = c(prevDim, prevDim))  # complete dims
        downcorner <- sparseMatrix(i = NULL, j = NULL, dims = c(gp[length(gp)]-prevDim-nrow(block), gp[length(gp)]-prevDim-nrow(block)))
        H <- cbind(H, bdiag(upcorner, block, downcorner))
      }
    }
    prevDim <- prevDim + nrow(block)
  }
  return(list("H" = H, "estimates" = estimates))
}





###################################
############### reml version
###################################
#, lik = TRUE, score = FALSE, finf = FALSE, Psi0, psi0
res_est <- function(XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H, psi1, order.par, inits = NULL) #, tolr = 1e-4
{
  # Define dimensions
  n <- length(Y)
  q <- ncol(Z)
  r <- ncol(H) / q # Assumes H = [H_1, ... , H_r], where H_j is q by q
  p <- ncol(XtX)
  rinds <- (1:r)[-order.par]
  if (is.null(inits)) {inits = c(1,rep(0.01, r-1))}

  objfun <- function(pars) { # pars = c(psi0, psi2,...), number of pars = r-1+1

    psi0 <- pars[1] # default pars[1] is the error scale
    Psi0 <- psi1 * H[,((order.par-1)*q+1):(order.par*q)]
    for (i in 1:length(rinds)) {
      Psi0 <- Psi0 + pars[i+1] * H[,((rinds[i]-1)*q+1):(rinds[i]*q)]        # pars: 2:r <-> H:(1:r)[-order.par]
    }

    # Loglikelihood to return
    ll <- NA

    # Score vector to return
    s_psi <- rep(NA, r + 1)

    # Fisher information to return
    I_psi <- matrix(NA, r + 1, r + 1)
    # print(paste("dims: Psi0 ", dim(Psi0), "q ",q, "dims: ZtZ", dim(ZtZ)))
    # Pre-compute A = (I_q + Psi0 Z'Z)^{-1} Psi0
    A <- Matrix::crossprod(Psi0, ZtZ) + Matrix::Diagonal(q) # q x q storage

    # Add likelihood term before overwriting, if(lik)
    ll <- Matrix::determinant(A, logarithm = TRUE)$modulus

    A <- Matrix::solve(A, Psi0)
    B <- XtZ %*% A # q x q

    # Create XtSiX
    U <- Matrix::forceSymmetric((1 / psi0) * (XtX - Matrix::tcrossprod(B, XtZ))) # p x p, now XtSiX
    U <- try(Matrix::chol(U)) # replace XtSiZ by its Cholesky root

    if(inherits(U,"try-error")){
      return(list("ll" = -Inf, "score" = s_psi, "finf" = I_psi, "beta" = rep(NA, p),
                  "I_b_inv_chol" = matrix(NA, p, p)))
    }

    # Create XtSiY for use in beta_tilde
    beta_tilde <- (1/ psi0) * (XtY - XtZ %*% Matrix::tcrossprod(A, YtZ)) # p x 1
    beta_tilde <- chol_solve(U, beta_tilde)

    # Replace Y by residuals
    Y <- Y - X %*% beta_tilde

    # n x 1 vector for storing \Sigma^{-1}e
    a <- (1 / psi0) * (Y - Z %*% (A %*% Matrix::crossprod(Z, Y))) # n x 1

    # if lik
    ll <- ll + 2 * sum(log(Matrix::diag(U)))
    ll <- ll + sum(Y * a) + n * log(2 * pi * psi0)
    ll <- -0.5 * ll


    # if score
    # Stochastic part of restricted score for psi
    s_psi[1] <- 0.5 * sum(a^2)
    v <- as.vector(Matrix::crossprod(Z, a)) # q x 1 vector storage
    s_psi[-1] <- 0.5 * colSums(matrix(as.vector(Matrix::crossprod(v, H)) * v,
                                      nrow = q))

    #############################################################################
    ## NOTHING BELOW SHOULD DEPEND ON Y.
    #############################################################################

    # if(finf & score){
    A <- Matrix::tcrossprod(A, ZtZ) # q x q, called M in manuscript

    s_psi[1] <- s_psi[1] - (0.5 / psi0) * n + (0.5 / psi0) * sum(Matrix::diag(A))

    I_psi[1, 1] <- (0.5 / psi0^2) * (n - 2 * sum(Matrix::diag(A)) +
                                       sum(Matrix::t(A) * A))

    E <- Matrix::crossprod(ZtZ, A) # q x q storage
    D <- XtZ %*% A # p x q storage
    XtSiZ <- (1 / psi0) * (XtZ - D) # p x q
    XtSi2Z <- (1 / psi0)^2 * (XtZ - 2 * D + D %*% A) # p x q


    C <- Matrix::tcrossprod(B, XtZ) # p x p storage, here XtZA ZtX
    G <- B %*% Matrix::tcrossprod(ZtZ, B) # p x q, here XtZA ZtZ AtZtX
    XtSi2X <- (1 / psi0)^2 * (XtX - 2 * C + G) # p x p
    XtSi3X <- (1 / psi0^3) * (XtX - 3 * C + 3 * G -
                                D %*% tcrossprod(A, B)) # p x p
    C <- chol_solve(U, XtSi2X)

    I_psi[1, 1] <- I_psi[1, 1] + 0.5 * sum(C * Matrix::t(C))

    s_psi[1] <- s_psi[1] + 0.5 * sum(Matrix::diag(C))


    I_psi[1, 1] <- I_psi[1, 1] - sum(Matrix::diag(chol_solve(U, XtSi3X)))
    # A (q x q), G (p x q) ARE FREE
    A <- (1 / psi0)^2 * (ZtZ - 2 * E + E %*% A) # ZtSi2Z right now
    E <-  (1/ psi0) * (ZtZ - E) # Now holds ZtSiZ
    D <- chol_solve(U, XtSiZ)
    A <- A - 2 * Matrix::crossprod(D, XtSi2Z) + Matrix::crossprod(XtSiZ, C %*% D)

    I_psi[-1, 1] <- 0.5 * colSums(matrix(Matrix::colSums(as.vector(A) * H), nrow = q))
    s_psi[-1] <- s_psi[-1] - 0.5 * colSums(matrix(Matrix::colSums(
      as.vector(E - Matrix::crossprod(XtSiZ, D)) * H), nrow = q))

    H2 <- Matrix::crossprod(XtSiZ, D %*% H) # Storage can be avoided by
    # muliply in loop
    # Has to come after H2 since H is overwritten
    H <- Matrix::crossprod(E, H) # = ZtSiZ %*% H
    for(ii in 1:r){
      idx1 <- ((ii - 1) * q + 1):(ii * q)
      for(jj in ii:r){
        idx2 <-  ((jj - 1) * q + 1):(jj * q)
        I_psi[jj + 1, ii + 1] <- 0.5 * sum(H[, idx1] * Matrix::t(H[, idx2])) -
          sum(H[, idx1] * Matrix::t(H2[, idx2])) + 0.5 *
          sum(H2[, idx1] * Matrix::t(H2[, idx2]))
      }
    }
    I_psi <- as.matrix(Matrix::forceSymmetric(I_psi, "L"))
    #I_psi <- as.matrix(Matrix::forceSymmetric(I_psi, uplo = "U"))#as.matrix(-Matrix::forceSymmetric(I_psi, uplo = "U"))
    # psi1 (par of interest) not need to calculate derivatives
    s_psi <- s_psi[-2]
    I_psi <- I_psi[-2,-2]

    result <- list(value = ll, gradient = s_psi, hessian = -I_psi)
    result
  }
  obj <- trust(objfun, parinit = inits, rinit = 1, rmax = 10, minimize = F)

  Psi0_opt <- psi1 * H[,((order.par-1)*q+1):(order.par*q)]
  for (i in 1:length(rinds)) {
    Psi0_opt <- Psi0_opt + obj$argument[i+1] * H[,((rinds[i]-1)*q+1):(rinds[i]*q)]        # pars: 2:r <-> H:(1:r)[-order.par]
  }
  Psi0_opt <- Psi0_opt / obj$argument[1]
  return(list(psi0 = obj$argument[1], Psi0 = Psi0_opt))
}

# XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H, psi1, order.par, inits = NULL, tolr = 1e-4
scoretest_par <- function(object, Hlist = NULL, par, order.par = 1, reml = T) {
  # num_cluster <- dim(ranef(object, condVar = TRUE)[[1]])[1]
  X <- getME(object, "X")
  Z <- getME(object, "Z")
  y <- getME(object, "y")

  covInfo <- covInfo(object)
  if (is.null(Hlist)) {
    H <- covInfo$H
  }       # could user specify Hlist and inits?
  inits <- covInfo$estimates[-(order.par+1)]

  if (reml) {
    XtX <- crossprod(X)
    XtY <- crossprod(X,y)
    XtZ <- crossprod(X,Z)
    ZtZ <- crossprod(Z)
    YtZ <- crossprod(y,Z)
    Psilist <- res_est(XtX, XtY, XtZ, ZtZ, YtZ, y, X, Z, H, psi1 = par, order.par = order.par, inits = inits)
    loglik <- res_ll(XtX, XtY, XtZ, ZtZ, YtZ, y, X, Z, H, Psi0 = Psilist$Psi0, psi0 = Psilist$psi0,
                     lik = TRUE, score = TRUE, finf = TRUE)
  } else {
    beta_hat <- fixef(object)
    e <- y - X %*% beta_hat
    ZtZXe <- Matrix::crossprod(Z, cbind(Z, X, e))
    Psilist <- est_par(Z, ZtZXe, e, H, par, order.par = order.par, inits = inits)
    #Psi0 <- Matrix::kronecker(Matrix::Diagonal(num_cluster), Psilist$Psi1 / Psilist$psi0)
    loglik <- loglik_psi(Z, ZtZXe, e, H, Psilist$Psi0, Psilist$psi0, loglik = TRUE,
                         score = TRUE, finf = TRUE, expected = TRUE)
  }
  finv <- Matrix::solve(loglik$finf)[2,2]
  test <- loglik$score[2]^2 * finv
  return(test)
}
