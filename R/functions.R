#' loglik_psi
#'
#' Calculates the log-likelihood, score vector and Fisher information matrix
#' for the covariance parameter vector \code{psi} in a linear mixed effects model.
#'
#' @param Z A matrix of fixed effects.
#' @param ZtZXe Pre-computed matrix with of Z'[Z, X, e]
#' @param e The residual vector.
#' @param H Matrix of derivatives of Psi with respect to elements of psi.
#'          Assumes H = [H_1, ... , H_r], where H_j is q by q.
#' @param Psi0 Coavariance matrix of random effects (Psi) divided by error
#'             variance psi0, with dimensions q by q.
#' @param psi0 The error variance.
#' @param loglik If \code{TRUE} (default), the log-likelihood will be calculated.
#' @param score If \code{TRUE} (default), the score vector will be calculated.
#' @param finf If \code{TRUE} (default), the Fisher information matrix will be calculated.
#' @param expected If \code{TRUE} (detault), return expected information;
#' otherwise observed.
#'
#' @return A list with components:
#' \item{ll}{The log-likelihood.}
#' \item{score}{The score vector.}
#' \item{finf}{The Fisher information matrix.}
#'
#' @import Matrix
#' @export
#' @useDynLib limestest, .registration=TRUE
loglik_psi <- function(Z, ZtZXe, e, H, Psi0, psi0, loglik = TRUE,
                       score = TRUE, finf = TRUE, expected = TRUE)
{
  # Define dimensions
  n <- length(e)
  q <- ncol(Psi0)
  r <- ncol(H) / q # Assumes H = [H_1, ... , H_r], where H_j is q by q
  p <- ncol(ZtZXe) - q - 1

  # loglikelihood to return
  ll <- NA

  # Score vector to return
  s_psi <- rep(NA, r + 1)

  # Fisher information to return
  I_psi <- matrix(NA, r + 1, r + 1)

  # Pre-compute Psi0ZtZ (columns 1:q), Pzi0ZtX (columns (q + 1):(q + p)),
  # and Psi0Zte (column q + p + 1)
  A <- Matrix::crossprod(Psi0, ZtZXe)

  # Add loglik term before overwriting
  if(loglik){
    ll <- -0.5 * Matrix::determinant(A[, 1:q] + Matrix::Diagonal(q))$modulus -
      0.5 * n * log(psi0)
  }

  # Matrix denoted M in manuscript is A[, 1:q]
  A <- Matrix::solve(A[, 1:q] + Matrix::Diagonal(q), A, sparse = TRUE)

  # Score for error variance psi_0
  # NB: REPLACE e by Sigma^{-1}e
  if(loglik | (finf & !expected)){
   e_save <- e
  }
  e <- (1 / psi0) * (e - Z %*% A[, q + p + 1]) # = Sigma^{-1}e
  if(loglik){
    ll <- ll  - 0.5 * sum(e * e_save)
  }

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

  if(!finf){
    # Compute -[ZtZ (I_q - M) * H_1, ..., ZtZ (I_q - M) * H_r] using
    # recycling, where * denotes elementwise multiplication
    # The "if" is because the calculation is a byproduct of a more expensive one
    # (B %*% H) done to get Fisher information
    H <- as.matrix(H) * as.vector(B)

    s_psi[-1] <- s_psi[-1] + (0.5 / psi0) * colSums(matrix(Matrix::colSums(H), nrow = q))
  } else{
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
    if(!expected){
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
          I_psi[ii + 1, jj + 1] <- I_psi[ii + 1, jj + 1] - (1 / psi0) *
            sum(crossprod(w[first_idx], B) * w[second_idx])
        }
      }
    }
  }
  I_psi <- Matrix::forceSymmetric(I_psi, uplo = "U")
  return(list("ll" = ll,  "score" = s_psi, "finf" = I_psi))
}

Psi_from_Hlist <- function(psi, Hlist)
{
  for(ii in seq_len(length(Hlist))){
    Hlist[[ii]] <- Hlist[[ii]] * psi[ii]
  }
  do.call(cbind, Hlist)
}

#' @export
uni_test_stat <- function(test_seq, test_idx, psi, psi0, Z, ZtZXe, e, Hlist)
{
  # first element of test_seq has to agree with psi[test_idx]
  m <- length(test_seq)
  test_stat <- rep(0, m)
  H <- do.call(cbind, Hlist)
  for(ii in seq_len(m)){
    if(ii > 1){
      # Update tested parameter
      psi[test_idx] <- test_seq[ii]
      # Do one-step Fisher scoring update (with previous Information)
      # for non-tested parameters
      Psi <- Psi_from_Hlist(psi, Hlist)
      s <- score_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                     Psi0 = Psi / psi0,
                     psi0 = psi0, finf = FALSE)
      psi[-test_idx] <- psi[-test_idx] +
        solve(score_inf$finf[-test_idx. -test_idx], s)
    }

    Psi <- Psi_from_Hlist(psi, Hlist)
    score_inf <- score_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                           Psi0 = Psi / psi0,
                           psi0 = psi0, finf = TRUE)

    eff_inf <- score_inf$finf[test_idx, test_idx] -
      sum(solve(score_inf$finf[-test_idx, -test_idx],
                score_inf$finf[-test_idx, test_idx]) *
            score_inf$finf[-test_idx, test_idx])

    test_stat[ii] <- score_inf$score[test_idx]^2 / eff_inf
  }
  test_stat
}


chol_solve <- function(U, b)
{
  Matrix::solve(U, Matrix::solve(Matrix::t(U), b))
}


#' Compute Restricted Likelihood, Score and Information
#'
#' Computes the restricted (residual) likelihood, score, and information matrix
#' for the variance parameter \code{psi} a linear mixed effects model
#'
#' @param XtX An n x p matrix of the crossproduct of the design matrix of fixed
#' effects (X) with itself.
#' @param XtY An n x 1 vector of the crossproduct of the design matrix of fixed
#' effects (X) with the response vector (Y).
#' @param XtZ An n x q matrix of the crossproduct of the design matrix of fixed
#' effects (X) with the design matrix of random effects (Z).
#' @param ZtZ A q x q matrix of the crossproduct of the design matrix of random
#' effects (Z) with itself.
#' @param YtZ A 1 x q matrix of the crossproduct of the response vector (Y) with
#' the design matrix of random effects (Z).
#' @param Y An n x 1 vector of the response variable.
#' @param X An n x p matrix of the design matrix of fixed effects.
#' @param Z An n x q matrix of the design matrix of random effects.
#' @param H A q x rq matrix, where H = [H_1, ..., H_r], with H_j being the
#' derivative of Psi with respect to psi_j
#' @param Psi0 The covariance matrix of the random effects (Psi) divided by the
#' error variance (psi0)
#' @param psi0 A scalar value of the error variance.
#' @param lik If \code{TRUE} (default), the log-likelihood will be computed.
#' @param score If \code{TRUE} (default), the score vector will be computed.
#' @param finf If \code{TRUE} (default), the Fisher information matrix will be
#' computed.
#'
#' @return A list with elements:
#' \describe{
#' \item{res_ll}{A scalar value of the restricted log-likelihood.}
#' \item{s_psi}{A (r + 1) x 1 vector of the restricted score of the variance parameters.}
#' \item{I_psi}{A (r + 1) x (r + 1) matrix of the restricted  information of the
#' variance parameters.}
#' }
#' @import Matrix
#' @export
res_ll <- function(XtX, XtY, XtZ, ZtZ, YtZ, Y, X, Z, H, Psi0, psi0, lik = TRUE, score = FALSE,
                   finf = FALSE)
{
  # Define dimensions
  n <- length(Y)
  q <- ncol(Psi0)
  r <- ncol(H) / q # Assumes H = [H_1, ... , H_r], where H_j is q by q
  p <- ncol(XtX)

  # Loglikelihood to return
  ll <- NA

  # Score vector to return
  s_psi <- rep(NA, r + 1)

  # Fisher information to return
  I_psi <- matrix(NA, r + 1, r + 1)
  # Pre-compute A = (I_q + Psi0 Z'Z)^{-1} Psi0
  A <- Matrix::crossprod(Psi0, ZtZ) + Matrix::Diagonal(q) # q x q storage

  # Add likelihood term before overwriting
  if(lik) ll <- Matrix::determinant(A, logarithm = TRUE)$modulus

  A <- Matrix::solve(A, Psi0)
  B <- XtZ %*% A # q x q

  # Create XtSiX
  U <- Matrix::forceSymmetric((1 / psi0) * (XtX - Matrix::tcrossprod(B, XtZ))) # p x p, now XtSiX
  U <- try(Matrix::chol(U)) # replace XtSiX by its Cholesky root
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

  if(lik){
    ll <- ll + 2 * sum(log(Matrix::diag(U)))
    ll <- ll + sum(Y * a) + n * log(2 * pi * psi0)
    ll <- -0.5 * ll
  }

  if(score){
    # Stochastic part of restricted score for psi
    s_psi[1] <- 0.5 * sum(a^2)
    v <- as.vector(Matrix::crossprod(Z, a)) # q x 1 vector storage
    s_psi[-1] <- 0.5 * colSums(matrix(as.vector(Matrix::crossprod(v, H)) * v,
                                      nrow = q))
  }
  #############################################################################
  ## NOTHING BELOW SHOULD DEPEND ON Y.
  #############################################################################

  if(finf & score){
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
    I_psi <- Matrix::forceSymmetric(I_psi, "L")
  } else if (score){
    A <- Matrix::tcrossprod(A, ZtZ) # q x q, called M in manuscript
    s_psi[1] <- s_psi[1] - (0.5 / psi0) * n + (0.5 / psi0) * sum(Matrix::diag(A))

    D <- XtZ %*% A # p x q storage

    C <- Matrix::tcrossprod(B, XtZ) # p x p storage
    C <- chol_solve(U, (1 / psi0)^2 * (XtX - 2 * C + B %*% Matrix::tcrossprod(ZtZ, B)))
    s_psi[1] <- s_psi[1] + 0.5 * sum(Matrix::diag(C))

    A <- Matrix::crossprod(ZtZ, A)
    A <- (1/ psi0) * (ZtZ - A)

    XtSiZ <- (1 / psi0) * (XtZ - D) # p x q
    D <- chol_solve(U, XtSiZ) # p x q
    v <- as.vector(A - Matrix::crossprod(XtSiZ, D)) # pq
    s_psi[-1] <- s_psi[-1] - 0.5 * colSums(matrix(Matrix::colSums(v * H), nrow = q))
  }
  return(list("ll" = ll[1], "score" = s_psi, "finf" = I_psi, "beta" = beta_tilde,
              "I_b_inv_chol" = U))
}

#' Get the covariance matrix of random effects
#'
#' Returns the covariance matrix of the random effects in a linear mixed model,
#' either as estimated by lme4::lmer or evaluated at a particular parameter value
#' supplied as an argument.
#'
#' @param lmerfit An lmerMod object from fitting a linear mixed model using lme4::lmer
#' @param psi Optional vector with covariance parameters (see details)
#'
#' @return A usually sparse covariance matrix of random effects of type dsCMatrix
#'
#' @details
#' If psi is not supplied, the estimated covariance matrix is returned. If psi is
#' supplied, the covariance matrix is calculated using those parameter values instead.
#' psi should be NULL (default) or a numeric vector of length getME(lmerfit, "m").
#' In the latter case the elements should be ordered as those in the vcov column of
#' as.data.frame(VarCorr(lmerfit), order = "lower.tri"). In particular, the last
#' element is the error variance.
#'
#'
#' @export
get_Psi <- function(lmerfit, psi = NULL){
  if(is.null(psi)){
    # Extract variances and covariances of random effects ordered as in the covmat
    # lower.tri may seem wrong since we are creating upper triangular Psi,
    # but appears to conform with the indexing in getME(, "Lind")
    psi <- as.data.frame(VarCorr(lmerfit), order = "lower.tri")$vcov
  }

  # Lambdat has the right structure, but not the same entries as Psi
  Psi_half <- lme4::getME(lmerfit, "Lambdat")

  # Separate error variance and covariance matrix for random effects
  psi0 <- psi[length(psi)]
  psi <- psi[-length(psi)]

  # Fill in the upper triangular part of Psi with the extracted elements
  param_idx <- lme4::getME(lmerfit, "Lind")
  Psi_half@x <- psi[param_idx]

  # Return symmetric matrix
  Matrix::forceSymmetric(Psi_half, uplo = "U")
}

# Make list of H matrices
get_H_list <- function(lmerfit)
{
  # Psi, and hence H, has the same structure as Lambda
  H <- lme4::getME(lmerfit, "Lambdat")
  param_idx <- lme4::getME(lmerfit, "Lind")

  # Replace values by parameter index
  H@x <- param_idx
  # Return list of H matrices
  lapply(seq_len(getME(lmerfit, "m")),
         function(i){
           M <- H
           M@x <- as.numeric(i == M@x)
           Matrix::forceSymmetric(Matrix::drop0(M), uplo = "U")
         })
}

# function which takes input of an lme4 model object and returns the
# log likelihood, score, and fisher information matrix with and without use of
# the restricted likelihood function
lmm_stuff <- function(fit, psiNull_error, psiNull_re, loglik = TRUE,
                      score = TRUE, finf = TRUE) {
  # Obtaining the matrix H
  Hlist <- get_Hlist(fit)

  # Forming the matrix Psi from the vector of null values
  Psi <- Psi_from_Hlist(Hlist, psiNull_re)
  psi0 <- psiNull_error
  # Dividing Psi by psi0 for use in the limestest functions
  Psi0 <- Psi/psi0

  # extract the rest of the design from the model
  X <- getME(fit, "X")
  y <- getME(fit, "y")
  Z <- getME(fit, "Z")

  # cross products
  XtX <- crossprod(X)
  ZtZ <- crossprod(Z)
  XtZ <- crossprod(X, Z)
  XtY <- crossprod(X, y)
  YtZ <- crossprod(y, Z)

  # residual log likelihood, score, and finf
  stuff_REML <- limestest::res_ll(XtX = XtX,
                                  XtY = XtY,
                                  XtZ = XtZ,
                                  ZtZ = ZtZ,
                                  YtZ = YtZ,
                                  Y = y,
                                  X = X,
                                  Z = Z,
                                  H = H,
                                  Psi0 = Psi0,
                                  psi0 = psi0,
                                  loglik,
                                  score,
                                  finf)

  # calculate residuals
  e <- y - X %*% stuff_REML$beta
  # cross product w/ residuals
  Zte <- crossprod(Z, e)

  # log likelihood, score, and finf of psi at the null values
  stuff <- limestest::loglik_psi(Z = Z,
                                 ZtZXe = cbind(ZtZ, t(XtZ), Zte),
                                 e = e,
                                 H = H,
                                 Psi0 = Psi0,
                                 psi0 = psi0,
                                 loglik,
                                 score,
                                 finf)
  return(list(stuff, stuff_REML))
}

# function which takes input of an lme4 model object, and null values of:
# psi0, the variance of the elements of the error vector and
# Psi, the covariance matrix of the random effects
# The function returns values of the score statistics with and without use of
# the restricted likelihood function.
#' @export
lmm_scorestat <- function(fit, psiNull_error, psiNull_re) {

  # get fisher info and score
  stuff2 <- lmm_stuff(fit, psiNull_error, psiNull_re)
  stuff <- stuff2[[1]]
  stuff_REML <- stuff2[[2]]

  # score statistics
  test_stat <- as.vector(crossprod(stuff$score, solve(stuff$finf, stuff$score)))
  test_stat_REML <- as.vector(crossprod(stuff_REML$score,
                                        solve(stuff_REML$finf, stuff_REML$score)))

  # return
  return(c(test_stat, test_stat_REML))
}

# Function for calculating the score statistic for a single parameter at a time
# with all other parameters fixed at their maximum likelihood estimate
# psi is the value of our single parameter at which to evaluate the score
# j is the index of the parameter in the vector of covariance parameters
scoreStatOneParam <- function(fit, psi, j) {
  # maximum likelihood estimate of the error variance
  psihat_error <- sigma(fit)^2

  # maximum likelihood estimate "psihat"

  ## OLD ##
  #Lambda <- getME(fit, "Lambda")
  #Psi <- Lambda %*% t(Lambda)
  #unPsiVals <- unique(Psi@x)
  #m <- getME(fit, "m")
  #while (length(unPsiVals) < m) {
  #  unPsiVals <- c(unPsiVals, 0)
  #}
  #psihat <- unPsiVals
  ##
  psihat <- as.data.frame(VarCorr(fit))[,"vcov"]
  # this was missing earlier: Lambda is scaled by the error variance so we need to do this
  psihat <- psihat*psihat_error

  # replacing the jth element with our value
  psihat[j] <- psi

  # evaluating the score
  allstuff <- lmm_stuff(fit, psihat_error, psihat)
  # choose REML or not right here: index 1 not REML, index 2 is REML
  stuff <- allstuff[[1]]
  # obtain the score for the given index
  # (first index of this vector is for the error variance)
  score <- stuff$score[j+1]

  # calculate efficient information (ei)
  I <- stuff[[3]][-1, -1]
  ei <- I[j,j] - I[j, -j] %*% solve(I[-j, -j], I[-j, j])
  ei <- ei[1,1]

  # divide by the sqrt of the efficient information to get score test stat
  scorestat <- score/sqrt(ei)
  return(scorestat)
}
