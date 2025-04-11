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
    # lower.tri despite creating upper triangular Psi since it appears to
    # with the indexing in getME(, "Lind")
    psi <- as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$vcov
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
  lapply(seq_len(lme4::getME(lmerfit, "m")),
         function(i){
           M <- H
           M@x <- as.numeric(i == M@x)
           Matrix::forceSymmetric(Matrix::drop0(M), uplo = "U")
         })
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
