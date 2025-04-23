#' Get the covariance matrix of random effects
#'
#' Returns the covariance matrix of the random effects in a linear mixed model,
#' either as estimated by lme4::lmer or evaluated at a particular parameter value
#' supplied as an argument.
#'
#' @param lmerfit An lmerMod object from fitting a linear mixed model using lme4::lmer
#' @param psimr Optional vector with covariance parameters, not including the error variance;
#   psi minus its r:th element (see details).
#'
#' @return A usually sparse covariance matrix of random effects of type dsCMatrix
#'
#' @details
#' If psimr is not supplied, the estimated covariance matrix is returned. If psimr is
#' supplied, the covariance matrix is calculated using those parameter values instead.
#' psimr should be NULL (default) or a numeric vector of length getME(lmerfit, "m").
#' In the latter case the elements should be ordered as those in the vcov column of
#' as.data.frame(VarCorr(lmerfit), order = "lower.tri"). The last
#' element in that column is is the error variance, which should be omitted.
#'
#'
#' @export
get_Psi_lmer <- function(lmerfit, psi_mr = NULL){
  if(is.null(psi_mr)){
    # Extract variances and covariances of random effects ordered as in the covmat
    # lower.tri despite creating upper triangular Psi since it appears to
    # with the indexing in getME(, "Lind")
    rm1 <- lme4::getME(lmerfit, "m")
    psi_mr <- as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$vcov[1:rm1]
  }

  # Lambdat has the right structure, but not the same entries as Psi
  Psi_half <- lme4::getME(lmerfit, "Lambdat")

  # Fill in the upper triangular part of Psi with the extracted elements
  param_idx <- lme4::getME(lmerfit, "Lind")
  Psi_half@x <- psi_mr[param_idx]

  # Return symmetric matrix
  Matrix::forceSymmetric(Psi_half, uplo = "U")
}

# Make list of H matrices
get_Hlist_lmer <- function(lmerfit)
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

get_precomp_lmer <- function(lmerfit, REML = NULL){

  if(is.null(REML)){
    # 0 indicates ML
    REML <- lme4::getME(lmerfit, "REML") != 0
  }

  X <- lme4::getME(lmerfit, "X")
  Z <- lme4::getME(lmerfit, "Z")

  if(REML){
    Y <- lme4::getME(lmerfit, "y")
    out <- list(XtY = as.vector(crossprod(X, Y)),
                ZtY = as.vector(crossprod(Z, Y)),
                XtX = as.matrix(crossprod(X)),
                XtZ = as.matrix(crossprod(X, Z)),
                ZtZ = methods::as(crossprod(Z), "generalMatrix"))
  } else{
    # Warning: These are not equal to residuals(lmerfit)
    e <- Y - X %*% lme4::getME(lmerfit, "beta")
    out <- list(e = as.vector(e),
                Zte = as.vector(crossprod(Z, e)),
                XtZ = as.matrix(crossprod(X, Z)),
                ZtZ = methods::as(crossprod(Z), "generalMatrix"))
  }
  # Return
  out
}

get_psi_hat_lmer <- function(lmerfit)
{
  as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$vcov
}

score_test_lmer <- function(lmerfit,
                            psi_null = NULL,
                            test_idx = NULL,
                            efficient = TRUE,
                            expected = TRUE,
                            profile = TRUE,
                            joint = TRUE)
{
  r_i <- lme4::getME(lmerfit, "m_i")
  r <- sum(r_i) + 1
  precomp <- get_precomp_lmer(lmerfit)

  Y <- lme4::getME(lmerfit, " y")
  X <- lme4::getME(lmerfit, "X")
  Z <- lme4::getME(lmerfit, "Z")
  Hlist <- get_Hlist_lmer(lmerfit)

  p <- ncol(X)
  n <- nrow(X)

  psi_hat <- get_psi_hat_lmer(lmerfit)

  last_idx <- lme4::getME(lmerfit, "Tp")

  REML <- lme4::getME(lmerfit, "REML") != 0

  if(is.null(test_idx)){
    # Test all RE parameters equal to zero either separately or jointly
    if(joint){
      psi_null <- psi_hat
      psi_null[-r] <- 0
      psi_null[r] <- stats::sigma(stats::lm(Y ~ 0 + X))^2
      if(!REML){
        psi_null[r] <- psi_null[r] * (n - p) / n
      }

      test_stat <- score_stat(psi = psi_null,
                             test_idx = 1:(r - 1),
                             b = NULL,
                             Y = Y,
                             X = X,
                             Z = Z,
                             Hlist = Hlist,
                             REML = REML,
                             expected = expected,
                             efficient = efficient,
                             signed = FALSE,
                             precomp = precomp)

      out <- matrix(c(test_stat, stats::pchisq(test_stat, df = r - 1, lower = F),
                      r - 1), nrow = 1, ncol = 3)

      colnames(out) <- c("stat", "pval", "df")
      rownames(out) <- "joint"
    } else{ # Separate tests
      out <- matrix(NA, nrow = r - 1, ncol = 3)
      param_idx <- 1
      psi_null <- psi_hat
      for(ii in seq_along(r_i)){ # Loop over terms
        term_idxs <- (last_idx[ii + 1] - r_i[ii] + 1):(last_idx[ii + 1])
        for(jj in seq_len(r_i[ii])){ # Loop over parameters within terms
          psi_null <- psi_hat
          psi_null[term_idxs] <- 0
          psi_null <- partial_min_psi(psi_start = psi_null,
                                      opt_idx = seq_len(r)[-param_idx],
                                      b = NULL,
                                      Y = Y,
                                      X = X,
                                      Z = Z,
                                      Hlist = Hlist,
                                      precomp = precomp,
                                      REML = REML,
                                      expected = expected)$psi_hat
          out[param_idx, 1] <- c(score_stat(psi = psi_null,
                                            test_idx = param_idx,
                                            b = NULL,
                                            Y = Y,
                                            X = X,
                                            Z = Z,
                                            Hlist = Hlist,
                                            REML = REML,
                                            expected = expected,
                                            efficient = efficient,
                                            signed = FALSE,
                                            precomp = precomp))
          out[param_idx, 2] <- stats::pchisq(out[param_idx, 1], df = 1, lower = F)
          out[param_idx, 3] <- 1
          # Increase parameter index
          param_idx <- param_idx + 1
        }
      }
      rownames(out) <- paste0("psi_", seq_len(r - 1))
      colnames(out) <- c("stat", "pval", "df")
    }
  } else if(!is.null(psi_null)){
    # Test the hypothesis
  } else{
    stop("Unable to test because null hypothesis parameter (psi_null) is missing")
  }
  out
}



