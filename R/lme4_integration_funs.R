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
get_Psi <- function(lmerfit, psimr = NULL){
  if(is.null(psimr)){
    # Extract variances and covariances of random effects ordered as in the covmat
    # lower.tri despite creating upper triangular Psi since it appears to
    # with the indexing in getME(, "Lind")
    rm1 <- lme4::getME(lmerfit, "m")
    psimr <- as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$vcov[1:rm1]
  }

  # Lambdat has the right structure, but not the same entries as Psi
  Psi_half <- lme4::getME(lmerfit, "Lambdat")

  # Fill in the upper triangular part of Psi with the extracted elements
  param_idx <- lme4::getME(lmerfit, "Lind")
  Psi_half@x <- psimr[param_idx]

  # Return symmetric matrix
  Matrix::forceSymmetric(Psi_half, uplo = "U")
}

# Make list of H matrices
get_Hlist <- function(lmerfit)
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

get_precomp <- function(lmerfit){
  # 0 indicates ML
  REML <- lme4::getME(lmerfit, "REML") != 0

  Y <- lme4::getME("y")
  X <- lme4::getME("X")
  Z <- lme4::getME("Z")

  if(!REML){
    b <- lme4::getME("beta")
    Y <- Y - X %*% b
  }

  list(XtX = crossprod(X),
       XtY = crossprod(X, Y),
       XtZ = crossprod(X, Z),
       ZtZ = crossprod(Z),
       YtZ = crossprod(Y, Z),
       Y = Y,
       X = X,
       Z = Z,
       Hlist = limestest:::get_Hlist(fit))
}

get_psi_hat(lmerfit)
{
  as.data.frame(VarCorr(lmerfit), order = "lower.tri")$vcov
}

lmer_score_test <- function(lmerfit,
                            psi_null = NULL,
                            test_idx = NULL,
                            efficient = TRUE,
                            expeted = TRUE,
                            profile = TRUE,
                            joint = TRUE)
{
  r_i <- lme4::getME(lmerfit, "m_i")
  r <- sum(r_i) + 1
  p <- ncol(precomp$X)
  n <- nrow(precomp$X)
  psi_hat <- get_psi_hat(lmerfit)
  precomp <- get_precomp(lmerfit)
  psi_start <- lme4::getME(lmerfit, "Tp")[-1]
  REML <- lme4::getME(lmerfit, "REML") != 0

  if(is.null(test_idx)){
    # Test all RE parameters equal to zero either separately or jointly
    if(joint){
      psi_null <- psi_hat
      psi_null[-r] <- 0
      psi_null[r] <- sigma(lm(precomp$Y ~ 0 + precomp$X))^2
      if(!REML){
        psi_null[r] <- psi_null[r] * (n - p) / n
      }
      teststat <- c(score_stat(psi = psi_null, test_idx = 1:(r - 1), precomp = precomp,
                 REML = REML, expected = expected,
                 efficient = efficient, signed = FALSE))
      out <- matrix(c(teststat, pchisq(teststat, df = r - 1, lower = F)), nrow = 1)
      colnames(out) <- c("stat", "pval", "df")
      rownames(out) <- "joint"
    } else{ # Separate tests
      out <- matrix(NA, nrow = r - 1, ncol = 3)
      param_idx <- 1
      psi_null <- psi_hat
      for(ii in seq_along(r_i)){ # Loop over terms
        term_idxs <- psi_start[ii]:(psi_start[ii] + r_i[ii] - 1)
        for(jj in seq_len(r_i[ii])){ # Loop over parameters within terms
          psi_null <- psi_hat
          psi_null[term_idxs] <- 0
          psi_null <- partial_min(opt_idx = seq_len(r)[-param_idx],
                                  precomp = precomp, psi_start = psi_null, REML = REML,
                                  expected = expected)$psi_hat
          out[param_idx, 1] <- c(score_stat(psi = psi_null, test_idx = param_idx,
                                            precomp = precomp,
                                            REML = REML, expected = expected,
                                            efficient = efficient, signed = FALSE))
          out[param_idx, 2] <- pchisq(out[param_idx, 1], df = 1, lower = F)
          out[param_idx, 3] <- 1
          # Increase parameter index
          param_idx <- param_idx + 1
        }
      }
    }
  } else if(!is.null(psi_null)){
    # Test the hypothesis
  } else{
    stop("Unable to test because null hypothesis parameter (psi_null) is missing")
  }

}



