#' Get the covariance matrix of random effects
#'
#' Returns the covariance matrix of the random effects vector in a linear mixed model,
#' either as estimated by lme4::lmer or evaluated at particular parameter values
#' supplied as an argument.
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using `lme4::lmer`.
#' @param psi_mr Optional numeric vector with covariance parameters, not including the 
#'   error variance; i.e., `psi` minus its r:th element, where r is the total number
#'   of covariance parameters (see details).
#'
#' @return A sparse symmetric covariance matrix (dsCMatrix) of the random effects vector.
#'
#' @details
#' If `psi_mr` is not supplied, the estimated covariance matrix from the fitted model 
#' is returned. If `psi_mr` is supplied, the covariance matrix is constructed using 
#' those parameter values instead.
#' 
#' When provided, `psi_mr` should be a numeric vector of length `getME(lmerfit, "m")`,
#' which equals r - 1 where r is the total number of covariance parameters (including
#' error variance). The elements should be ordered as in the `vcov` column of
#' `as.data.frame(VarCorr(lmerfit), order = "lower.tri")`, with the last element 
#' (error variance) omitted.
#'
#' @export
get_Psi_lmer <- function(lmerfit, psi_mr = NULL){
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }
  
  rm1 <- lme4::getME(lmerfit, "m")
  
  if(is.null(psi_mr)){
    # Extract variances and covariances of random effects ordered as in the covmat
    # lower.tri despite creating upper triangular Psi since it appears consistent
    # with the indexing in getME(, "Lind")
    psi_mr <- as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$vcov[1:rm1]
  } else {
    # Validate psi_mr
    if (!is.numeric(psi_mr)) {
      stop("psi_mr must be a numeric vector")
    }
    if (length(psi_mr) != rm1) {
      stop("psi_mr must have length ", rm1, " (getME(lmerfit, 'm'))")
    }
  }

  # Lambdat has the right structure, but not the same entries as Psi
  Psi_half <- lme4::getME(lmerfit, "Lambdat")

  # Fill in the upper triangular part of Psi with the extracted elements
  param_idx <- lme4::getME(lmerfit, "Lind")
  Psi_half@x <- psi_mr[param_idx]

  # Return symmetric matrix
  Matrix::forceSymmetric(Psi_half, uplo = "U")
}

#' Get structure matrices for covariance parameterization
#'
#' Extracts the list of structure matrices (H matrices) from an lme4 fit that
#' determine how covariance parameters map to the covariance matrix structure.
#' These matrices are used in likelihood computations where the covariance matrix
#' is expressed as a linear combination: Psi = sum(psi[i] * H[[i]]).
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using
#'   `lme4::lmer`.
#'
#' @return A list of sparse symmetric matrices (dsCMatrix), one for each
#'   covariance parameter (excluding error variance). The length of the list
#'   equals `getME(lmerfit, "m")`, which is r - 1 where r is the total number
#'   of covariance parameters including error variance.
#'
#' @details
#' Each matrix in the returned list is an indicator matrix showing which elements
#' of the random effects covariance matrix are associated with each parameter.
#' The i-th matrix has 1s in positions determined by the i-th covariance parameter
#' and 0s elsewhere.
#'
#' @keywords internal
get_Hlist_lmer <- function(lmerfit)
{
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }
  
  # Psi, and hence H, has the same structure as Lambdat
  H <- lme4::getME(lmerfit, "Lambdat")
  param_idx <- lme4::getME(lmerfit, "Lind")

  # Replace values by parameter index
  H@x <- param_idx
  
  # Create list of indicator matrices, one for each covariance parameter
  # Each matrix has 1s where that parameter appears, 0s elsewhere
  lapply(seq_len(lme4::getME(lmerfit, "m")),
         function(i){
           M <- H
           M@x <- as.numeric(i == M@x)
           Matrix::forceSymmetric(Matrix::drop0(M), uplo = "U")
         })
}

#' Get precomputed quantities from lme4 fit
#'
#' Extracts and computes quantities from an lme4 fit that can be reused in
#' likelihood calculations to avoid redundant computations. The quantities
#' computed depend on whether REML or ML estimation is used.
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using
#'   `lme4::lmer`.
#' @param REML Logical indicating whether to compute quantities for REML
#'   (\code{TRUE}) or ML (\code{FALSE}). If \code{NULL} (default), uses the
#'   estimation method from the fitted model.
#'
#' @return A list containing precomputed cross-products and other quantities:
#'   \itemize{
#'     \item For REML: \code{XtY}, \code{ZtY}, \code{XtX}, \code{XtZ}, \code{ZtZ}
#'     \item For ML: \code{e} (residuals), \code{Zte}, \code{XtX}, \code{XtZ}, \code{ZtZ}
#'   }
#'
#' @details
#' The residuals \code{e} in the ML case are computed as Y - X * beta where beta
#' comes from the fitted model. These may differ slightly from \code{residuals(lmerfit)}
#' due to how lme4 computes residuals internally.
#'
#' @keywords internal
get_precomp_lmer <- function(lmerfit, REML = NULL){
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }
  
  if(is.null(REML)){
    # 0 indicates ML, non-zero indicates REML
    REML <- lme4::getME(lmerfit, "REML") != 0
  } else {
    if (!is.logical(REML) || length(REML) != 1) {
      stop("REML must be a single logical value")
    }
  }

  X <- lme4::getME(lmerfit, "X")
  Z <- lme4::getME(lmerfit, "Z")
  Y <- lme4::getME(lmerfit, "y")
  
  if(REML){
    out <- list(XtY = as.vector(crossprod(X, Y)),
                ZtY = as.vector(crossprod(Z, Y)),
                XtX = as.matrix(crossprod(X)),
                XtZ = as.matrix(crossprod(X, Z)),
                ZtZ = methods::as(crossprod(Z), "generalMatrix"))
  } else{
    # Compute residuals from fitted fixed effects
    # Note: may differ from residuals(lmerfit) due to lme4's internal computation
    e <- Y - X %*% lme4::getME(lmerfit, "beta")
    out <- list(e = as.vector(e),
                Zte = as.vector(crossprod(Z, e)),
                XtX = as.matrix(crossprod(X)),
                XtZ = as.matrix(crossprod(X, Z)),
                ZtZ = methods::as(crossprod(Z), "generalMatrix"))
  }
  out
}

#' Extract estimated covariance parameters from lme4 fit
#'
#' Extracts all estimated variance and covariance parameters from a fitted
#' linear mixed model, including the error variance.
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using
#'   `lme4::lmer`.
#'
#' @return A numeric vector containing all estimated covariance parameters,
#'   ordered as in the `vcov` column of 
#'   `as.data.frame(VarCorr(lmerfit), order = "lower.tri")`. The last element
#'   is the error variance. The vector has length r, where r is the total
#'   number of covariance parameters.
#'
#' @details
#' This function extracts the full vector of covariance parameter estimates
#' (often denoted psi or theta in the package), including:
#' \itemize{
#'   \item Variances and covariances of random effects
#'   \item Error variance (last element)
#' }
#' 
#' The ordering follows lme4's internal parameterization with "lower.tri" ordering.
#'
#' @keywords internal
get_psi_hat_lmer <- function(lmerfit)
{
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }
  
  # Extract variance components
  vcov_vec <- as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$vcov
  
  if (!is.numeric(vcov_vec) || length(vcov_vec) == 0) {
    stop("Failed to extract variance components from lmerfit")
  }
  
  vcov_vec
}

score_test_lmer <- function(lmerfit,
                            psi_null = NULL,
                            test_idx = NULL,
                            efficient = TRUE,
                            expected = TRUE,
                            profile = TRUE)
{

  r_i <- lme4::getME(lmerfit, "m_i")
  r <- sum(r_i) + 1
  # Default to testing zero random effect variances and unit error variance
  if(is.null(psi_null)){
    psi_null <- c(rep(0, r - 1), 1)
  }

  stopifnot(is.numeric(psi_null) && length(psi_null) > 0)

  # Test all random effects by default
  if(is.null(test_idx)){
    test_idx <- seq_len(r - 1)
  } else{
    test_idx <- sort(unique(test_idx))
  }

  stopifnot(is.numeric(test_idx) && all(test_idx == floor(test_idx)) && length(test_idx) > 0)
  k <- length(test_idx)


  if(psi_null[r] == 0){
    stop("Vanishing error variance is not permitted")
  }

  precomp <- get_precomp_lmer(lmerfit)
  Y <- lme4::getME(lmerfit, "y")
  X <- lme4::getME(lmerfit, "X")
  Z <- lme4::getME(lmerfit, "Z")
  Hlist <- get_Hlist_lmer(lmerfit)
  n <- length(Y)
  p <- ncol(X)
  REML <- lme4::getME(lmerfit, "REML") != 0

  # Profile
  if(profile && all(psi_null[-r] == 0) && all(seq_len(r - 1) %in% test_idx)){
    # Testing null of Psi = 0, i.e., no random effects is a
    # special case where partial maximizer has closed form solution from lm
    psi_null[r] <- stats::sigma(stats::lm(Y ~ 0 + X))^2
    if(!REML){ # Degree of freedom correction only for REML
      psi_null[r] <- psi_null[r] * (n - p) / n
    }
  } else if(profile){
    psi_null <- partial_min_psi(psi_start = psi_null,
                                opt_idx = seq_len(r)[-test_idx],
                                b = NULL,
                                Y = Y,
                                X = X,
                                Z = Z,
                                Hlist = Hlist,
                                precomp = precomp,
                                REML = REML,
                                expected = expected)$psi_hat
  }

  test_stat <- score_stat(psi = psi_null,
                          test_idx = test_idx,
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
  # Return
  c("stat" = test_stat,
    "p_val" = stats::pchisq(test_stat, df = k, lower = FALSE),
    "df" = k)
}

test_all_lmer <- function(lmerfit,
                           psi_null = NULL,
                           test_idx = NULL,
                           efficient = TRUE,
                           expected = TRUE)
{
  r_i <- lme4::getME(lmerfit, "m_i")
  r <- sum(r_i) + 1

  # Default to testing zero random effect variances and unit error variance
  if(is.null(psi_null)){
    psi_null <- c(rep(0, r - 1), 1)
  }

  # Test all random effects by default
  if(is.null(test_idx)){
    test_idx <- seq_len(r - 1)
  } else{
    test_idx <- sort(unique(test_idx))
  }

  if(psi_null[r] == 0){
    stop("Vanishing error variance is not permitted")
  }

  precomp <- get_precomp_lmer(lmerfit)
  Y <- lme4::getME(lmerfit, "y")
  X <- lme4::getME(lmerfit, "X")
  Z <- lme4::getME(lmerfit, "Z")
  Hlist <- get_Hlist_lmer(lmerfit)
  psi_hat <- get_psi_hat_lmer(lmerfit)
  REML <- lme4::getME(lmerfit, "REML") != 0
  k <- length(test_idx)

  # Is parameter a variance or covaraince?
  is_var_param <- is.na(as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$var2)
  stopifnot(all(psi_null[is_var_param] >= 0))

  # Prepare to loop over parameters

  # Used to determine which term a parameter belongs to
  last_idx <- lme4::getME(lmerfit, "Tp")

  out <- matrix(NA, nrow = k, ncol = 3) # Store test results
  param_idx <- 1 # Counts which parameter test is for
  test_num <- 1  # Counts the tests
  for(ii in seq_along(r_i)){ # Loop over terms
    # Index for parameters corresponding to term
    term_idxs <- (last_idx[ii + 1] - r_i[ii] + 1):(last_idx[ii + 1])

    # Covariance matrix dimension for ith term
    dim_i <- as.integer(0.5 * (-1 + sqrt(1 + 8 * r_i[ii])))

    for(jj in seq_len(r_i[ii])){ # Loop over parameters within terms
      if(param_idx %in% test_idx){
        # Create starting point for profile optimization that is guaranteed
        # to be valid
        psi_start <- c(rep(0, r - 1), 1) # Start model with no random effects
        # psi_start[term_idxs][!is_var_param[term_idxs]] <- 0 # Set off-diagonal to zero
        psi_start[param_idx] <- psi_null[param_idx] # Set null hypothesis value
        if(!is_var_param[param_idx]){
         # If the null hypothesis is for a covariance, set corresponding variances
         # to ensure positive semi-definite starting value for Psi. This works
         # because other off-diagonal elements were set to zero
         row_col <- get_row_col_ltri(jj, n = dim_i)
         var_idx1 <- get_idx_ltri(row = row_col[2], col = row_col[2], n = dim_i)
         var_idx2 <- get_idx_ltri(row = row_col[1], col = row_col[1], n = dim_i)
         psi_start[term_idxs][c(var_idx1, var_idx2)] <-
           pmax(psi_start[term_idxs][c(var_idx1, var_idx2)], abs(psi_start[param_idx]))
        }
        out[test_num, ] <- score_test_lmer(lmerfit = lmerfit,
                                           psi_null = psi_start,
                                           test_idx = param_idx,
                                           efficient = efficient,
                                           expected = expected,
                                           profile = TRUE)
        test_num <- test_num + 1

      }
      param_idx <- param_idx + 1
    }
  }
  out
}



