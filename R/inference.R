Psi_from_Hlist <- function(psi_mr, Hlist)
{
  rm1 <- length(Hlist)
  if(rm1 >= 1){
    q <- ncol(Hlist[[1]])
    Psi <- Matrix::sparseMatrix(i = seq_len(0), j = seq_len(0), x = 0, dims = c(q, q))
    for(ii in seq_len(rm1)){
      Psi <- Psi + Hlist[[ii]] * psi_mr[ii]
    }
  } else{
    Psi <- 0
  }
  Psi
}

maximize_loglik <- function(start_val, opt_idx, Y, X, Z, Hlist, expected = TRUE,
                             REML = TRUE, precomp = NULL, ...)
{
  r <- length(Hlist) + 1
  p <- ncol(X)
  assertthat::assert_that((REML && length(start_val) == r) ||
                            (!REML && length(start_val) == p + r),
                          msg = "start_val should be length p + r for ML and r for REML")
  if(is.null(precomp)) {
    precomp <- list("XtX" = crossprod(X),
                    "XtZ" = as.matrix(crossprod(X, Z)),
                    "ZtZ" = methods::as(crossprod(Z), "generalMatrix"))
    if(REML) {
      precomp$XtY <- as.vector(crossprod(X, Y))
      precomp$ZtY <- as.vector(crossprod(Z, Y))
    }
  }

  #############################################################################
  # Define the objective function to be minimized
  #############################################################################
  if(!REML){
    obj_fun <- function(x) {
      theta <- start_val
      theta[opt_idx] <- x
      psi <- theta[(p + 1):(r + p)]
      H <- do.call(cbind, Hlist)
      Psi <- Psi_from_H_cpp(psi_mr = psi[-r], H = H)
      b <- if (p == 0) NULL else theta[1:p]
      e <- if (p == 0) Y else Y - X %*% b
      ll_things <- loglik(Psi_r = Psi / psi[r],
                          psi_r = psi[r],
                          H = H,
                          e = e,
                          X = X,
                          Z = Z,
                          XtX = precomp$XtX,
                          XtZ = precomp$XtZ,
                          ZtZ = precomp$ZtZ,
                          get_val = TRUE,
                          get_score = TRUE,
                          get_inf = TRUE,
                          expected = expected)
      list("value" = -ll_things$value, "gradient" = -ll_things$score[opt_idx],
           "hessian" = as.matrix(ll_things$inf_mat[opt_idx, opt_idx]))
    }
  } else{
    obj_fun <- function(x) {
      psi <- start_val
      psi[opt_idx] <- x
      H <- do.call(cbind, Hlist)
      Psi <- Psi_from_H_cpp(psi_mr = psi[-r], H = H)
      ll_things <- loglik_res(Psi_r = Psi / psi[r],
                              psi_r = psi[r],
                              H = H,
                              Y = Y,
                              X = X,
                              Z = Z,
                              XtX = precomp$XtX,
                              XtZ = precomp$XtZ,
                              ZtZ = precomp$ZtZ,
                              XtY = precomp$XtY,
                              ZtY = precomp$ZtY,
                              get_val = TRUE,
                              get_score = TRUE,
                              get_inf = TRUE)
      list("value" = -ll_things$value, "gradient" = -ll_things$score[opt_idx],
           "hessian" = as.matrix(ll_things$inf_mat[opt_idx, opt_idx]))
    }
  }

  #############################################################################
  # Do minimization
  #############################################################################
  fit <- trust::trust(objfun = obj_fun, parinit = start_val[opt_idx], rinit  = 1,
                      rmax = 100, ...)
  # Return results
  start_val[opt_idx] <- fit$argument
  names(start_val) <- if(REML) paste0("psi", 1:r) else c(paste0("b", 1:p), paste0("psi", 1:r))
  list("arg" = start_val, "value" = -fit$value, "conv" = fit$converged,
       "iter" = fit$iterations)
}


score_stat <- function(theta, test_idx, Y, X, Z, Hlist, REML = TRUE,
                       expected = TRUE, efficient = TRUE, signed = FALSE,
                       precomp = NULL)
{
  if(!expected && REML){
    warning("Observed information not available for restricted likelihood; using
            expected.")
  }
  test_idx <- unique(test_idx)
  p <- ncol(X)
  r <- length(Hlist) + 1

  psi <- if(REML || p < 1) theta else theta[-seq_len(p)]
  b <- if(!REML && p >= 1) theta[seq_len(p)] else NULL

  ll_things <- loglikelihood(psi = psi,
                             b = b,
                             Y = Y,
                             X = X,
                             Z = Z,
                             Hlist = Hlist,
                             REML = REML,
                             get_val = FALSE,
                             get_score = TRUE,
                             get_inf = TRUE,
                             get_beta = (!REML && p>=1),
                             expected = expected,
                             precomp = precomp)

  inf_mat <- ll_things$inf_mat[test_idx, test_idx, drop = F]

  if(efficient && (r > 1)){
    inf_mat <- inf_mat - crossprod(ll_things$inf_mat[-test_idx, test_idx, drop = F],
                                   solve(ll_things$inf_mat[-test_idx, -test_idx, drop = F],
                                         ll_things$inf_mat[-test_idx, test_idx, drop = F]))
  }

  if(signed){
    ed <- eigen(inf_mat, symmetric = TRUE)
    inf_mat <- ed$vectors %*% (sqrt(ed$values) * t(ed$vectors))

    # Inefficient since eigen decomposition available, but may be more stable
    test_stat <- solve(inf_mat, ll_things$score[test_idx])
  } else{
    test_stat <- crossprod(ll_things$score[test_idx], solve(inf_mat,
                                                             ll_things$score[test_idx]))
  }
  as.vector(test_stat)
}

#' Score test statistic with nuisance parameters
#'
#' Computes the score test statistic when nuisance parameters handled by
#' unconstrained optimization.
#'
#' @param theta_null Vector of parameter values under the null hypothesis to
#'   search around. Has to be a valid parameter in the null hypothesis set. If
#'   \code{REML = FALSE}, this should be of length \eqn{p + r} with
#'   \code{theta_null = c(beta, psi)}. If \code{REML = TRUE}, this should be of
#'   length \eqn{r} with \code{theta_null = psi}. 
#' @param test_idx Integer specifying which element of \code{theta_null}
#'   to test.
#' @param max_radius Numeric value or vector of length 2 specifying the radius
#'   around \code{theta_null[test_idx]} to search. If length 1, search is
#'   symmetric. If length 2, \code{max_radius[1]} specifies the lower radius and
#'   \code{max_radius[2]} the upper radius. If 0 (default), only theta_null is
#' considered @param num_points Integer specifying the number of points to
#' evaluate in the range
#'   defined by \code{max_radius}. Default is 100.
#' @param Y Vector of length \eqn{n} of responses.
#' @param X Matrix of size \eqn{n \times p} of fixed effect predictors.
#' @param Z Sparse \eqn{n \times q} random effect design matrix.
#' @param Hlist List of matrices determining how \eqn{\psi} is mapped to
#' \eqn{\Psi}
#'   (see \code{?loglikelihood}).
#' @param REML Logical. If \code{TRUE}, use restricted likelihood; otherwise use
#'   regular likelihood. Default is \code{TRUE}.
#' @param expected Logical. If \code{TRUE}, use expected information matrix;
#' otherwise
#'   use observed. Default is \code{TRUE}.
#' @param efficient Logical. If \code{TRUE}, use the efficient information.
#' Default is \code{TRUE}.
#' @param signed Logical. If \code{TRUE}, return signed score test statistic;
#' otherwise
#'   return squared (chi-squared type) statistic. Default is \code{TRUE}.
#' @param precomp Optional list of pre-computed quantities (see
#' \code{?loglikelihood}).
#'   If \code{NULL}, these will be computed internally.
#'
#' @return Named numeric vector of score test statistics. Names correspond to the
#'   values of the test parameter at which the statistic was evaluated.
#'
#' @details This function is useful for constructing confidence intervals or
#' performing hypothesis tests that account for uncertainty in nuisance
#' parameters. For each value in the specified range around
#' \code{theta_null[test_idx]}, the function: \enumerate{ \item Fixes the test
#' parameter at that value \item Optimizes the remaining (nuisance) parameters
#'   \item Computes the score test statistic at the resulting parameter values
#' }
#'
#' The search proceeds in two directions from \code{theta_null[test_idx]}: first
#' decreasing values, then increasing values. This allows for efficient
#' warm-starting of the optimization at each step.
score_nuisance <- function(theta_null, test_idx, max_radius = 0, num_points = 1e2,
                           Y, X, Z, Hlist, REML = TRUE, expected = TRUE,
                           efficient = TRUE, signed = TRUE, precomp = NULL) {
  # Argument checking
  assertthat::assert_that(is.numeric(theta_null), length(theta_null) > 0,
                          msg = "theta_null should be a numeric vector of positive length")

  assertthat::assert_that(is.numeric(test_idx), length(test_idx) == 1,
                          test_idx >= 1, test_idx <= length(theta_null),
                          msg = "test_idx should be a single integer between 1 and length(theta_null)")

  assertthat::assert_that(is.numeric(max_radius), length(max_radius) %in% c(1, 2),
                          all(max_radius >= 0),
                          msg = "max_radius should be a numeric value or vector of length 1 or 2 with non-negative values")

  assertthat::assert_that(is.numeric(num_points), length(num_points) == 1,
                          num_points >= 1,
                          msg = "num_points should be a positive integer")

  assertthat::assert_that(is.vector(Y, mode = "numeric"), length(Y) > 0,
                          msg = "Y should be a numeric vector of positive length")

  assertthat::assert_that(is.matrix(X), nrow(X) == length(Y),
                          msg = "X should be a matrix with nrow(X) == length(Y)")

  assertthat::assert_that(is(Z, "sparseMatrix"), ncol(Z) >= 1, nrow(Z) == length(Y),
                          msg = "Z should be a sparse matrix with nrow(Z) == length(Y)")

  assertthat::assert_that(is.list(Hlist), length(Hlist) >= 1,
                          all(sapply(Hlist, methods::is, "sparseMatrix")),
                          msg = "Hlist should be a list of sparse matrices")

  assertthat::assert_that(is.logical(REML), is.logical(expected),
                          is.logical(efficient), is.logical(signed),
                          msg = "REML, expected, efficient, and signed should all be logical")

  assertthat::assert_that(is.null(precomp) || is.list(precomp),
                          msg = "precomp should be NULL or a list")
    
  p <- ncol(X)
  # If max_radius not supplied, can use information matrix to make default
  if (length(max_radius) == 2) {
    lwr <- theta_null[test_idx] - max_radius[1]
    upr <- theta_null[test_idx] + max_radius[2]
  } else { # Assume length == 1
    lwr <- theta_null[test_idx] - max_radius
    upr <- theta_null[test_idx] + max_radius
  }
  if (lwr == upr) {
    null_values <- theta_null[test_idx]
    start_idx <- 1
  } else {
    null_values <- unique(sort(c(seq(lwr, upr, length.out = num_points),
      theta_null[test_idx])))
    start_idx <- which(null_values == theta_null[test_idx])
  }
  num_null <- length(null_values) # Can be 1, num_points, or num_points + 1
  # Start searching to the left of start_idx, including start_idx
  stat_vals <- rep(0, num_null)
  theta_tilde <- theta_null
  d <- length(theta_tilde)
  b <- if(!REML && p > 0) theta_tilde[1:p] else NULL
  if (is.null(precomp)) precomp <- get_precomp(Y = Y, X = X, Z = Z, b = b,
    REML = REML)

  # Evaluate test-statistic at null_values[start_idx - ii + 1]
  for(ii in 1:start_idx) {
    # Starting value for optimization is parameter vector with null
    # fixed and nuisance parameters at solutions at previous iteration
    # Should be valid for small enough step size.
    theta_tilde[test_idx] <- null_values[start_idx - ii + 1]
    theta_tilde <- maximize_loglik(start_val = theta_tilde,
                                    opt_idx = seq(d)[-test_idx],
                                    Y = Y,
                                    X = X,
                                    Z = Z,
                                    Hlist = Hlist,
                                    expected = TRUE,
                                    REML = REML,
                                    precomp = precomp)$arg
    # Update residual and relevant entries of precompute
    if(!REML && p > 0) {
      precomp$e <- Y - X %*% theta_tilde[1:p]
      precomp$Zte <- as.vector(crossprod(Z, precomp$e))
    }
    # Store the solution from start_idx to use when searching other direction
    if (ii == 1) {
      theta_tilde_start <- theta_tilde
    } 
    # Store test_statistic value
    stat_vals[start_idx - ii + 1] <- score_stat(theta = theta_tilde, 
                                                test_idx = test_idx, 
                                                Y = Y,
                                                X = X,
                                                Z = Z,
                                                Hlist = Hlist, 
                                                REML = REML,
                                                expected = expected,
                                                efficient = efficient,
                                                signed = signed,
                                                precomp = precomp) 
  }
  # Search to the right of start_idx
  if(start_idx < num_null) {
    theta_tilde <- theta_tilde_start
    # Evaluate test-statistic at null_values[ii], ii > start_idx
    for(ii in (start_idx + 1):num_null) {
      theta_tilde[test_idx] <- null_values[ii]
      theta_tilde <- maximize_loglik(start_val = theta_tilde,
                                    opt_idx = seq(d)[-test_idx],
                                    Y = Y,
                                    X = X,
                                    Z = Z,
                                    Hlist = Hlist,
                                    expected = TRUE,
                                    REML = REML,
                                    precomp = precomp)$arg
    # Update residual and relevant entries of precompute
    if(!REML && p > 0) {
      precomp$e <- Y - X %*% theta_tilde[1:p]
      precomp$Zte <- as.vector(crossprod(Z, precomp$e))
    }
    # Store test_statistic value
    stat_vals[ii] <- score_stat(theta = theta_tilde, 
                                test_idx = test_idx, 
                                Y = Y,
                                X = X,
                                Z = Z,
                                Hlist = Hlist, 
                                REML = REML,
                                expected = expected,
                                efficient = efficient,
                                signed = signed,
                                precomp = precomp) 
    }
  }
  # Return
  names(stat_vals) <- null_values
  stat_vals
}
