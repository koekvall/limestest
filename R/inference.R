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

#' Maximize log-likelihood with respect to specified parameters
#'
#' Optimizes a subset of parameters in a linear mixed effects model while
#' holding others fixed. Uses the trust region algorithm for optimization.
#'
#' @param start_val Numeric vector of starting parameter values. If
#'   \code{REML = FALSE}, this should be of length \eqn{p + r} with
#'   \code{start_val = c(beta, psi)}, where \eqn{p} is the number of fixed
#'   effects and \eqn{r} is the number of variance parameters. If
#'   \code{REML = TRUE}, this should be of length \eqn{r} with
#'   \code{start_val = psi}.
#' @param opt_idx Integer vector specifying which elements of \code{start_val}
#'   to optimize. All other parameters are held fixed at their starting values.
#'   Must have length > 0 and contain unique positive integers not exceeding
#'   \code{length(start_val)}.
#' @param Y Numeric vector of length \eqn{n} containing the response values.
#' @param X Numeric matrix of size \eqn{n \times p} containing fixed effect
#'   predictors.
#' @param Z Sparse matrix of size \eqn{n \times q} containing the random effect
#'   design matrix.
#' @param Hlist List of sparse matrices determining how \eqn{\psi} is mapped to
#'   the covariance matrix \eqn{\Psi} (see \code{?loglikelihood}).
#' @param expected Logical. If \code{TRUE}, use expected Fisher information
#'   matrix; otherwise use observed information. Default is \code{TRUE}.
#' @param REML Logical. If \code{TRUE}, use restricted maximum likelihood;
#'   otherwise use maximum likelihood. Default is \code{TRUE}.
#' @param precomp List or \code{NULL} containing precomputed quantities to speed
#'   up computation (see \code{?get_precomp}). If \code{NULL}, quantities are
#'   computed internally.
#' @param ... Additional arguments passed to \code{\link[trust]{trust}} optimizer,
#'   such as tolerance settings or iteration limits.
#'
#' @return A list with components:
#'   \item{arg}{Numeric vector of optimized parameter values. Parameters not in
#'     \code{opt_idx} retain their starting values from \code{start_val}.}
#'   \item{value}{The maximized log-likelihood value.}
#'   \item{conv}{Logical indicating whether the optimization converged.}
#'   \item{iter}{Integer giving the number of iterations performed.}
#'
#' @details
#' This function optimizes a subset of parameters in a linear mixed effects model
#' using the trust region algorithm from the \code{trust} package.
#'
#' The optimization uses both gradient (score) and Hessian (information matrix)
#' information for efficient convergence. A warning is issued if the optimization
#' does not converge.
#'
#'
#' @seealso \code{\link{score_stat}}, \code{\link{score_nuisance}},
#'   \code{\link{loglikelihood}}, \code{\link[trust]{trust}}
maximize_loglik <- function(start_val, opt_idx, Y, X, Z, Hlist, expected = TRUE,
                             REML = TRUE, precomp = NULL, ...) {
  # Argument checking
  assertthat::assert_that(is.vector(start_val, mode = "numeric"), length(start_val) > 0,
                          msg = "start_val should be a numeric vector of positive length")
  
  assertthat::assert_that(is.vector(opt_idx, mode = "numeric"), length(opt_idx) > 0,
                          all(opt_idx == floor(opt_idx)), all(opt_idx > 0),
                          msg = "opt_idx should be a vector of positive integers with length > 0")
  
  assertthat::assert_that(is.vector(Y, mode = "numeric"), length(Y) > 0,
                          msg = "Y should be a numeric vector of positive length")
  n <- length(Y)
  
  assertthat::assert_that(is.matrix(X), nrow(X) == n,
                          msg = "X should be a matrix with nrow(X) == length(Y)")
  
  assertthat::assert_that(is(Z, "sparseMatrix"), nrow(Z) == n, ncol(Z) > 0,
                          msg = "Z should be a sparse matrix with nrow(Z) == length(Y) and ncol(Z) > 0")
  
  assertthat::assert_that(is.list(Hlist), length(Hlist) > 0,
                          all(sapply(Hlist, methods::is, "sparseMatrix")),
                          msg = "Hlist should be a list of sparse matrices")
  
  assertthat::assert_that(is.logical(expected), length(expected) == 1,
                          msg = "expected should be a single logical value")
  
  assertthat::assert_that(is.logical(REML), length(REML) == 1,
                          msg = "REML should be a single logical value")
  
  assertthat::assert_that(is.null(precomp) || is.list(precomp),
                          msg = "precomp should be NULL or a list")
  
  r <- length(Hlist) + 1
  p <- ncol(X)
  expected_length <- if(REML) r else p + r
  assertthat::assert_that(length(start_val) == expected_length,
                          msg = paste0("start_val should have length ", expected_length,
                                      " (", if(REML) "r" else "p + r", 
                                      " for REML = ", REML, ")"))
  
  assertthat::assert_that(max(opt_idx) <= length(start_val),
                          msg = "opt_idx values must not exceed length(start_val)")
  
  assertthat::assert_that(length(unique(opt_idx)) == length(opt_idx),
                          msg = "opt_idx should not contain duplicate values")
  
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
      e <- if (p == 0) Y else Y - X %*% theta[1:p]
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
  
  if (!fit$converged) {
    warning("Optimization did not converge. Results may be unreliable. ",
            "Iterations: ", fit$iterations)
  }
  
  # Return results
  start_val[opt_idx] <- fit$argument
  names(start_val) <- if(REML) paste0("psi", 1:r) else c(paste0("b", 1:p), paste0("psi", 1:r))
  list("arg" = start_val, "value" = -fit$value, "conv" = fit$converged,
       "iter" = fit$iterations)
}

#' Score test statistic
#'
#' Computes the score test statistic for testing hypotheses about parameters
#' in a linear mixed effects model. The test statistic can be computed with or
#' without efficient information (accounting for nuisance parameters), and can
#' return either the quadratic form (unsigned) or the signed root statistic.
#'
#' @param theta Numeric vector of parameter values at which to evaluate the
#'   score test statistic. If \code{REML = FALSE}, this should be of length
#'   \eqn{p + r} with \code{theta = c(beta, psi)}, where \eqn{p} is the number
#'   of fixed effects and \eqn{r} is the number of variance parameters. If
#'   \code{REML = TRUE}, this should be of length \eqn{r} with \code{theta = psi}.
#' @param test_idx Integer vector specifying which elements of \code{theta}
#'   are being tested. These are the parameters constrained by the null hypothesis.
#' @param Y Numeric vector of length \eqn{n} containing the response values.
#' @param X Numeric matrix of size \eqn{n \times p} containing fixed effect
#'   predictors.
#' @param Z Sparse matrix of size \eqn{n \times q} containing the random effect
#'   design matrix.
#' @param Hlist List of sparse matrices determining how \eqn{\psi} is mapped to
#'   the covariance matrix \eqn{\Psi} (see \code{?loglikelihood}).
#' @param REML Logical. If \code{TRUE}, use restricted maximum likelihood;
#'   otherwise use maximum likelihood. Default is \code{TRUE}.
#' @param expected Logical. If \code{TRUE}, use expected Fisher information
#'   matrix; otherwise use observed information. Default is \code{TRUE}.
#'   Note: Observed information is not available for REML.
#' @param efficient Logical. If \code{TRUE}, use efficient information that
#'   accounts for estimation of nuisance parameters. If \code{FALSE}, use
#'   the full information matrix without adjustment. Default is \code{TRUE}.
#' @param signed Logical. If \code{TRUE}, return the signed root statistic
#'   (vector). If \code{FALSE}, return the quadratic form statistic (scalar
#'   for single parameter tests). Default is \code{FALSE}.
#' @param known_idx Integer vector or \code{NULL} specifying which elements of
#'   \code{theta} (other than those in \code{test_idx}) have known values and
#'   should be held fixed (not treated as nuisance parameters to be profiled over).
#'   If \code{NULL} (default), all parameters not in \code{test_idx} are treated
#'   as nuisance parameters. Must not overlap with \code{test_idx}.
#' @param precomp List or \code{NULL} containing precomputed quantities to speed
#'   up computation (see \code{?get_precomp}). If \code{NULL}, all quantities
#'   are computed from scratch.
#'
#' @return Numeric value or vector containing the score test statistic. If
#'   \code{signed = FALSE}, returns a scalar (the quadratic form). If
#'   \code{signed = TRUE}, returns a vector (the signed root statistic).
#'   Under the null hypothesis, the squared statistic asymptotically follows
#'   a chi-squared distribution with degrees of freedom equal to
#'   \code{length(test_idx)}.
#'
#' @details
#' The score test statistic is computed as:
#' \deqn{T = S_t^T I_{tt}^{-1} S_t}
#' where \eqn{S_t} is the score vector for the test parameters and \eqn{I_{tt}}
#' is the information matrix (possibly efficient if \code{efficient = TRUE}).
#'
#' When \code{efficient = TRUE} and there are nuisance parameters (parameters
#' not in \code{test_idx} or \code{known_idx}), the efficient information is:
#' \deqn{I_{tt}^{eff} = I_{tt} - I_{tn} I_{nn}^{-1} I_{nt}}
#' where subscripts \eqn{t} denote test parameters and \eqn{n} denote nuisance
#' parameters.
#'
#' When \code{known_idx} is specified, those parameters are treated as fixed and
#' known (not as nuisance parameters). This is useful when some parameters have
#' been estimated separately or are constrained to specific values.
#'
#' The function performs checks for numerical stability, including the condition
#' number of the information matrix. Warnings are issued if potential numerical
#' problems are detected.
#'
#' @seealso \code{\link{score_nuisance}}, \code{\link{loglikelihood}}
#'
#' @export
score_stat <- function(theta, test_idx, Y, X, Z, Hlist, REML = TRUE,
                       expected = TRUE, efficient = TRUE, signed = FALSE, known_idx = NULL,
                       precomp = NULL)
{
  # Argument checking
  assertthat::assert_that(is.vector(theta, mode = "numeric"), length(theta) > 0,
                          msg = "theta should be a numeric vector of positive length")
  
  assertthat::assert_that(is.vector(test_idx, mode = "numeric"), length(test_idx) > 0,
                          all(test_idx == floor(test_idx)), all(test_idx > 0),
                          msg = "test_idx should be a vector of positive integers")
  
  assertthat::assert_that(is.null(known_idx) || 
                          (is.vector(known_idx, mode = "numeric") && length(known_idx) >= 0 &&
                           all(known_idx == floor(known_idx)) && all(known_idx > 0)),
                          msg = "known_idx should be NULL or a vector of positive integers")
  
  assertthat::assert_that(is.vector(Y, mode = "numeric"), length(Y) > 0,
                          msg = "Y should be a numeric vector of positive length")
  n <- length(Y)
  
  assertthat::assert_that(is.matrix(X), nrow(X) == n,
                          msg = "X should be a matrix with nrow(X) == length(Y)")
  
  assertthat::assert_that(is(Z, "sparseMatrix"), nrow(Z) == n, ncol(Z) > 0,
                          msg = "Z should be a sparse matrix with nrow(Z) == length(Y) and ncol(Z) > 0")
  
  assertthat::assert_that(is.list(Hlist), length(Hlist) > 0,
                          all(sapply(Hlist, methods::is, "sparseMatrix")),
                          msg = "Hlist should be a list of sparse matrices")
  
  assertthat::assert_that(is.logical(REML), length(REML) == 1,
                          msg = "REML should be a single logical value")
  
  assertthat::assert_that(is.logical(expected), length(expected) == 1,
                          msg = "expected should be a single logical value")
  
  assertthat::assert_that(is.logical(efficient), length(efficient) == 1,
                          msg = "efficient should be a single logical value")
  
  assertthat::assert_that(is.logical(signed), length(signed) == 1,
                          msg = "signed should be a single logical value")
  
  assertthat::assert_that(is.null(precomp) || is.list(precomp),
                          msg = "precomp should be NULL or a list")
  
  p <- ncol(X)
  r <- length(Hlist) + 1
  expected_length <- if(REML) r else p + r
  assertthat::assert_that(length(theta) == expected_length,
                          msg = paste0("theta should have length ", expected_length,
                                      " (", if(REML) "r" else "p + r", 
                                      " for REML = ", REML, ")"))
  
  assertthat::assert_that(max(test_idx) <= length(theta),
                          msg = "test_idx values must not exceed length(theta)")
  
  if (!is.null(known_idx)) {
    assertthat::assert_that(max(known_idx) <= length(theta),
                            msg = "known_idx values must not exceed length(theta)")
    
    assertthat::assert_that(length(intersect(test_idx, known_idx)) == 0,
                            msg = "test_idx and known_idx should not overlap")
  }
  
  if(!expected && REML){
    warning("Observed information not available for restricted likelihood; using
            expected.")
  }
  
  # Remove duplicates from index vectors
  test_idx <- unique(test_idx)
  if (!is.null(known_idx) && length(known_idx) > 0) {
    known_idx <- unique(known_idx)
  } else {
    known_idx <- NULL  # Treat empty vector as NULL
  }

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
    # Check condition of information matrix
    cond <- tryCatch(
      kappa(ll_things$inf_mat, exact = TRUE),
      error = function(e) Inf
    )
    if (cond > 1e12) {
      warning("Information matrix is poorly conditioned (condition number: ",
              format(cond, scientific = TRUE), "). Results may be unreliable.")
    }
  
  inf_mat <- ll_things$inf_mat[test_idx, test_idx, drop = FALSE]

  # Use efficient information only if there are nuisance parameters
  # Combine test_idx and known_idx to exclude from nuisance parameters
  exclude_idx <- if (is.null(known_idx)) test_idx else c(test_idx, known_idx)
  
  if (efficient && (length(ll_things$score) > length(exclude_idx))) {
    A_nt <- ll_things$inf_mat[-exclude_idx, test_idx, drop = FALSE]
    I_nn <- ll_things$inf_mat[-exclude_idx, -exclude_idx, drop = FALSE]
    inf_mat <- inf_mat - crossprod(A_nt, solve(I_nn, A_nt))
  }

  if (signed) {
    ed <- eigen(inf_mat, symmetric = TRUE)
    # Numerical floor to avoid taking sqrt of tiny/negative values
    ev <- pmax(ed$values, .Machine$double.eps)
    inf_root <- ed$vectors %*% (sqrt(ev) * t(ed$vectors))
    test_stat <- solve(inf_root, ll_things$score[test_idx])
  } else {
    test_stat <- crossprod(ll_things$score[test_idx], solve(inf_mat,
                                    ll_things$score[test_idx]))
  }
  as.vector(test_stat)
}

#' Score test statistic with nuisance parameters
#'
#' Computes the score test statistic over a range of values for a test parameter,
#' while optimizing nuisance parameters at each point. This is useful for
#' constructing confidence intervals and profile likelihood-based inference.
#'
#' @param theta_start Numeric vector of parameter values at which to center the
#'   search. Must be a valid parameter vector. If \code{REML = FALSE}, this
#'   should be of length \eqn{p + r} with \code{theta_start = c(beta, psi)},
#'   where \eqn{p} is the number of fixed effects and \eqn{r} is the number of
#'   variance parameters. If \code{REML = TRUE}, this should be of length \eqn{r}
#'   with \code{theta_start = psi}.
#' @param test_idx Integer specifying which element of \code{theta_start} to
#'   test. This parameter will be fixed at a range of values (null hypotheses)
#'   while other parameters are optimized.
#' @param max_radius Numeric value or vector of length 2 specifying the radius
#'   around \code{theta_start[test_idx]} to search. If length 1, the search is
#'   symmetric. If length 2, \code{max_radius[1]} specifies the lower radius and
#'   \code{max_radius[2]} specifies the upper radius. If 0 (default), only
#'   \code{theta_start[test_idx]} is evaluated.
#' @param num_points Integer specifying the number of points to evaluate in the
#'   range defined by \code{max_radius}. Default is 100.
#' @param Y Numeric vector of length \eqn{n} containing the response values.
#' @param X Numeric matrix of size \eqn{n \times p} containing fixed effect
#'   predictors.
#' @param Z Sparse matrix of size \eqn{n \times q} containing the random effect
#'   design matrix.
#' @param Hlist List of sparse matrices determining how \eqn{\psi} is mapped to
#'   the covariance matrix \eqn{\Psi} (see \code{?loglikelihood}).
#' @param REML Logical. If \code{TRUE}, use restricted maximum likelihood;
#'   otherwise use maximum likelihood. Default is \code{TRUE}.
#' @param expected Logical. If \code{TRUE}, use expected Fisher information
#'   matrix; otherwise use observed information. Default is \code{TRUE}.
#' @param efficient Logical. If \code{TRUE}, use efficient information that
#'   accounts for estimation of nuisance parameters. Default is \code{TRUE}.
#' @param signed Logical. If \code{TRUE}, return the signed root statistic;
#'   otherwise return the squared (chi-squared type) statistic. Default is
#'   \code{TRUE}.
#' @param known_idx Integer vector or \code{NULL} specifying which elements of
#'   \code{theta_start} (other than \code{test_idx}) have known values and should
#'   be held fixed during optimization. These parameters will not be treated as
#'   nuisance parameters. If \code{NULL} (default), all parameters except
#'   \code{test_idx} are optimized. Must not overlap with \code{test_idx}.
#' @param precomp List or \code{NULL} containing precomputed quantities to speed
#'   up computation (see \code{?get_precomp}). If \code{NULL}, all quantities
#'   are computed internally.
#' @param ... Additional arguments passed to \code{\link[trust]{trust}} optimizer
#'   used in \code{\link{maximize_loglik}}.
#'
#' @return Named numeric vector of score test statistics. Names correspond to the
#'   values of the test parameter at which the statistic was evaluated. Under the
#'   null hypothesis, the squared statistic asymptotically follows a chi-squared
#'   distribution with 1 degree of freedom.
#'
#' @details
#' This function performs profile-likelihood-based inference by fixing the test
#' parameter at various values and optimizing over nuisance parameters. For each
#' value in the specified range around \code{theta_start[test_idx]}, the function:
#' \enumerate{
#'   \item Fixes the test parameter at that value
#'   \item Fixes any known parameters specified in \code{known_idx}
#'   \item Optimizes the remaining (nuisance) parameters using
#'         \code{\link{maximize_loglik}}
#'   \item Computes the score test statistic at the resulting parameter values
#'         using \code{\link{score_stat}}
#' }
#'
#' The search proceeds in two directions from \code{theta_start[test_idx]}: first
#' decreasing values, then increasing values. This allows for efficient
#' warm-starting of the optimization at each step, using the previous solution
#' as the starting point.
#'
#' When \code{efficient = TRUE}, the score test accounts for uncertainty in the
#' nuisance parameters. When \code{known_idx} is specified, those parameters are
#' treated as fixed and known, not as nuisance parameters to be profiled over.
#' @export
score_nuisance <- function(theta_start, test_idx, max_radius = 0, num_points = 1e2,
                           Y, X, Z, Hlist, REML = TRUE, expected = TRUE,
                           efficient = TRUE, signed = TRUE, known_idx = NULL,
                           precomp = NULL, ...) {
  # Argument checking
  assertthat::assert_that(is.vector(theta_start, mode = "numeric"), length(theta_start) > 0,
                          msg = "theta_start should be a numeric vector of positive length")

  assertthat::assert_that(is.vector(test_idx, mode = "numeric"), length(test_idx) == 1,
                          test_idx >= 1, test_idx <= length(theta_start),
                          msg = "test_idx should be a single integer between 1 and length(theta_start)")

  assertthat::assert_that(is.vector(max_radius, mode = "numeric"), length(max_radius) %in% c(1, 2),
                          all(max_radius >= 0),
                          msg = "max_radius should be a numeric value or vector of length 1 or 2 with non-negative values")

  assertthat::assert_that(is.vector(num_points, mode = "numeric"), length(num_points) == 1,
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

  assertthat::assert_that(is.null(known_idx) || 
                          (is.numeric(known_idx) && length(known_idx) >= 0 &&
                           all(known_idx == floor(known_idx)) && all(known_idx > 0)),
                          msg = "known_idx should be NULL or a vector of positive integers")

  assertthat::assert_that(is.null(precomp) || is.list(precomp),
                          msg = "precomp should be NULL or a list")

  p <- ncol(X)
  r <- length(Hlist) + 1
  expected_length <- if(REML) r else p + r
  assertthat::assert_that(length(theta_start) == expected_length,
                          msg = paste0("theta_start should have length ", expected_length,
                                      " (", if(REML) "r" else "p + r", 
                                      " for REML = ", REML, ")"))

  if (!is.null(known_idx)) {
    assertthat::assert_that(max(known_idx) <= length(theta_start),
                            msg = "known_idx values must not exceed length(theta_start)")
    
    assertthat::assert_that(!(test_idx %in% known_idx),
                            msg = "test_idx and known_idx should not overlap")
  }

  p <- ncol(X)
  # If max_radius not supplied, could use information matrix to make default
  if (length(max_radius) == 2) {
    lwr <- theta_start[test_idx] - max_radius[1]
    upr <- theta_start[test_idx] + max_radius[2]
  } else { # Assume length == 1
    lwr <- theta_start[test_idx] - max_radius
    upr <- theta_start[test_idx] + max_radius
  }
  if (lwr == upr) {
    null_values <- theta_start[test_idx]
    start_idx <- 1
  } else {
    null_values <- unique(sort(c(seq(lwr, upr, length.out = num_points),
      theta_start[test_idx])))
    start_idx <- which(null_values == theta_start[test_idx])
  }
  num_null <- length(null_values) # Can be 1, num_points, or num_points + 1
  # Start searching to the left of start_idx, including start_idx
  stat_vals <- rep(0, num_null)
  theta_tilde <- theta_start
  d <- length(theta_tilde)
  b <- if (!REML && p > 0) theta_tilde[1:p] else NULL
  if (is.null(precomp)) precomp <- get_precomp(Y = Y, X = X, Z = Z, b = b,
    REML = REML)

  # Determine which parameters to optimize (exclude test and known parameters)
  exclude_idx <- if (is.null(known_idx)) test_idx else c(test_idx, known_idx)
  opt_idx <- seq(d)[-exclude_idx]
  
  # Evaluate test-statistic at null_values[start_idx - ii + 1]
  for (ii in 1:start_idx) {
    # Starting value for optimization is parameter vector with null
    # fixed and nuisance parameters at solutions at previous iteration
    # Should be valid for small enough step size.
    theta_tilde[test_idx] <- null_values[start_idx - ii + 1]
    theta_tilde <- maximize_loglik(start_val = theta_tilde,
                                   opt_idx = opt_idx,
                                   Y = Y,
                                   X = X,
                                   Z = Z,
                                   Hlist = Hlist,
                                   expected = expected,
                                   REML = REML,
                                   precomp = precomp, ...)$arg
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
                                                known_idx = known_idx,
                                                precomp = precomp)
  }
  # Search to the right of start_idx
  if(start_idx < num_null) {
    theta_tilde <- theta_tilde_start
    if(!REML && p > 0) {
      precomp$e <- Y - X %*% theta_tilde[1:p]
    }
    # Evaluate test-statistic at null_values[ii], ii > start_idx
    for(ii in (start_idx + 1):num_null) {
      theta_tilde[test_idx] <- null_values[ii]
      theta_tilde <- maximize_loglik(start_val = theta_tilde,
                                    opt_idx = opt_idx,
                                    Y = Y,
                                    X = X,
                                    Z = Z,
                                    Hlist = Hlist,
                                    expected = expected,
                                    REML = REML,
                                    precomp = precomp, ...)$arg
    # Update residual and relevant entries of precompute
    if(!REML && p > 0) {
      precomp$e <- Y - X %*% theta_tilde[1:p]
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
                                known_idx = known_idx,
                                precomp = precomp)
    }
  }
  # Return
  names(stat_vals) <- null_values
  stat_vals
}
