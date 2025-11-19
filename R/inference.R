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


score_nuisance <- function(theta_null, test_idx, max_radius = 0, num_points = 1e2,
                           Y, X, Z, Hlist, REML = TRUE, expected = TRUE, efficient = TRUE, signed = TRUE,
                           precomp = NULL) {
    # Insert argument checking here

    p <- ncol(X)
    n <- nrow(X)
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
