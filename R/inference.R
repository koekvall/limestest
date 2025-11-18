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
      b <- if(p == 0) NULL else theta[1:p]
      e <- if(p == 0) Y else Y - X %*% b
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
  list("arg" = start_val, "value" = -fit$value, "conv" = fit$converged,
       "iter" = fit$iterations)
}


score_stat <- function(psi, test_idx, b = NULL, Y, X, Z, Hlist, REML = TRUE,
                       expected = TRUE, efficient = TRUE, signed = FALSE,
                       precomp = NULL)
{
  if(!expected && REML){
    warning("Observed information not available for restricted likelihood; using
            expected.")
  }

  test_idx <- unique(test_idx)
  r <- length(psi)
  k <- length(test_idx)
  stopifnot(all(test_idx %in% seq_len(r)), k <= r)

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
                             get_beta = FALSE,
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
