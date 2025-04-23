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


partial_min_psi <- function(psi_start, opt_idx, b = NULL, Y, X, Z, Hlist,
                        REML = TRUE, expected = TRUE, precomp = NULL, ...)
{
  if(!is.null(b) & REML){
    warning("Coefficient vector supplied but not used by restricted likelihood")
    b <- NULL
  }
  #############################################################################
  # Define the objective function to be minimized
  #############################################################################
  obj_fun <- function(x){
    psi_arg <- psi_start
    psi_arg[opt_idx] <- x
    ll_things <- loglikelihood(psi = psi_arg,
                               b = b,
                               Y = Y,
                               X = X,
                               Z = Z,
                               Hlist = Hlist,
                               REML = REML,
                               get_val = TRUE,
                               get_score = TRUE,
                               get_inf = TRUE,
                               get_beta = FALSE,
                               expected = expected,
                               precomp = precomp)
    list("value" = -ll_things$value, "gradient" = -ll_things$score[opt_idx],
         "hessian" = as.matrix(ll_things$inf_mat[opt_idx, opt_idx]))
  }

  #############################################################################
  # Do minimization
  #############################################################################
  browser()
  fit <- trust::trust(objfun = obj_fun, parinit = psi_start[opt_idx], rinit  = 1,
                      rmax = 100, ...)
  # Return results
  psi_start[opt_idx] <- fit$argument
  list("psi_hat" = psi_start, "value" = -fit$value, "conv" = fit$converged,
       "iter" = fit$iterations)
}

score_stat <- function(psi, test_idx, b = NULL, Y, X, Z, Hlist, REML = TRUE,
                       expected = TRUE, efficient = TRUE, signed = FALSE,
                       precomp = NULL)
{
  if(!expected & REML){
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

  if(efficient & (r > 1)){
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
