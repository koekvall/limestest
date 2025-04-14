partial_min <- function(opt_idx, precomp, psi_start, b = NULL, REML = TRUE,
                           expected = TRUE, ...)
{
  if(!is.null(b) & REML){
    warning("Coefficient vector supplied but not used by restricted likelihood")
    b <- NULL
  }
  H <- do.call(cbind, precomp$Hlist)
  r <- length(psi_start)
  #############################################################################
  # Define the objective function to be minimized
  #############################################################################
  obj_fun <- function(x){
    psi_arg <- psi_start
    psi_arg[opt_idx] <- x
    ll_things <- limestest:::loglikelihood(psi = psi_arg,
                                           b = b,
                                           precomp = precomp,
                                           REML = REML,
                                           expected = expected)
    list("value" = -ll_things$value, "gradient" = -ll_things$score,
         "hessian" = ll_things$infmat)
  }

  #############################################################################
  # Do minimization
  #############################################################################
  fit <- trust::trust(objfun = obj_fun, parinit = psi_start[opt_idx], rinit  = 1,
                      rmax = 100, ...)
  # Return results
  psi_start[opt_idx] <- fit$argument
  list("psihat" = psi_start, "value" = -fit$value, "conv" = fit$converged,
       "iter" = fit$iterations)
}

Psi_from_Hlist <- function(psi, Hlist)
{
  for(ii in seq_len(length(Hlist))){
    Hlist[[ii]] <- Hlist[[ii]] * psi[ii]
  }
  do.call(cbind, Hlist)
}

score_stat <- function(psi, test_idx, precomp, REML = TRUE, expected = TRUE,
                       efficient = TRUE, signed = FALSE)
{
  if(!expected & REML){
    warning("Observed information not available for restricted likelihood; using
            expected.")
  }
  test_idx <- unique(test_idx)
  r <- length(psi)
  k <- length(test_idx)
  stopifnot(all(test_idx %in% seq_len(r)), k <= r)

  psi_r <- psi[r]
  psi <- psi[-r]

  Psi_r <- (1 / psi_r) * Psi_from_Hlist(psi = psi, Hlist = precomp$Hlist)
  H <- do.call(cbind, precomp$Hlist)

  if(REML){
    ll_things <- res_ll(XtX = precomp$XtX,
                        XtY = precomp$XtY,
                        XtZ = precomp$XtZ,
                        ZtZ = precomp$ZtZ,
                        YtZ = precomp$YtZ,
                        Y = precomp$Y,
                        X = precomp$X,
                        Z = precomp$Z,
                        H = H,
                        Psi_r = Psi_r,
                        psi_r = psi_r,
                        lik = FALSE,
                        score = TRUE,
                        finf = TRUE)
  } else{
    ll_things <- loglik_psi(Z = precomp$Z,
                            ZtZXe = cbind(precomp$ZtZ, precomp$ZtX, precomp$Y),
                            e = precomp$Y,
                            H = H,
                            Psi_r = Psi_r,
                            psi_r = psi_r,
                            loglik = FALSE,
                            score = TRUE,
                            finf = TRUE,
                            expected = expected)
  }


  infmat <- ll_things$finf[test_idx, test_idx]

  if(efficient & (k > 1)){
    infmat <- infmat - crossprod(ll_things$finf[-test_idx, test_idx],
                                 solve(ll_things$finf[-test_idx, -test_idx],
                                       ll_things$finf[-test_idx, test_idx]))
  }

  if(signed){
    ed <- eigen(infmat, symmetric = TRUE)
    infmat <- ed$vectors %*% (sqrt(ed$values) * t(ed$vectors))

    # Inefficient since eigen decomposition available, but may be more stable
    test_stat <- solve(infmat, ll_things$score[test_idx])
  } else{
    test_stat <- crossprod(ll_things$score[test_idx], solve(infmat,
                                                           ll_things$score[test_idx]))
  }
  test_stat
}

#' @export
uni_test_stat <- function(test_seq, test_idx, psi, psi_r, Z, ZtZXe, e, Hlist)
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
                     Psi_r = Psi / psi_r,
                     psi_r = psi_r, finf = FALSE)
      psi[-test_idx] <- psi[-test_idx] +
        solve(score_inf$finf[-test_idx. -test_idx], s)
    }

    Psi <- Psi_from_Hlist(psi, Hlist)
    score_inf <- score_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                           Psi_r = Psi / psi_r,
                           psi_r = psi_r, finf = TRUE)

    eff_inf <- score_inf$finf[test_idx, test_idx] -
      sum(solve(score_inf$finf[-test_idx, -test_idx],
                score_inf$finf[-test_idx, test_idx]) *
            score_inf$finf[-test_idx, test_idx])

    test_stat[ii] <- score_inf$score[test_idx]^2 / eff_inf
  }
  test_stat
}
