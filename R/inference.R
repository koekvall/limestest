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

  psi0 <- psi[r]
  psi <- psi[-r]

  Psi0 <- (1 / psi0) * Psi_from_Hlist(psi = psi, Hlist = precomp$Hlist)
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
                        Psi0 = Psi0,
                        psi0 = psi0,
                        lik = FALSE,
                        score = TRUE,
                        finf = TRUE)
  } else{
    ll_things <- loglik_psi(Z = precomp$Z,
                            ZtZXe = cbind(precomp$ZtZ, precomp$ZtX, precomp$Y),
                            e = precomp$Y,
                            H = H,
                            Psi0 = Psi0,
                            psi0 = psi0,
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
    test_stat <- solve(infmat, ll_things$score[test_idx]) # Inefficient but stable
  } else{
    test_stat <- crossprod(ll_things$score[test_idx], solve(infmat,
                                                           ll_things$score[test_idx]))
  }
  test_stat
}

#' @export
uni_test_stat <- function(test_seq, test_idx, psi, psi0, Z, ZtZXe, e, Hlist)
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
                     Psi0 = Psi / psi0,
                     psi0 = psi0, finf = FALSE)
      psi[-test_idx] <- psi[-test_idx] +
        solve(score_inf$finf[-test_idx. -test_idx], s)
    }

    Psi <- Psi_from_Hlist(psi, Hlist)
    score_inf <- score_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                           Psi0 = Psi / psi0,
                           psi0 = psi0, finf = TRUE)

    eff_inf <- score_inf$finf[test_idx, test_idx] -
      sum(solve(score_inf$finf[-test_idx, -test_idx],
                score_inf$finf[-test_idx, test_idx]) *
            score_inf$finf[-test_idx, test_idx])

    test_stat[ii] <- score_inf$score[test_idx]^2 / eff_inf
  }
  test_stat
}
