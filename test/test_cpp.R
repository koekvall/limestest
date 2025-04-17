
test.R <- loglik_psi(Z = Z, ZtZXe = ZtZXe, e = e, H = H,
                    Psi0 = Psi_test / psi0_test, psi0 = psi0_test, loglik = T,
                    score = T, finf = T, expected = T)
test.Rcpp <- limestest:::loglik_psiRcpp(Z = Z, e = e, H = H,
                     Psi0 = Psi_test / psi0_test, psi0 = psi0_test, loglik = T,
                     score = T, finf = T, expected = T)

restest.R <- res_ll(XtX = crossprod(X), XtY = crossprod(X,y), XtZ = crossprod(X,Z),
                    ZtZ = crossprod(Z), YtZ = crossprod(y,Z),
                    Y=y, X=X, Z = Z, H = H, Psi0 = Psi_test / psi0_test, psi0 = psi0_test, lik = T,
                    score = T, finf = T)

restest.Rcpp <- limestest:::res_llRcpp(X=X, Y=y, Z=Z, H = H,
                            Psi0 = Psi_test / psi0_test, psi0 = psi0_test, lik = T,
                            score = T, finf = T)

all.equal(restest.R$ll, restest.Rcpp$ll)
all.equal(restest.R$score, restest.Rcpp$score)
all.equal(as.vector(restest.R$beta), restest.Rcpp$beta)
all.equal(matrix(restest.R$I_b_inv_chol@x, nrow=5), restest.Rcpp$I_b_inv_chol)
