#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// Compute projection onto intersection of set with restricted elements
// and set of symmetric matrices with eigenvalues bounded below.
//
// Arguments:
//
// X: symmetric matrix to project
// restr_idx: vector of indexes of X to restrict (in column-major order)
// restr: vector of values to restrict elements of X to (column major order)
// eps: scalar (double) lower bound on the eigenvalues of matrices in the set
//   projected on.
// tol: scalar (double) tolerance for terminating the algorithm.
// maxit: positive integer with maximum number of iterations for algorithm.
arma::mat project_rcpp(arma::mat X, const arma::uvec restr_idx,
                       const arma::vec restr, const double eps, const double tol, uint maxit)
{
  const uint d = X.n_cols;
  const double Inf = std::numeric_limits<double>::infinity();

  arma::mat Xk(d, d);
  if(restr_idx.n_elem == 0){
    arma::vec e(d);
    X = arma::symmatu(X);
    arma::eig_sym(e, X, X);
    Xk = X * arma::diagmat(arma::clamp(e, eps, Inf)) * X.t();
  } else{
    arma::mat Xkm1 = X;
    arma::mat Pkm1(d, d, arma::fill::zeros);
    arma::mat Qkm1(d, d, arma::fill::zeros);
    double obj_old = 0.0;
    for (size_t kk = 0; kk < maxit; kk++) {
      if(kk >= maxit - 1){
        Rcpp::warning("Projection algorithm reached max. iter. \n");
      }
      arma::mat Ykm1 = Xkm1 + Pkm1;
      Ykm1.elem(restr_idx - 1) = restr; // C++ index start at zero
      arma::mat Pk = Xkm1 + Pkm1 - Ykm1;
      arma::mat U = Ykm1 + Qkm1;
      arma::vec e(d);
      U = arma::symmatu(U);
      arma::eig_sym(e, U, U);
      Xk = U * arma::diagmat(arma::clamp(e, eps, Inf)) * U.t();
      arma::mat Qk = Ykm1 + Qkm1 - Xk;
      double obj_new = 0.5 * arma::accu(arma::square((Xk - X)));
      if(std::abs(obj_new - obj_old) < tol){
        break;
      }
      // prepare next iteration
      Xkm1 = Xk;
      Pkm1 = Pk;
      Qkm1 = Qk;
      obj_old = obj_new;
    }
  }
  return Xk;
}

//' Psi_from_H_cpp
//'
//' Constructs the covariance matrix of the random effects
//'
//' @param psi_mr A vector of
//' @param H Matrix of derivatives of Psi with respect to elements of psi.
//'        Assumes H = [H_1, ... , H_r], where H_j is q by q.
//' @return The covariance matrix Psi
// [[Rcpp::export]]
Eigen::SparseMatrix<double> Psi_from_H_cpp(const Eigen::Map<Eigen::VectorXd> psi_mr,
                                           const Eigen::MappedSparseMatrix<double> H) { //
  int q = H.rows();
  int rm1 = H.cols() / q;
  Eigen::SparseMatrix<double> Psi(q, q);
  for (int ii = 0; ii < rm1; ii++) {
     Psi += psi_mr(ii) * H.middleCols(ii * q, q);
  }
  return Psi;
}

//' loglik_psi_cpp
//'
//' Calculates the log-likelihood, score vector and Fisher information matrix
//' for the variance parameter vector \code{psi} in a linear mixed effects model.
//' This function is implemented in C++ and is faster than the equivalent function \code{loglik_psi}
//'
//' @param ZtZ A matrix of fixed effects.
//' @param e The residual vector.
//' @param H Matrix of derivatives of Psi with respect to elements of psi.
//'        Assumes H = [H_1, ... , H_r], where H_j is q by q.
//' @param Psi_r Coavariance matrix of random effects (Psi) divided by error
//'        variance psi_r, with dimensions q by q.
//' @param psi_r The error variance
//' @param loglik If \code{TRUE} (default), the log-likelihood will be calculated.
//' @param score If \code{TRUE} (default), the score vector will be calculated.
//' @param finf If \code{TRUE} (default), the Fisher information matrix will be calculated.
//' @param expected If \code{TRUE} (detault), return expected information;
//'        otherwise observed.
//'
//' @return A list with components:
//' \item{ll}{The log-likelihood.}
//' \item{score}{The score vector.}
//' \item{finf}{The Fisher information matrix.}
//'
//' @useDynLib limestest, .registration=TRUE
// [[Rcpp::export]]

Rcpp::List loglik_psi_cpp(const Eigen::MappedSparseMatrix<double> ZtZ,
                          const Eigen::Map<Eigen::MatrixXd> XtZ,
                          const Eigen::Map<Eigen::VectorXd> Zte,
                          const Eigen::MappedSparseMatrix<double> Z,
                          Eigen::VectorXd e,
                          Eigen::SparseMatrix<double> H,
                          const Eigen::MappedSparseMatrix<double> Psi_r,
                          const double psi_r,
                          const bool get_val = true,
                          const bool get_score = true,
                          const bool get_inf = true,
                          const bool expected = true) {

  // Define dimensions
  int n = e.size();
  int q = Psi_r.cols();
  int rm1 = H.cols() / q;
  int r = rm1 + 1;
  // loglikelihood to return
  double ll = NA_REAL;
  // Score vector to return
  Eigen::VectorXd s_psi(r);
  // Fisher information to return
  Eigen::MatrixXd I_psi(r, r);

  Eigen::SparseMatrix<double> A = Psi_r * ZtZ;
  Eigen::VectorXd Ae = Psi_r * Zte;

  // solver for Psi_rZtZ + I_q
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

  Eigen::SparseMatrix<double> Id_q(q, q);
  Id_q = Eigen::MatrixXd::Identity(q, q).sparseView();
  solver.compute(A + Id_q);

  // Add loglik term before overwriting
  if (get_val) {
    ll = -0.5 * solver.logAbsDeterminant() - 0.5 * n * log(psi_r);
  }
  // Matrix denoted M in manuscript is A[, 1:q]
  A = solver.solve(A);
  Ae = solver.solve(Ae);

  Eigen::VectorXd e_save(n);
  if(get_val | (get_inf & !expected)){
     e_save = e;
  }
  // update e
  e = (1.0 / psi_r) * (e - Z * Ae); // = Sigma^{-1}e
  if (get_val) {
    ll = ll - 0.5 * e.dot(e_save);
  }
  double trace_M = A.diagonal().sum();

  s_psi(rm1) = 0.5 * e.dot(e) - (0.5 / psi_r) * (n - trace_M);

  // v is Z'Sigma^{-1}e, w is e'Sigma^{-1}Z[H1...Hr]
  Eigen::VectorXd v = (Z.transpose() * e);
  Eigen::VectorXd w = H.transpose() * v;

  for (int ii = 0; ii < rm1; ii++) {
    s_psi(ii) = 0.5 * w.middleRows(ii * q, q).dot(v);
  }

  // B = Z'Z (M - I_q) in paper notation, sparse
  Eigen::SparseMatrix<double> B = A;
  B.diagonal().array() -= 1;// -= Eigen::VectorXd::Constant(q,1);
  B = ZtZ * B;

  if (!get_inf) {
    // Compute -[ZtZ (I_q - M) * H_1, ..., ZtZ (I_q - M) * H_r] using recycling
    // The "if" is because the calculation is a byproduct of a more expensive one
    // (B %*% H) done to get Fisher information
    for (int ii = 0; ii < rm1; ii++) {
      s_psi(ii) += (0.5 / psi_r) * H.middleCols(ii * q, q).cwiseProduct(B).sum();
    }
  } else {
    H = B * H;
    I_psi(rm1, rm1) = (0.5 / (psi_r * psi_r)) * (n - 2 * trace_M +
      Eigen::SparseMatrix<double>(A.transpose()).cwiseProduct(A).sum());

    // M = M-I
    A.diagonal().array() -= 1;
    Eigen::SparseMatrix<double> H1(q, q);
    Eigen::SparseMatrix<double> H2(q, q);
    for (int ii = 0; ii < rm1; ii++) {
      H1 = H.middleCols(ii * q, q);
      I_psi(rm1, ii) = I_psi(ii, rm1) = (0.5 / (psi_r * psi_r)) *
        A.cwiseProduct(H1).sum();
        s_psi(ii) += (0.5 / psi_r) * H1.diagonal().sum();
      for (int jj = ii; jj < rm1; jj++) {
        H2 = H.middleCols(jj * q, q).transpose();
        I_psi(ii, jj) = I_psi(jj, ii) = (0.5 / (psi_r * psi_r)) * H1.cwiseProduct(H2).sum();
      }
    }
    if (!expected) {
      I_psi *= -1.0;
      // u = Sigma^{-2}e. Some calculations could be saved from before
      Eigen::VectorXd u = (1.0 / (psi_r * psi_r)) * (e_save + Z * (-2.0 * Ae + (A + Id_q) * Ae));
      I_psi(rm1, rm1) += e.dot(u);
      v = Z.transpose() * u;

      for (int ii = 0; ii < rm1; ii++) {
        I_psi(ii, rm1) +=  w.middleRows(ii * q, q).dot(v);
        I_psi(rm1, ii) = I_psi(ii, rm1);
        for (int jj = ii; jj < rm1; jj++) {
          I_psi(ii, jj) -= (1 / psi_r) * (B * w.middleRows(ii * q, q)).dot(w.middleRows(jj * q, q));
          I_psi(jj, ii) = I_psi(ii, jj);
        }
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("value") = ll,
                            Rcpp::Named("score") = s_psi,
                            Rcpp::Named("inf_mat") = I_psi);
}


//' Compute Restricted Likelihood, Score and Information
//'
//' Computes the restricted (residual) likelihood, score, and information matrix
//' for the variance parameter \code{psi} a linear mixed effects model
//'
//' @param X An n x p matrix of the design matrix of fixed effects.
//' @param Y An n x 1 vector of the response variable.
//' @param Z An n x q matrix of the design matrix of random effects.
//' @param H A q x rq matrix, where H = [H_1, ..., H_r], with H_j being the
//' derivative of Psi with respect to psi_j
//' @param Psi_r The covariance matrix of the random effects (Psi) divided by the
//' error variance (psi_r)
//' @param psi_r A scalar value of the error variance.
//' @param lik If \code{TRUE} (default), the log-likelihood will be computed.
//' @param score If \code{TRUE} (default), the score vector will be computed.
//' @param finf If \code{TRUE} (default), the Fisher information matrix will be
//' computed.
//'
//' @return A list with elements:
//' \describe{
//' \item{ll}{A scalar value of the restricted log-likelihood.}
//' \item{score}{A (r + 1) x 1 vector of the restricted score of the variance parameters.}
//' \item{finf}{A (r + 1) x (r + 1) matrix of the restricted  information of the
//' variance parameters.}
//' \item{beta}{Maximum likelihood estimate for the fixed effects parameter}
//' \item{I_b_inv_chol}{Inverse of the Cholesky root of the information matrix for the fixed effects parameter beta}
//' }
//' @import Matrix
// [[Rcpp::export]]
Rcpp::List res_ll_cpp(Eigen::VectorXd Y,
                      const Eigen::Map<Eigen::MatrixXd> X,
                      const Eigen::MappedSparseMatrix<double> Z,
                      const Eigen::Map<Eigen::MatrixXd> XtY,
                      const Eigen::Map<Eigen::MatrixXd> ZtY,
                      const Eigen::Map<Eigen::MatrixXd> XtX,
                      const Eigen::Map<Eigen::MatrixXd> XtZ,
                      const Eigen::MappedSparseMatrix<double> ZtZ,
                      Eigen::SparseMatrix<double> H,
                      const Eigen::MappedSparseMatrix<double> Psi_r,
                      const double psi_r,
                      const bool get_val = true,
                      const bool get_score = true,
                      const bool get_inf = true)
{
;
  // Define dimensions
  int n = Y.size();
  int q = Psi_r.cols();
  int rm1 = H.cols() / q;
  int r = rm1 + 1;
  int p = X.cols();
  // loglikelihood to return
  double ll = NA_REAL;
  // Score vector to return
  Eigen::VectorXd s_psi(r);
  // Information matrix to return
  Eigen::MatrixXd I_psi(r, r);


  // Pre-compute (I_q + Psi_r Z'Z)^{-1} Psi_r
  // solver for Psi_rZtZ + I_q
  Eigen::SparseMatrix<double> Id_q(q, q);
  Id_q = Eigen::MatrixXd::Identity(q, q).sparseView();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(Psi_r * ZtZ + Id_q);

  if (get_val) {
    ll = solver.logAbsDeterminant();
  }
  Eigen::SparseMatrix<double> A = solver.solve(Psi_r);
  Eigen::MatrixXd B = XtZ * A;

  //Create XtSiX
  Eigen::MatrixXd U = (1.0 / psi_r) * (XtX - B * XtZ.transpose());  //p*p

  // Force symmetric
   U = U.selfadjointView<Eigen::Upper>();
  // llt decomposition
  Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> llt(U);

  if (llt.info() != Eigen::Success) {
    return Rcpp::List::create(
      Rcpp::Named("ll") = -R_PosInf,
      Rcpp::Named("score") = s_psi,
      Rcpp::Named("finf") = I_psi,
      Rcpp::Named("beta") = Eigen::VectorXd::Constant(p, NA_REAL),
      Rcpp::Named("I_b_inv_chol") = Eigen::MatrixXd::Constant(p, p, NA_REAL));
  }

  // Create XtSiY
  // (1/ psi_r) * (XtY - XtZ %*% A %*% ZtY)
  Eigen::VectorXd beta_tilde = (1.0 / psi_r) * (XtY -
    XtZ * (A * ZtY));
  beta_tilde = llt.solve(beta_tilde);

  // Replace response with residuals
  Y = Y - X * beta_tilde;

  // n x 1 vector for storing \Sigma^{-1}e
  Eigen::VectorXd a = (1.0 / psi_r) * (Y - Z * (A * ZtY));


  if (get_val) {
    ll += 2.0 * llt.matrixLLT().diagonal().array().log().sum();
    ll += Y.dot(a) + n * log(2.0 * M_PI * psi_r);
    ll *= -0.5;
  }

  if (get_score) {
    // Stochastic part of the restricted score for psi
    s_psi(rm1) = 0.5 * a.dot(a);
    Eigen::VectorXd v = Z.transpose() * a;
    for (int ii = 0; ii < rm1; ii++) {
      s_psi(ii) = 0.5 * v.dot(H.middleCols(ii * q, q) * v);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  // NOTHING BELOW SHOULD DEPEND ON Y
  /////////////////////////////////////////////////////////////////////////////
  if (get_inf) {
    A = A * ZtZ; // q x q, called M in manuscript
    s_psi(rm1) -= ((0.5 / psi_r) * n - (0.5 / psi_r) * A.diagonal().sum());

    Eigen::SparseMatrix<double> E = A.transpose(); // q x q sparse
    I_psi(rm1, rm1) = (0.5 / (psi_r * psi_r)) * (n - 2.0 * A.diagonal().sum() + E.cwiseProduct(A).sum());

    E = ZtZ * A;
    Eigen::MatrixXd D = XtZ * A; // p x q storage
    Eigen::MatrixXd XtSiZ = (1.0 / psi_r) * (XtZ - D);
    Eigen::MatrixXd XtSi2Z = (1.0 / (psi_r * psi_r)) * (XtZ - 2.0 * D + D * A);

    Eigen::MatrixXd C = B * XtZ.transpose(); // p x p storage, here XtZA ZtX
    Eigen::MatrixXd G = B * ZtZ * B.transpose(); // p x q, here XtZA ZtZ AtZtX
    Eigen::MatrixXd XtSi2X = (1.0 / (psi_r * psi_r)) * (XtX - 2 * C + G);
    Eigen::MatrixXd XtSi3X = (1.0 / (psi_r * psi_r * psi_r)) *
      (XtX - 3 * C + 3 * G - D * (A * B.transpose()));
    C = llt.solve(XtSi2X);

    I_psi(rm1, rm1) += 0.5 * C.transpose().cwiseProduct(C).sum();
    s_psi(rm1) += 0.5 * C.diagonal().sum();
    I_psi(rm1, rm1) -= llt.solve(XtSi3X).diagonal().sum();

    A = (1.0 / (psi_r * psi_r)) * (ZtZ - 2.0 * E + E * A); // ZtSi2Z right now
    E = (1.0 / psi_r) * (ZtZ - E); // Now holds ZtSiZ
    D = llt.solve(XtSiZ);

    // It is possible this loop can be moved into the loop below,
    // but may need additional storage of H
    Eigen::SparseMatrix<double> H_b1(q, q);
    for(int ii = 0; ii < rm1; ii++) {
      H_b1 = H.middleCols(ii * q, q);
      s_psi(ii) -= 0.5 * (E - XtSiZ.transpose() * D).cwiseProduct(H_b1).sum();
    }

    A = A - 2.0 * D.transpose() * XtSi2Z + XtSiZ.transpose() * (C * D);
    // H is overwritten in preparation for information calculation
    Eigen::MatrixXd H2 = XtSiZ.transpose() * (D * H);// Dense
    H = E.transpose() * H;

    Eigen::SparseMatrix<double> H_b2(q, q);
    Eigen::MatrixXd H2_b1(q, q);
    Eigen::MatrixXd H2_b2(q, q);

    for (int ii = 0; ii < rm1; ii++) {
      H_b1 = H.middleCols(ii * q, q);
      H2_b1 = H2.middleCols(ii * q, q);

      I_psi(ii, rm1) = I_psi(rm1, ii) = 0.5 * A.cwiseProduct(H_b1).sum();
      for (int jj = ii; jj < rm1; jj++) {
        H_b2 = H.middleCols(jj * q, q).transpose();
        H2_b2 = H2.middleCols(jj * q, q).transpose();
         //Continue here
        I_psi(ii, jj) = 0.5 * H_b1.cwiseProduct(H_b2).sum() -
          H_b1.cwiseProduct(H2_b2).sum() + 0.5 * H2_b1.cwiseProduct(H2_b2).sum();
        I_psi(jj, ii) = I_psi(ii, jj);
      }
    }
  } else if (get_score) {
    A = A * ZtZ;
    s_psi(rm1) = s_psi(rm1) - (0.5 / psi_r) * n + (0.5 / psi_r) * A.diagonal().sum();
    Eigen::MatrixXd D = XtZ * A;
    Eigen::MatrixXd C = B * XtZ.transpose();
    C = llt.solve((1.0 / (psi_r * psi_r)) * (XtX - 2 * C + B * (ZtZ * B.transpose())));
    s_psi(rm1) += 0.5 * C.diagonal().sum();
    A = ZtZ * A;
    A = (1.0 / psi_r) * (ZtZ - A);
    Eigen::MatrixXd XtSiZ = (1.0 / psi_r) * (XtZ - D);    // p*q
    D = llt.solve(XtSiZ);
    Eigen::MatrixXd V = A - XtSiZ.transpose() * D;
    for (int ii = 0; ii < rm1; ii++) {
      s_psi(ii) -= 0.5 * V.cwiseProduct(H.middleCols(ii * q, q)).sum();
    }
  }
  Eigen::MatrixXd lltU = llt.matrixU();
  return Rcpp::List::create(
    Rcpp::Named("value") = ll,
    Rcpp::Named("score") = s_psi,
    Rcpp::Named("inf_mat") = I_psi,
    Rcpp::Named("beta") = beta_tilde,
    Rcpp::Named("I_b_inv_chol") = lltU);
}

