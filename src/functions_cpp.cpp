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
//' @param psi_mr A vector of covariance parameter (see ?loglikelihood)
//' @param H Sparse matrix of derivatives of Psi with respect to elements of psi,
//'        \eqn{H = [H_1, \dots , H_{r - 1}]}, where \eqn{H_j = \partial \Psi / \partial \psi_j}.
//' @return The covariance matrix \eqn{\Psi}
// [[Rcpp::export]]
Eigen::SparseMatrix<double> Psi_from_H_cpp(const Eigen::Map<Eigen::VectorXd> psi_mr,
                                           const Eigen::MappedSparseMatrix<double> H) { //
  int q = H.rows();
  int rm1 = H.cols() / q;
  Eigen::SparseMatrix<double> Psi(q, q);
  for (int ii = 0; ii < rm1; ii++) {
     if(psi_mr(ii) != 0.0){ // avoid initializing elements that are zero anyway
       Psi += psi_mr(ii) * H.middleCols(ii * q, q);
     }
  }
  return Psi;
}

//' Log-likelihood using RcppEigen
//'
//' Computes the log-likelihood, score vector, and information matrix
//' for the covariance parameter vector in a linear mixed effects model.
//'
//' @param e Vector of length \eqn{n} of errors, or residuals, \eqn{e = Y - X \beta}.
//' @param Z Sparse \eqn{n \times q} random effect design matrix of class \code{dgCMatrix}
//' @param Zte Precomputed vector \code{crossprod(Z, e)} of class \code{numeric}
//' @param XtZ Precomputed matrix \code{crossprod(X, Z)} of class \code{matrix}
//' @param ZtZ Precomputed matrix \code{crossprod(Z)} of class \code{dgCMatrix}
//' @param Psi_r The \eqn{q\times q} covariance matrix of random effects (\eqn{\Psi}) divided by error
//'        variance, \eqn{\Psi_r = \Psi / \psi_r}.
//' @param psi_r The error variance \eqn{\psi_r > 0}.
//' @param H Sparse \eqn{q \times (qr - q)} matrix of horizontally concatenated
//'        derivatives of \eqn{\Psi} (see details) of class \code{dgCMatrix}.
//' @param get_val If \code{TRUE}, the value of the loglikelihood is computed
//' @param get_score If \code{TRUE} the score vector is calculated.
//' @param get_inf If \code{TRUE}, an information matrix is calculated.
//' @param expected If \code{TRUE}, the expected information is calculated; otherwise
//' the observed, or negative Hessian of the loglikelihood.
//'
//' @return A list with components:
//' \item{value}{The value of the log-likelihood}
//' \item{score}{The score, or gradient of the log-likelihood, for \eqn{\psi}}
//' \item{inf_mat}{The information matrix for \eqn{\psi}}
//'
//' @details The model is \deqn{Y = X\beta + Z U + E,} where \eqn{U \sim N_q(0, \Psi)}
//' and \eqn{E \sim N_n(0, \psi_r I_n)}. The first \eqn{r - 1} elements of \eqn{\psi}
//' parameterize \eqn{\Psi}, while the \eqn{r}th and last element is the error
//' variance. It is assumed that \eqn{H_j = \partial \Psi / \partial \psi_j} is
//' a (usually sparse) matrix of zeros and ones, \eqn{j \in \{1, \dots, r - 1\}},
//' and that \eqn{\Psi = \sum_{j = 1}^{r - 1}\psi_j H_j}. Thus, \eqn{\psi_1, \dots, \psi_{r - 1}}
//' are variances and covariances of random effects.
//' The argument matrix \code{H} is \eqn{H = [H_1, \dots, H_{r - 1}]}.
//'
//' The fixed effects  \eqn{\beta} affect the likelihood only through the
//' precomputed \eqn{e = Y - X\beta}.
//'
//' The score and information for \eqn{\beta} are not computed.
//'
//' @useDynLib limestest, .registration=TRUE
//' @import Matrix
// [[Rcpp::export]]

Rcpp::List loglik(
                  const Eigen::MappedSparseMatrix<double> Psi_r,
                  const double psi_r,
                  Eigen::SparseMatrix<double> H,
                  Eigen::VectorXd e,
                  const Eigen::Map<Eigen::MatrixXd> X,
                  const Eigen::MappedSparseMatrix<double> Z,
                  const Eigen::Map<Eigen::MatrixXd> XtX,
                  const Eigen::Map<Eigen::MatrixXd> XtZ,
                  const Eigen::MappedSparseMatrix<double> ZtZ,
                  const bool get_val = true,
                  const bool get_score = true,
                  const bool get_inf = true,
                  const bool expected = true) {
  // Define dimensions
  int p = X.cols();
  int n = e.size();
  int q = Psi_r.cols();
  int r = H.cols() / q + 1;

  // Initialize returns
  double ll = NA_REAL;
  Eigen::VectorXd S(p + r);
  Eigen::MatrixXd I(p + r, p + r);

  // Initialize identity matrix
  Eigen::SparseMatrix<double> Id_q(q, q);
  Id_q = Eigen::MatrixXd::Identity(q, q).sparseView();

  // Compute matrices B and A in manuscript
  Eigen::SparseMatrix<double> B = Psi_r * ZtZ + Id_q;
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(B);
  bool stop_early = solver.info() != Eigen::Success;

  // Add loglik term before overwriting
  if (get_val) {
    ll = -0.5 * solver.logAbsDeterminant() - 0.5 * (double)n * log(psi_r);
    stop_early = stop_early || (solver.info() != Eigen::Success);
  }

  // This is A = (I_q + Psi_r ZtZ)^{-1}Psi_r
  Eigen::SparseMatrix<double> A = solver.solve(Psi_r);

  stop_early = stop_early || (solver.info() != Eigen::Success);

  if(stop_early || (psi_r <= 0.0)){
    ll = -R_PosInf;
    return Rcpp::List::create(Rcpp::Named("value") = ll,
                              Rcpp::Named("score") = S,
                              Rcpp::Named("inf_mat") = I);
  }


  Eigen::VectorXd etilde = (1.0 / psi_r) * (e - Z * (A * (Z.transpose() * e)));

  if (get_val) {
    ll = ll - 0.5 * etilde.dot(e);
  }

  // Get score for beta
  S.head(p) = X.transpose() * etilde;

  // Get score for psi
  Eigen::SparseMatrix<double> C = Id_q - A * ZtZ;

  S(p + r - 1) = 0.5 * etilde.dot(etilde) - (1 / psi_r) * (n - q +
    C.diagonal().sum());

  // v is Z'Sigma^{-1}e, w is e'Sigma^{-1}Z[H1...Hr]
  Eigen::VectorXd v = (Z.transpose() * etilde);
  Eigen::VectorXd w = H.transpose() * v;
  for (int ii = 0; ii < r - 1; ii++) {
    S(p + ii) = 0.5 * w.middleRows(ii * q, q).dot(v);
  }

  // Overwriting B with Z' \Sigma^{-1} Z
  B = (1 / psi_r) * ZtZ * C;

  if (!get_inf) {
    for (int ii = 0; ii < r - 1; ii++) {
      S(p + ii) += (0.5 / psi_r) * H.middleCols(ii * q, q).cwiseProduct(B).sum();
    }
  } else {
    I(p + r - 1, p + r - 1) = (0.5 / (psi_r * psi_r)) *
      (n - q + C.cwiseProduct(C.transpose()).sum());

    // Replace H by Z'\Sigma^{-1}Z H
    H = B * H;
    Eigen::SparseMatrix<double> H1(q, q);
    Eigen::SparseMatrix<double> H2(q, q);

    for (int jj = 0; jj < r - 1; jj++) {
      H1 = H.middleCols(jj * q, q);
      I(p + jj, p + r - 1) = (0.5 / psi_r) * H1.cwiseProduct(C).sum();
      S(p + jj) += (0.5 / psi_r) * H1.diagonal().sum();
      for (int ii = 0; ii <= jj; ii++) {
        H2 = H.middleCols(ii * q, q).transpose();
        I(p + ii, p + jj) = (0.5 / (psi_r * psi_r)) * H1.cwiseProduct(H2).sum();
      }
    }
    if (!expected) {
      I.bottomRightCorner(q, q) = -I.bottomRightCorner(q, q);

      // Replace e by Sigma^{-2}e = \Sigma^{-1}\tilde{e} = \check{e}
      e = (1 / psi_r) * (etilde - Z * (A * v));
      Eigen::VectorXd u = Z.transpose() * e;

      I.topRightCorner(p, 1) = X.transpose() * e;
      I(p + r - 1, p + r - 1) += e.dot(e);

      for (int jj = 0; jj < r - 1; jj++) {
        I.block(0, p + jj, p, 1) = (1 / psi_r) * (XtZ *(C * w.middleRows(jj * q, q)));
        I(p + jj, r - 1) += w.middleRows(jj * q, q).dot(u);
        for (int ii = 0; ii <= jj; ii++) {
          I(p + ii, p + jj) +=  (B * w.middleRows(ii * q, q)).dot(w.middleRows(jj * q, q));
        }
      }
    }
  }
  I = I.selfadjointView<Eigen::Upper>();
  return Rcpp::List::create(Rcpp::Named("value") = ll,
                            Rcpp::Named("score") = S,
                            Rcpp::Named("inf_mat") = I);
}

//' Restricted log-likelihood using RcppEigen
//'
//' Computes the restricted log-likelihood, score vector, and information matrix
//' for the covariance parameter vector in a linear mixed effects model.
//'
//' @param Y Vector of length \eqn{n} of responses, of class \code{numeric}.
//' @param X Matrix of size \eqn{n \times p} of predictors, of class \code{matrix}
//' @param Z Sparse \eqn{n \times q} random effect design matrix of class \code{dgCMatrix}
//' @param XtY Precomputed vector \code{crossprod(X, Y)} of class \code{numeric}
//' @param ZtY Precomputed vector \code{crossprod(Z, Y)} of class \code{numeric}
//' @param XtX Precomputed matrix \code{crossprod(X)} of class \code{matrix}
//' @param XtZ Precomputed matrix \code{crossprod(X, Z)} of class \code{matrix}
//' @param ZtZ Precomputed matrix \code{crossprod(Z)} of class \code{dgCMatrix}
//' @param Psi_r The \eqn{q\times q} covariance matrix of random effects (\eqn{\Psi}) divided by error
//'        variance, \eqn{\Psi_r = \Psi / \psi_r}.
//' @param psi_r The error variance \eqn{\psi_r > 0}.
//' @param H Sparse \eqn{q \times (qr - q)} matrix of horizontally concatenated
//'        derivatives of \eqn{\Psi} (see details) of class \code{dgCMatrix}.
//' @param get_val If \code{TRUE}, the value of the loglikelihood is computed
//' @param get_score If \code{TRUE} the score vector is calculated.
//' @param get_inf If \code{TRUE}, an information matrix is calculated.
//' the observed, or negative Hessian of the loglikelihood.
//'
//' @return A list with components:
//' \item{value}{The value of the restricted log-likelihood}
//' \item{score}{The restricted score, or gradient of the restriced log-likelihood, for \eqn{\psi}}
//' \item{inf_mat}{The restricted information matrix for \eqn{\psi}}
//' \item{beta}{Partial maximizer of the regular likelihood in \eqn{\beta},
//'   \eqn{\tilde{\beta} = (X' \Sigma^{-1} X)' X' \Sigma^{-1}Y},
//'   where \eqn{\Sigma = Z\Psi Z' + \psi_r I_n}}
//' \item{I_b_inv_chol}{Cholesky root of the expected inverse information matrix
//'   for \eqn{\beta}, \eqn{I(\beta; \psi) = X' \Sigma^{-1} X}}
//'
//' @details The model is \deqn{Y = X\beta + Z U + E,} where \eqn{U \sim N_q(0, \Psi)}
//' and \eqn{E \sim N_n(0, \psi_r I_n)}. The first \eqn{r - 1} elements of \eqn{\psi}
//' parameterize \eqn{\Psi}, while the \eqn{r}th and last element is the error
//' variance. It is assumed that \eqn{H_j = \partial \Psi / \partial \psi_j} is
//' a (usually sparse) matrix of zeros and ones, \eqn{j \in \{1, \dots, r - 1\}},
//' and that \eqn{\Psi = \sum_{j = 1}^{r - 1}\psi_j H_j}. Thus, \eqn{\psi_1, \dots, \psi_{r - 1}}
//' are variances and covariances of random effects.
//' The argument matrix \code{H} is \eqn{H = [H_1, \dots, H_{r - 1}]}.
//'
//'
// [[Rcpp::export]]
Rcpp::List res_ll_cpp(Eigen::VectorXd Y,
                      const Eigen::Map<Eigen::MatrixXd> X,
                      const Eigen::MappedSparseMatrix<double> Z,
                      const Eigen::Map<Eigen::MatrixXd> XtY,
                      const Eigen::Map<Eigen::MatrixXd> ZtY,
                      const Eigen::Map<Eigen::MatrixXd> XtX,
                      const Eigen::Map<Eigen::MatrixXd> XtZ,
                      const Eigen::MappedSparseMatrix<double> ZtZ,
                      const Eigen::MappedSparseMatrix<double> Psi_r,
                      const double psi_r,
                      Eigen::SparseMatrix<double> H,
                      const bool get_val = true,
                      const bool get_score = true,
                      const bool get_inf = true)
{
  // Define dimensions
  int n = Y.size();
  int q = Psi_r.cols();
  int rm1 = H.cols() / q;
  int r = rm1 + 1;
  int p = X.cols();
  // loglikelihood to return
  double ll = NA_REAL;
  // Score vector to return
  Eigen::VectorXd s_psi = Eigen::VectorXd::Zero(r);
  // Information matrix to return
  Eigen::MatrixXd I_psi =  Eigen::MatrixXd::Zero(r, r);


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
      Rcpp::Named("value") = -R_PosInf,
      Rcpp::Named("score") = s_psi,
      Rcpp::Named("inf_mat") = I_psi,
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
  Eigen::VectorXd a = (1.0 / psi_r) * (Y - Z * (A * (Z.transpose() * Y)));


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
    A = A - 2.0 * D.transpose() * XtSi2Z + XtSiZ.transpose() * (C * D);

    // It is possible this loop can be moved into the loop below,
    // but may need additional storage of H
    Eigen::SparseMatrix<double> H_b1(q, q);
    for(int ii = 0; ii < rm1; ii++) {
      H_b1 = H.middleCols(ii * q, q);
      s_psi(ii) -= 0.5 * (E - XtSiZ.transpose() * D).cwiseProduct(H_b1).sum();
      I_psi(ii, rm1) = 0.5 * A.cwiseProduct(H_b1).sum();
    }

    // H is overwritten in preparation for information calculation
    Eigen::MatrixXd H2 = XtSiZ.transpose() * (D * H);// Dense
    H = E.transpose() * H;

    Eigen::SparseMatrix<double> H_b2(q, q);
    Eigen::MatrixXd H2_b1(q, q);
    Eigen::MatrixXd H2_b2(q, q);
    for (int jj = 0; jj < rm1; jj++) {
      H_b1 = H.middleCols(jj * q, q);
      H2_b1 = H2.middleCols(jj * q, q);
      for (int ii = 0; ii <= jj; ii++) {
        H_b2 = H.middleCols(ii * q, q).transpose();
        H2_b2 = H2.middleCols(ii * q, q).transpose();
         //Continue here
        I_psi(ii, jj) = 0.5 * H_b1.cwiseProduct(H_b2).sum() -
          H_b1.cwiseProduct(H2_b2).sum() + 0.5 * H2_b1.cwiseProduct(H2_b2).sum();
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

  U = llt.matrixU();
  I_psi = I_psi.selfadjointView<Eigen::Upper>();
  return Rcpp::List::create(
    Rcpp::Named("value") = ll,
    Rcpp::Named("score") = s_psi,
    Rcpp::Named("inf_mat") = I_psi,
    Rcpp::Named("beta") = beta_tilde,
    Rcpp::Named("I_b_inv_chol") = U);
}

