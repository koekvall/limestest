#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]


// Eigen::SparseMatrix<double> sparseCrossprod(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B) {
//   return A.transpose()*B;
// }


// [[Rcpp::export]]
Eigen::SparseMatrix<double> get_PsiRcpp(Eigen::Map<Eigen::VectorXd> psi, Eigen::MappedSparseMatrix<double> H) { //
  int q = H.rows();
  int r = H.cols() / q;
  Eigen::SparseMatrix<double> Psi(10,10);
  for (int ii = 0; ii < r; ii++) {
    Psi += psi(ii) * H.middleCols(ii * q, (ii+1) * q);
  }
  return Psi;
}


// [[Rcpp::export]]
List loglik_psiRcpp(Eigen::MappedSparseMatrix<double> Z,
                                           Eigen::MappedSparseMatrix<double> ZtZXe,
                               Eigen::Map<Eigen::VectorXd> e,
                               Eigen::MappedSparseMatrix<double> H,
                               Eigen::MappedSparseMatrix<double> Psi0,
                               double psi0,
                               bool loglik = true,
                               bool score = true,
                               bool finf = true,
                               bool expected = true) {

  // Define dimensions
  int n = e.size();
  int q = Psi0.cols();
  int r = H.cols() / q;
  // loglikelihood to return
  double ll = NA_REAL;
  // Score vector to return
  Eigen::VectorXd s_psi(r + 1);
  // Fisher information to return
  Eigen::MatrixXd I_psi(r + 1, r + 1);

  // Pre-compute Psi0ZtZ (columns 1:q), Pzi0ZtX (columns (q + 1):(q + p)),
  // and Psi0Zte (column q + p + 1)
  Eigen::SparseMatrix<double> A = Psi0 * ZtZXe;

  // solver for Psi0ZtZ + I_q
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A.leftCols(q) + Eigen::MatrixXd::Identity(q,q).sparseView());

  // Add loglik term before overwriting
  if (loglik) {
    ll = -0.5 * solver.logAbsDeterminant() - 0.5 * n * log(psi0);
  }

  // Matrix denoted M in manuscript is A[, 1:q]
  A = solver.solve(A).eval();

  Eigen::VectorXd e_save(n) ;
  if(loglik | (finf & !expected)){
    e_save = e;
  }
  e = (1.0 / psi0) * (e - Z * A.rightCols(1)).eval(); // = Sigma^{-1}e
  if (loglik) {
    ll = ll - 0.5 * e.dot(e_save);
  }

  // find trace of M
  double trace_M = 0.;
  for (int i = 0; i < q; i++) {
    trace_M += A.coeffRef(i,i);
  }

  s_psi(0) = 0.5 * e.dot(e) - (0.5 / psi0) * (n - trace_M);

  // v is Z'Sigma^{-1}e, w is e'Sigma^{-1}Z[H1...Hr]
  Eigen::SparseMatrix<double> v = (Z.transpose() * e).sparseView(); // v is not sparse in the example
  Eigen::SparseMatrix<double> w = v.transpose() * H;

  for (int i = 1; i <= r; i++) {
    Eigen::SparseMatrix<double> temp = w.middleCols((i-1)*q,q) * v;
    s_psi(i) = 0.5 * temp.coeffRef(0,0);
  }

  // B = Z'Z (M - I_q) in paper notation, sparse
  Eigen::SparseMatrix<double> B = A.leftCols(q);
  B.diagonal() -= Eigen::VectorXd::Constant(q,1);
  B = ZtZXe.leftCols(q) * B;

  if (!finf) {
    // Compute -[ZtZ (I_q - M) * H_1, ..., ZtZ (I_q - M) * H_r] using recycling
    // The "if" is because the calculation is a byproduct of a more expensive one
    // (B %*% H) done to get Fisher information
    for (int i = 1; i <= r; i++) {
      s_psi(i) += (0.5 / psi0) * H.middleCols((i-1)*q,q).cwiseProduct(B).sum();
    }
  } else {
    Eigen::SparseMatrix<double> BH = B * H;
    I_psi(0, 0) = (0.5 / (psi0 * psi0)) * (n - 2 * trace_M +
      Eigen::SparseMatrix<double>(A.leftCols(q).transpose()).cwiseProduct(A.leftCols(q)).sum());

    // M = M-I
    for (int i = 0; i < q; i++) {
      A.coeffRef(i,i) -= 1;
    }

    for (int i = 1; i <= r; i++) {
      I_psi(0, i) = I_psi(i, 0) = (0.5 / (psi0*psi0)) * A.leftCols(q).cwiseProduct(BH.middleCols((i-1)*q,q)).sum();
      s_psi(i) += (0.5 / psi0) * H.middleCols((i-1)*q,q).cwiseProduct(B).sum();
      for (int j = i; j <= r; j++) {
        I_psi(i, j) = I_psi(j, i) = (0.5 / (psi0*psi0)) *
          Eigen::SparseMatrix<double>(BH.middleCols((i-1)*q,q).transpose()).cwiseProduct(BH.middleCols((j-1)*q,q)).sum();
      }
    }
    if (!expected) {
      I_psi = -I_psi.eval();

      // u = Sigma^{-2}e. Some calculations could be saved from before
      Eigen::VectorXd u = (1.0 / (psi0 * psi0)) * (e_save + Z * (-2 * A.rightCols(1) + (A.leftCols(q) + Eigen::MatrixXd::Identity(q,q).sparseView()) * A.rightCols(1)));
      I_psi(0, 0) += e.dot(u);
      Eigen::VectorXd Zu = (Z.transpose() * u);

      for (int i = 1; i <= r; i++) {
        Eigen::MatrixXd temp = w.middleCols((i-1)*q,q) * Zu;
        I_psi(0, i) += temp(0,0);
        I_psi(i, 0) += temp(0,0);
        for (int j = i; j <= r; j++) {
          I_psi(i, j) -= (1 / psi0) * (w.middleCols((i-1)*q,q) * B).cwiseProduct(w.middleCols((j-1)*q,q)).sum();
          if (i != j) {
            I_psi(j, i) -= (1 / psi0) * (w.middleCols((i-1)*q,q) * B).cwiseProduct(w.middleCols((j-1)*q,q)).sum();
          }
        }
      }
    }
  }
  return List::create(Named("ll") = ll, Named("score") = s_psi, Named("finf") = I_psi);
}



// [[Rcpp::export]]
List res_llRcpp(Eigen::Map<Eigen::MatrixXd> X,
                Eigen::Map<Eigen::VectorXd> Y,
                Eigen::MappedSparseMatrix<double> Z,
                Eigen::MappedSparseMatrix<double> H,
                Eigen::MappedSparseMatrix<double> Psi0,
                double psi0,
                bool lik = true,
                bool score = false,
                bool finf = false) {
  // Define dimensions
  int n = Y.size();
  int q = Psi0.cols();
  int r = H.cols() / q;
  int p = X.cols();

  // loglikelihood to return
  double ll = NA_REAL;
  // Score vector to return
  Eigen::VectorXd s_psi(r + 1);
  // Fisher information to return
  Eigen::MatrixXd I_psi(r + 1, r + 1);

  Eigen::SparseMatrix<double> ZtZ = Z.transpose() * Z;
  Eigen::MatrixXd XtZ = X.transpose() * Z;
  Eigen::MatrixXd XtX = X.transpose() * X;

  // Pre-compute A = (I_q + Psi0 Z'Z)^{-1} Psi0
  Eigen::SparseMatrix<double> A = Psi0 * ZtZ + Eigen::MatrixXd::Identity(q,q).sparseView();
  // solver for Psi0ZtZ + I_q
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);

  if (lik) {
    ll = solver.logAbsDeterminant();
  }

  A = solver.solve(Psi0);
  Eigen::MatrixXd B = XtZ * A;

  //Create XtSiX
  Eigen::MatrixXd U = (1.0 / psi0) * (XtX - B * XtZ.transpose());
  // force U to be symmetric
  for (int i=0; i<p-1; i++) {
    for (int j=i+1; j<p; i++) {
      U(j,i) = U(i,j);
    }
  }
  Eigen::LLT<Eigen::MatrixXd> llt(U);

//   if (llt.info() != Eigen::Success) {
//     return List::create(
//       Named("ll") = -R_PosInf,
//       Named("score") = s_psi,
//       Named("finf") = I_psi,
//       Named("beta") = Eigen::VectorXd::Constant(p, NA_REAL),
//       Named("I_b_inv_chol") = Eigen::MatrixXd::Constant(p, p, NA_REAL));
//   }

  // Create XtSiY for use in beta_tilde
  Eigen::VectorXd beta_tilde = (1.0 / psi0) * (X.transpose()*Y - XtZ * A * Z.transpose() * Y);
  beta_tilde = llt.solve(beta_tilde);

  // Replace Y by residuals
  Y = (Y - X * beta_tilde).eval();

  // n x 1 vector for storing \Sigma^{-1}e
  Eigen::VectorXd a = (1.0 / psi0) * (Y - Z * (A * (Z.transpose() * Y)));

  if (lik) {
    ll += 2 * log(llt.matrixL().determinant());//.diagonal().array().log().sum();
    ll += Y.dot(a) + n * log(2 * M_PI * psi0);
    ll *= -0.5;
 }

  if (score) {
    // Stochastic part of the restricted score for psi
    s_psi(0) = 0.5 * a.dot(a);
    Eigen::VectorXd v = Z.transpose() * a;
    for (int i = 1; i <= r; i++) {
      s_psi(i) = 0.5 * v.dot(H.middleCols((i-1)*q,q) * v);
    }
  }

  if (finf && score) {
    A = A * ZtZ; // q x q, called M in manuscript

    s_psi(0) -= 0.5 / psi0 * n + 0.5 / psi0 * A.diagonal().sum();

    I_psi(0, 0) = 0.5 / (psi0 * psi0) * (n - 2 * A.diagonal().sum() + Eigen::SparseMatrix<double>(A.transpose()).cwiseProduct(A).sum());

    Eigen::SparseMatrix<double> E = ZtZ * A; // q x q storage
    Eigen::MatrixXd D = XtZ * A; // p x q storage
    Eigen::MatrixXd XtSiZ = (1.0 / psi0) * (XtZ - D);
    Eigen::MatrixXd XtSi2Z = (1.0 / (psi0 * psi0)) * (XtZ - 2 * D + D * A);

    Eigen::MatrixXd C = B * XtZ.transpose(); // p x p storage, here XtZA ZtX
    Eigen::MatrixXd G = B * ZtZ * B.transpose(); // p x q, here XtZA ZtZ AtZtX
    Eigen::MatrixXd XtSi2X = (1.0 / (psi0 * psi0)) * (XtX - 2 * C + G);
    Eigen::MatrixXd XtSi3X = (1.0 / (psi0 * psi0 * psi0)) * (XtX - 3 * C + 3 * G - D * (A * B.transpose()));
    C = llt.solve(XtSi2X);

    I_psi(0, 0) += 0.5 * C.transpose().cwiseProduct(C).sum();
    s_psi(0) += 0.5 * C.diagonal().sum();

    I_psi(0, 0) -= llt.solve(XtSi3X).diagonal().sum();

    A = (1.0 / (psi0 * psi0)) * (ZtZ - 2 * E + E * A); // ZtSi2Z right now
    E = (1.0 / psi0) * (ZtZ - E); // Now holds ZtSiZ
    D = llt.solve(XtSiZ);
    A = A - 2 * D.transpose() * XtSi2Z + XtSiZ.transpose() * (C * D);
    Eigen::MatrixXd H2 = XtSiZ.transpose() * (D * H);

    for (int i = 1; i <= r; i++) {
      I_psi(i, 0) = I_psi(0, i) = 0.5 * A.cwiseProduct(H.middleCols((i-1)*q,q)).sum();
      s_psi(i) -= 0.5 * (E - XtSiZ.transpose() * D).cwiseProduct(H.middleCols((i-1)*q,q)).sum();
      for (int j = i; j <= r; j++) {
        I_psi(i, j) = I_psi(j, i) = 0.5 * Eigen::SparseMatrix<double>(H.middleCols((i-1)*q,q).transpose()).cwiseProduct(H.middleCols((j-1)*q,q)).sum() -
          H2.middleCols((j-1)*q,q).transpose().cwiseProduct(H.middleCols((i-1)*q,q)).sum() +
          H2.middleCols((i-1)*q,q).transpose().cwiseProduct(H2.middleCols((j-1)*q,q)).sum();
      }
    }
  } else if (score) {
    A = A * ZtZ;
    s_psi(0) -= (0.5 / psi0) * n + (0.5 / psi0) * A.diagonal().sum();
    Eigen::MatrixXd D = XtZ * A;
    Eigen::MatrixXd C = B * XtZ.transpose();
    C = llt.solve((1.0 / (psi0 * psi0)) * (XtX - 2 * C + B * (ZtZ * B.transpose())));
    s_psi(0) += 0.5 * C.diagonal().sum();
    A = ZtZ * A;
    A = (1/ psi0) * (ZtZ - A).eval();
    Eigen::MatrixXd XtSiZ = (1 / psi0) * (XtZ - D);    // p*q
    D = llt.solve(XtSiZ);
    Eigen::MatrixXd V = A - XtSiZ.transpose() * D;
    for (int i = 1; i <= r; i++) {
      s_psi(i) -= 0.5 * V.cwiseProduct(H.middleCols((i-1)*q,q)).sum();
    }
  }
  return List::create(
    Named("ll") = ll, Named("score") = s_psi, Named("finf") = I_psi,
    Named("beta") = beta_tilde);//,
    //Named("I_b_inv_chol") = llt);
}








// // [[Rcpp::export]]
// Eigen::VectorXd uni_test_statRcpp(Eigen::Map<Eigen::VectorXd> test_seq, int test_idx, Eigen::VectorXd psi, double psi0, Eigen::MatrixXd Z, Eigen::MatrixXd ZtZXe, Eigen::VectorXd e, Eigen::MatrixXd H) {
//   int m = test_seq.size();
//   Eigen::VectorXd test_stat(m);
//
//   for (int ii = 0; ii < m; ii++) {
//     if (ii > 0) {
//       psi(test_idx) = test_seq(ii);
//       Eigen::VectorXd Psi = get_PsiRcpp(psi, H);
//       Eigen::VectorXd s = score_psi(Z, ZtZXe, e, H, Psi / psi0, psi0, false);
//       Eigen::VectorXd updated_psi = psi;
//       updated_psi.segment(0, test_idx) = psi.segment(0, test_idx) - s.segment(0, test_idx);
//       updated_psi.segment(test_idx + 1, psi.size() - test_idx - 1) = psi.segment(test_idx + 1, psi.size() - test_idx - 1) - s.segment(test_idx + 1, psi.size() - test_idx - 1);
//       Eigen::MatrixXd finf = score_psi(Z, ZtZXe, e, H, Psi / psi0, psi0, true);
//       Eigen::VectorXd s_new = s - finf * s;
//       updated_psi = updated_psi + finf.ldlt().solve(s_new);
//       psi = updated_psi;
//     }
//
//     Eigen::VectorXd Psi = getPsi(psi, H);
//     Eigen::MatrixXd score_inf = score_psi(Z, ZtZXe, e, H, Psi / psi0, psi0, true);
//     double eff_inf = score_inf(test_idx, test_idx) - (score_inf.block(0, 0, psi.size(), psi.size()).ldlt().solve(score_inf.block(0, test_idx, psi.size(), 1)) * score_inf.block(0, test_idx, psi.size(), 1)).sum();
//     test_stat(ii) = std::pow(score_inf(test_idx, 0), 2) / eff_inf;
//   }
//
//   return test_stat;
// }
