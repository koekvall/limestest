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
                    Eigen::MappedSparseMatrix<double> ZtZXe, //X is not needed
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

  // Eigen::SparseMatrix<double> ZtZ = Z.transpose() * Z;
  // Eigen::VectorXd Zte = Z.transpose() * e;

  // Pre-compute Psi0ZtZ (columns 1:q), Pzi0ZtX (columns (q + 1):(q + p)),
  // and Psi0Zte (column q + p + 1)
  Eigen::SparseMatrix<double> A = Psi0 * ZtZXe;
  // Eigen::VectorXd PZe = Psi0 * Zte;

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
      Eigen::VectorXd Zu = (Z.transpose() * u); // just use u

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

