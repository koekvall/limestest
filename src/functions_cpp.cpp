#include <RcppEigen.h>
#include <Eigen/SparseCholesky>


// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<double> res_ll_cpp(const Eigen::Map<Eigen::MatrixXd>& XtX,
                  const Eigen::Map<Eigen::MatrixXd>& XtY,
                  const Eigen::Map<Eigen::MatrixXd>& XtZ,
                  const Eigen::MappedSparseMatrix<double>& ZtZ,
                  const Eigen::Map<Eigen::MatrixXd>& YtZ,
                  Eigen::Map<Eigen::VectorXd> Y,
                  const Eigen::Map<Eigen::MatrixXd>& X,
                  const Eigen::MappedSparseMatrix<double>& Z,
                  const Eigen::MappedSparseMatrix<double>& H,
                  const Eigen::MappedSparseMatrix<double>& Psi0,
                  const Eigen::Map<Eigen::VectorXd> psi0,
                  const bool& lik = true,
                  const bool& score = false,
                  const bool& finf = false)
{
  // Define dimensions
  const int n = Y.size();
  const int q = Psi0.cols();
  const int r = H.cols() / q;
  const int p = XtX.cols();

  // Loglikelihood to return
  double ll = NA_REAL;

  // Score vector to return
  Eigen::VectorXd s_psi(r + 1);
  s_psi.fill(NA_REAL);

  // Fisher information to return
  Eigen::MatrixXd I_psi(r + 1, r + 1);
  I_psi.fill(NA_REAL);

  // Pre-compute A = (I_q + Psi0 Z'Z)^{-1} Psi0
  Eigen::SparseMatrix<double> A = Psi0 * ZtZ + Eigen::MatrixXd::Identity(q, q);
  // Eigen::LLT<Eigen::SparseMatrix<double>> A_chol(A);
  return A;
}
