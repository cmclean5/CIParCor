// -----------------------------------------
// Reference Code:
// https://github.com/jaredhuling/rfunctions
// -----------------------------------------

#include "geninv.h"

using namespace Rcpp;
using namespace RcppEigen;

//computes X'WX where W is diagonal (input w as vector)
MatrixXd xtwx(const MatrixXd& xx, const MatrixXd& ww) {
  const int n(xx.cols());
  MatrixXd AtWA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint() * ww.asDiagonal()));
  return (AtWA);
}

//computes X'SX where S is not diagonal (input ss as matrix)
MatrixXd xtsx(const MatrixXd& xx, const MatrixXd& ss) {
  const int n(xx.cols());
  MatrixXd AtSA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint() * ss));
  return (AtSA);
}

MatrixXd xtx(const MatrixXd& xx) {
  const int n(xx.cols());
  MatrixXd AtA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint()));
  return (AtA);
}

MatrixXd xxt(const MatrixXd& xx) {
  const int m(xx.rows());
  MatrixXd AtA(MatrixXd(m, m).setZero().
    selfadjointView<Lower>().rankUpdate(xx));
  return (AtA);
}

double det(const MatrixXd& xx) {
  return (xx.determinant());
}

double ldet(const MatrixXd& xx) {
  const VectorXd Dvec(xx.ldlt().vectorD());
  return(Dvec.array().log().sum());
}


MatrixXd conjugate_gradient(const MatrixXd& A, const VectorXd& b, int maxit, double tol)
{

  const int n(A.cols());
  VectorXd x(n);
  VectorXd r(n);
  VectorXd p(n);
  VectorXd Ap(n);
  x.fill(0);
  
  double rsold;
  double rsnew;
  double alpha;
  int iters = maxit; 
  
  r = b; 
  p = r;
  rsold = r.squaredNorm();
  
  for (int i = 0; i < maxit; i++) {
    Ap = A * p;
    alpha = rsold / (p.transpose() * Ap);
    x = x + (alpha * p.array()).matrix();
    r = r - (alpha * Ap.array()).matrix();
    rsnew = r.squaredNorm();
    if (sqrt(rsnew) < tol) {
      iters = i;
      break;
    }
    p = r + ((rsnew / rsold) * p.array()).matrix();
    rsold = rsnew;
  }
  return(x);
}

//------
// Note short guide to RcppEigen:
// https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
// Convert from Covariance matrix to Correlation matrix...
// this is not very RcppEigen'sih yet. 
//------
RcppExport SEXP cov2cor(SEXP sCov, double sign){

  using namespace Rcpp;
  using namespace RcppEigen;
  try{
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    
    Eigen::Map<MatrixXd> Cov(as<Map<MatrixXd> >(sCov)); 
    const int n(Cov.rows());
    const int m(Cov.cols());
    const int mn(std::min(n, m));    

    VectorXd sigma(VectorXd(mn).setZero());
    MatrixXd Cor(MatrixXd(mn, mn).setZero());

    for(int i=0; i<mn; i++){
      sigma(i) = Cov.array()((i*mn)+i);
    }
    
    for(int i=0; i < sigma.size(); i++){
      sigma(i) = std::sqrt(sigma(i));
    }
    
    for(int i=0; i<mn; i++){
      for(int j=0; j<mn; j++){
	Cor.array()((i*mn)+j) = sign * Cov.array()((i*mn)+j) / sigma(j) / sigma(i);
      }
    }    

    return wrap(Cor);

  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
    
}

// -----------------------------------------
// Moore-Penrose Generalized Inverse
// Reference:
// https://arxiv.org/ftp/arxiv/papers/0804/0804.4809.pdf
// -----------------------------------------
RcppExport SEXP geninv(SEXP GG)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::Lower;
    const Eigen::Map<MatrixXd> G(as<Map<MatrixXd> >(GG));
    const int n(G.rows());
    const int m(G.cols());
    const int mn(std::min(n, m));
    
    bool transp(false);
    double tol(1.0e-10);
    MatrixXd A(MatrixXd(mn, mn));
    MatrixXd L(MatrixXd(mn, mn).setZero());
    
    
    
    if (n < m) {
      transp = true;
      A = xxt(G);
    } else {
      A = xtx(G);
    }

    int r = 0;
    for (int k = 0; k < mn; k++) {
      r++;
      
      if (r == 1) {
        L.block(k, r - 1, mn - k, 1) = A.block(k, k, mn - k, 1);
      } else {
        L.block(k, r - 1, mn - k, 1) = A.block(k, k, mn - k, 1) - 
                L.block(k, 0, mn - k, r - 1) * L.block(k, 0, 1, r - 1).adjoint();
      }
      
      if (L(k, r - 1) > tol) {
        L.block(k, r - 1, 1, 1) = L.block(k, r - 1, 1, 1).array().sqrt();
        if (k + 1 < mn) {
          L.block(k + 1, r - 1, mn - k - 1, 1) = L.block(k + 1, r - 1, mn - k - 1, 1) / L(k, r - 1);
        }
      } else {
        r--;
      }
    }

    MatrixXd M(MatrixXd(r, r));
    M = xtx(L.block(0, 0, mn, r)).inverse();

    MatrixXd Y(MatrixXd(m, n));
    
    if (transp) {
      Y = G.adjoint() * L.block(0, 0, mn, r) * M * M * L.block(0, 0, mn, r).adjoint();
    } else {
      Y = L.block(0, 0, mn, r) * M * M * L.block(0, 0, mn, r).adjoint() * G.adjoint();
    }

    return wrap(Y);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall

}

