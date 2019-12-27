// -----------------------------------------
// Reference Code:
// https://github.com/jaredhuling/rfunctions
// -----------------------------------------


#ifndef _GENINV_H
#define _GENINV_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/LU>
#include <Eigen/SparseCholesky>
#include <vector> 
#include <functional> 
#include <algorithm> 
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Lower;


//computes X'WX where W is diagonal (input w as vector)
MatrixXd xtwx(const MatrixXd& xx, const MatrixXd& ww);

//computes X'SX where S is not diagonal (input ss as matrix)
MatrixXd xtsx(const MatrixXd& xx, const MatrixXd& ss);

//computes X'X 
MatrixXd xtx(const MatrixXd& xx);

//computes XX' 
MatrixXd xxt(const MatrixXd& xx);

//computes det(XX)
double det(const MatrixXd& xx);

//computes log det(XX)
double ldet(const MatrixXd& xx);

//solve Ax = b using conjugate gradient
MatrixXd conjugate_gradient(const MatrixXd& A, const VectorXd& b, int maxit, double tol);

RcppExport SEXP cov2cor(SEXP, double=1.0);

RcppExport SEXP geninv(SEXP);

#endif
