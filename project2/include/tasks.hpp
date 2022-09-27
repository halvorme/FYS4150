#ifndef __tasks_hpp__
#define __tasks_hpp__


#include <armadillo>
#include "utils.hpp"
#include "algo.hpp"

int problem2(const int N, const arma::mat A, const double d, const double a, arma::vec& eigval_anal, arma::mat& eigvec_anal);

int problem3();

int problem4(const arma::mat& A, const double eps, const int maxiter, arma::vec& eigval_jac, 
				arma::mat& eigvec_jac, arma::vec& eigval_anal, arma::mat& eigvec_anal);

int problem5(const double eps, const int maxiter, const int N_max);

int problem6(const double eps, const int maxiter);

#endif