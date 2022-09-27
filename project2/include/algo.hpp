#ifndef __algo_hpp__
#define __algo_hpp__


#include <armadillo>
#include "utils.hpp"

// Calculates the exact eigenvalues and -vectors from the analytic expressions
arma::vec analEigval(const int N, const double d, const double a);
arma::mat analEigvec(const int N, const double d, const double a);


// Performs one Jacobi rotation, on the symmetric matrix 'A'.
// Takes the matrix 'A', rotation matrix 'R' and a set of indices, 'k', 'l'
// Modifies 'A' and 'R' according to the Jacobi method. 
// The entry A(k,l) is the one set to zero.
int jacobi_rotate(arma::mat& A, arma::mat& R, const int k, const int l);


// Implements Jacobi's method to find eigenvalues and -vectors of the matrix 'A'.
// Calls on 'jacobi_rotate()' until all off-diagonal elements are smaller than 'eps'.
// If the algorithm coverges in less than 'maxiter' iterations it writes the result 
// to 'eigenvalues' and 'eigenvectors'. 
// Number of iterations are written to the integer "iterations"
int jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iter, bool& converged);


#endif