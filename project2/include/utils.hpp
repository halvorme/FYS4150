#ifndef __utils_hpp__
#define __utils_hpp__


#include <armadillo>

// Finds the biggest off-diagonal element of the symmetric matrix A.
// Modifies 'k' and 'l' to the indices of the biggest element
// Returns the absolute value of the biggest element
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);


// Initialises a tridiagonal matrix of size n*n, with a single value on each of the three diagonals. 
// The entries on the subdiagonal get the value 'l', on the main diagonal 'd' and on the superdiagonal 'u'.
arma::mat triMat(const int n, const double l, const double d, const double u);

// Initialises a tridiagonal matrix with 'l' on the subdiagonal, 'd' on the diagonal and 'u' on the superdiagonal
arma::mat triMat(const arma::vec* l, const arma::vec* d, const arma::vec* u);


// Writes the three first eigenvectors to file
int printVec(std::string filename, arma::mat& eigvec, double h);


#endif 