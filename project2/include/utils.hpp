#ifndef __utils_hpp__
#define __utils_hpp__


#include <armadillo>

// Finds the biggest off-diagonal element of the symmetric matrix A.
// - The matrix indices of the max element are returned by writing to the  
//   int references k and l (row and column, respectively)
// - The value of the max element A(k,l) is returned as the function
//   return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);


// Initialises a tridiagonal matrix of size n*n, with a single value on each of the three diagonals. 
// The entries on the subdiagonal get the value 'l', on the main diagonal 'd' and on the superdiagonal 'u'.
arma::mat triMat(const int n, const double l, const double d, const double u);

arma::mat triMat(const arma::vec* l, const arma::vec* d, const arma::vec* u);


#endif 