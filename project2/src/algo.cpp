#include <armadillo>
#include <cmath>
#include "algo.hpp"

arma::vec exactEigval(const int N, const double d, const double a)
{
	arma::vec eigval(N);

	for (int i = 0; i < N; i++)
	{
		eigval(i) = d + 2*a * std::cos((M_PI * (i+1))/(N+1));
	}

	return eigval;
}


arma::mat exactEigvec(const int N, const double d, const double a)
{
	arma::mat eigvec(N,N);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			eigvec(j,i) = std::sin((M_PI*(i+1)*(j+1))/(N+1));
		}
	}

	return arma::normalise(eigvec);
}

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);
