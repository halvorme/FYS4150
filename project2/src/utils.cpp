#include <armadillo>
#include "utils.hpp"

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l)
{
	int N = A.n_rows;

	double maxval = 0.;

	for (int i = 0; i < N; i++)
	{
		for (int j = i+1; j < N; j++)
		{
			double testval = std::abs(A(i,j));
			if (testval > maxval)
			{
				maxval = testval;
				k = i;
				l = j;
			}
		}
	}

	return maxval;
}


arma::mat triMat(const int n, const double l, const double d, const double u)
{
	arma::mat A = d * arma::mat(n, n, arma::fill::eye);

	for (int i = 0; i < n-1; i++)
	{
		A(i, i+1) = u;
		A(i+1, i) = l;
	}

	return A;
}


arma::mat triMat(const arma::vec* l, const arma::vec* d, const arma::vec* u)
{
	int n = arma::size(*d)[0];
	arma::mat A = arma::mat(n,n, arma::fill::zeros);

	for (int i = 0; i < n-1; i++)
	{
		A(i, i) = (*d)(i);
		A(i, i+1) = (*u)(i);
		A(i+1, i) = (*l)(i);
	}
	A(n-1, n-1) = (*d)(n-1);

	return A;
}

