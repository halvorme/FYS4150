#include <armadillo>
#include <iomanip>
#include "utils.hpp"

// Finds the biggest off-diagonal element of the symmetric matrix A.
// Modifies 'k' and 'l' to the indices of the biggest element
// Returns the absolute value of the biggest element
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l)
{
	int N = A.n_rows;

	double maxval = 0.;

	// Runs over upper triangle
	for (int i = 0; i < N; i++)
	{
		for (int j = i+1; j < N; j++)
		{
			double testval = std::abs(A(i,j));

			// Check if biggest
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


// Initialises a tridiagonal matrix of size n*n, with a single value on each of the three diagonals. 
// The entries on the subdiagonal get the value 'l', on the main diagonal 'd' and on the superdiagonal 'u'.
arma::mat triMat(const int n, const double l, const double d, const double u)
{
	// Initialises A with main diagonal
	arma::mat A = d * arma::mat(n, n, arma::fill::eye);

	// Fills sub- and superdiagonal
	for (int i = 0; i < n-1; i++)
	{
		A(i, i+1) = u;
		A(i+1, i) = l;
	}

	return A;
}


// Initialises a tridiagonal matrix with 'l' on the subdiagonal, 'd' on the diagonal and 'u' on the superdiagonal
arma::mat triMat(const arma::vec* l, const arma::vec* d, const arma::vec* u)
{
	// Initialises A with zeros.
	int n = arma::size(*d)[0];
	arma::mat A = arma::mat(n,n, arma::fill::zeros);

	// Fills the three diagonals of A
	for (int i = 0; i < n-1; i++)
	{
		A(i, i) = (*d)(i);
		A(i, i+1) = (*u)(i);
		A(i+1, i) = (*l)(i);
	}
	A(n-1, n-1) = (*d)(n-1);

	return A;
}


// Writes the three first eigenvectors to file
int printVec(std::string filename, arma::mat& eigvec, double h)
{
	// Sets format of output
	int prec = 6;
	int width = prec + 10;

	int N = eigvec.n_cols;

	// Opens file
	std::ofstream ofile;
	ofile.open("data/" + filename);

	// Writes first endpoint
	for (int i = 0; i < 4; i++)
	{
		ofile << std::setw(width) << std::setprecision(prec) << std::scientific << 0.;
	}
	ofile << std::endl;

	// Writes eigenvectors to file
	for (int i = 0; i < N; i++)
	{
		ofile << std::setw(width) << std::setprecision(prec) << std::scientific << h*(i+1);
		for (int j = 0; j < 3; j++)
		{
			ofile << std::setw(width) << std::setprecision(prec) << std::scientific << eigvec(i,j);
		}
		ofile << std::endl;
	}

	// Writes second endpoint
	ofile << std::setw(width) << std::setprecision(prec) << std::scientific << 1.;
	for (int i = 0; i < 3; i++)
	{
		ofile << std::setw(width) << std::setprecision(prec) << std::scientific << 0.;
	}
	ofile << std::endl;

	// Closes file
	ofile.close();

	return 0;
}