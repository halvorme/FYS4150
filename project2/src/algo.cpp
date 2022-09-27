#include <armadillo>
#include <cmath>
#include "algo.hpp"


// Calculates the exact eigenvalues from the analytic expression
arma::vec analEigval(const int N, const double d, const double a)
{
	arma::vec eigval(N);

	for (int i = 0; i < N; i++)
	{
		eigval(i) = d + 2*a * std::cos((M_PI * (i+1))/(N+1));
	}

	return eigval;
}


// Calculates the exact eigenvectors from the analytic expression
arma::mat analEigvec(const int N, const double d, const double a)
{
	arma::mat eigvec(N,N);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			eigvec(j,i) = std::sin((M_PI*(i+1)*(j+1))/(N+1));
		}
	}
	// Normalises vectors before they are returned
	return arma::normalise(eigvec);
}


// Performs one Jacobi rotation, on the symmetric matrix 'A'.
// Takes the matrix 'A', rotation matrix 'R' and a set of indices, 'k', 'l'
// Modifies 'A' and 'R' according to the Jacobi method. 
// The entry A(k,l) is the one set to zero.
int jacobi_rotate(arma::mat& A, arma::mat& R, const int k, const int l)
{
	const int N = A.n_cols;

	double a_kl = A(k,l);
	double a_kk = A(k,k);
	double a_ll = A(l,l);

	double a_ki, a_li, r_ik, r_il;

	// tan, cos and sin
	double t, s, c;

	// Finds the trig functions of the rotation angle
	double tau = (a_ll - a_kk)/(2.*a_kl);
	if (tau >= 0.)
	{
		t = 1./(tau + std::sqrt(1. + tau*tau));
	}
	else
	{
		t = -1./(std::sqrt(1. + tau*tau) - tau);
	}
	c = 1./std::sqrt(1. + t*t);
	s = c*t;

	// Updates A
	A(k,k) = c*c*a_kk + s*s*a_ll - 2.*c*s*a_kl;
	A(l,l) = c*c*a_ll + s*s*a_kk + 2.*c*s*a_kl;
	A(k,l) = 0;
	A(l,k) = 0;

	for (int i = 0; i < N; i++)
	{
		if (i != k && i !=l)
		{
			a_ki = A(k,i);
			a_li = A(l,i);
			A(k,i) = c*a_ki - s*a_li;
			A(i,k) = A(k,i);
			A(l,i) = c*a_li + s*a_ki;
			A(i,l) = A(l,i);
		}
		// Updates R
		r_ik = R(i,k);
		r_il = R(i,l);
		R(i,k) = c*r_ik - s*r_il;
		R(i,l) = c*r_il + s*r_ik;
	}

	return 0;
}


// Implements Jacobi's method to find eigenvalues and -vectors of the matrix 'A'.
// Calls on 'jacobi_rotate()' until all off-diagonal elements are smaller than 'eps'.
// If the algorithm coverges in less than 'maxiter' iterations it writes the result 
// to 'eigenvalues' and 'eigenvectors'. 
// Number of iterations are written to the integer "iterations"
int jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, 
						arma::mat& eigenvectors, const int maxiter, int& iter, 
						bool& converged)
{
	int N = A.n_cols;

	// Initialise eigenvector matrix with correct size
	eigenvectors = arma::mat(N, N);

	// Initialise copy of A to work with and rotation matrix R
	arma::mat B = A;
	arma::mat R = arma::mat(N, N, arma::fill::eye);
	int k, l;

	double maxval = max_offdiag_symmetric(B, k, l);
	iter = 0;

	// Runs Jacobi's method until B is diagonal or the maximum number of iterations is reached
	while (maxval > eps && iter < maxiter)
	{
		jacobi_rotate(B, R, k, l);

		maxval = max_offdiag_symmetric(B, k, l);
		iter++;
	}

	// Checks if Jacobi's method was completed
	if (maxval > eps){converged = false;}
	else{converged = true;}

	if (converged)
	{
		// Sorts the eigenvalues and vectors
		arma::uvec order = arma::sort_index(B.diag());
		eigenvalues = sort(B.diag());

		for (int i = 0; i < N; i++)
		{
			eigenvectors.col(i) = arma::sign(R(0, order(i))) * R.col(order(i));
		}
	}
	else
	{
		// Returns error if method did not complete
		std::cout << "Error: Jacobi's method did not converge in " << maxiter << " iterations." << std::endl;
	}

	return 0;
}
