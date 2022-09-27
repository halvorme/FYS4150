#include <armadillo>
#include <cmath>

#include "tasks.hpp"

int problem2(const int N, const arma::mat A, const double d, const double a, arma::vec& eigval_anal, arma::mat& eigvec_anal)
{
	// Finds the eigenvalues and eigenvectors of A, using Armadillos in-built function 'eig_sym'
	arma::vec eigval;
	arma::mat eigvec;

	arma::eig_sym(eigval, eigvec, A);

	// Fixes the sign of the eigenvectors
	for (int i = 0; i < N; i++)
	{
		eigvec.col(i) = arma::sign(eigvec(0,i)) * eigvec.col(i);
	}

	// Finds the eigenvalues and -vectors from the analytic expressions
	eigval_anal = analEigval(N, d, a);
	eigvec_anal = analEigvec(N, d, a);

	// Checks if the eigenvalues and -vectors are equal
	bool approxval = arma::approx_equal(eigval_anal, eigval, "absdiff", 1e-9);
	bool approxvec = arma::approx_equal(eigvec_anal, eigvec, "absdiff", 1e-9);



	// The rest of the function outputs the results of Problem 2 to terminal

	std::cout << "Problem 2" << std::endl;

	std::cout << "A:" << std::endl;
	A.print();
	std::cout << std::endl;

	std::cout << "Eigenvalues of A: " << std::endl;
	eigval.print();
	std::cout << std::endl;

	std::cout << "Eigenvectors of A (coloumns):" << std::endl;
	eigvec.print();
	std::cout << std::endl;

	std::cout << "The analytic eigenvalues are equal those of Armadillo: ";
	if (approxval){std::cout << "True";}
	else {std::cout << "False";}
	std::cout << std::endl;

	std::cout << "The analytic eigenvectors are equal those of Armadillo: ";
	if (approxvec){std::cout << "True";}
	else {std::cout << "False";}
	std::cout << std::endl << std::endl << std::endl;

	return 0;
}


// Test of the function 'max_offdiag_symmetric'
int problem3()
{
	int N = 4;

	// Initialise A
	arma::mat A(N, N, arma::fill::eye);
	A(0,3) = .5;
	A(1,2) = -.7;
	A(3,0) = A(0,3);
	A(2,1) = A(2,1);

	int k = 0;
	int l = 0;

	// Find maximal off-diagonal value and set (k,l) to its position
	double maxval = max_offdiag_symmetric(A, k, l);

	// Prints results to terminal

	std::cout << "Problem 3" << std::endl << std::endl;

	std::cout << "Maximal entry: ("  << k << "," << l << ")" << std::endl;
	std::cout << "Maximal absolute value: " << maxval << std::endl;

	std::cout << std::endl << std::endl;

	return 0;
}


int problem4(const arma::mat& A, const double eps, const int maxiter, arma::vec& eigval_jac, 
				arma::mat& eigvec_jac, arma::vec& eigval_anal, arma::mat& eigvec_anal)
{
	bool converged;
	int iter;

	jacobi_eigensolver(A, eps, eigval_jac, eigvec_jac, maxiter, iter, converged);

	if (!converged)
	{
		return 1;
	}
	

	bool approxval = arma::approx_equal(eigval_anal, eigval_jac, "absdiff", 1e-9);
	bool approxvec = arma::approx_equal(eigvec_anal, eigvec_jac, "absdiff", 1e-9);


	std::cout << "Problem 4 - Jacobi's method" << std::endl << std::endl;

	std::cout << "Number of iterations: " << iter << std::endl << std::endl;;
	
	std::cout << "Eigenvalues:" << std::endl;
	eigval_jac.print();
	std::cout << std::endl;

	std::cout << "Eigenvectors:" << std::endl;
	eigvec_jac.print();
	std::cout << std::endl;

	std::cout << "The eigenvalues from Jacobi's method are equal to the analytic ones: ";
	if (approxval){std::cout << "True";}
	else {std::cout << "False";}
	std::cout << std::endl;

	std::cout << "The eigenvectors from Jacobi's method are equal to the analytic ones: ";
	if (approxvec){std::cout << "True";}
	else {std::cout << "False";}
	std::cout << std::endl << std::endl << std::endl;

	return 0;
}