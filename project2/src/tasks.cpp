#include <armadillo>
#include <cmath>

#include "tasks.hpp"

int problem2()
{
	std::cout << "Problem 2:" << std::endl;
	int N = 6;

	// Initialise A
	arma::mat A = triMat(N, -1., 2., -1.);

	std::cout << "A:" << std::endl;
	A.print();
	std::cout << std::endl;

	// Finds the eigenvalues and eigenvectors of A
	arma::vec eigval;
	arma::mat eigvec;

	arma::eig_sym(eigval, eigvec, A);

	for (int i = 0; i < N; i++)
	{
		if (eigvec(0,i)<0)
		{
			eigvec.col(i) = -eigvec.col(i);
		}
	}

	std::cout << "Eigenvalues of A: ";
	for (int i = 0; i < N-1; i++)
	{
		std::cout << eigval(i) << ", ";
	}
	std::cout << eigval(N-1) << std::endl;

	std::cout << "Eigenvectors of A (coloumns):" << std::endl;
	eigvec.print();
	std::cout << std::endl;

	// Calculates the analytic eigenvalues and -vectors
	arma::vec B = exactEigval(N, 2., -1.);
	arma::mat C = exactEigvec(N, 2., -1.);

	// Checks if the eigenvalues and -vectors are equal
	bool approxval = arma::approx_equal(B, eigval, "absdiff", 1e-14);
	bool approxvec = arma::approx_equal(C, eigvec, "absdiff", 1e-14);

	std::cout << "Equal eigenvalues: ";
	if (approxval){std::cout << "True";}
	else {std::cout << "False";}
	std::cout << std::endl;

	std::cout << "Equal eigenvectors: ";
	if (approxvec){std::cout << "True";}
	else {std::cout << "False";}
	std::cout << std::endl << std::endl;

	return 0;
}


int problem3()
{
	std::cout << "Problem 3:" << std::endl;
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

	std::cout << "Maximal entry: ("  << k << "," << l << ")" << std::endl;
	std::cout << "Maximal value: " << maxval << std::endl;

	std::cout << std::endl;
	
	return 0;
}