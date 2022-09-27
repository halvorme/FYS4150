#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "tasks.hpp"
#include "utils.hpp"

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


// Tests the implementation of Jacobi's method
int problem4(const arma::mat& A, const double eps, const int maxiter, arma::vec& eigval_jac, 
				arma::mat& eigvec_jac, arma::vec& eigval_anal, arma::mat& eigvec_anal)
{
	bool converged;
	int iter;

	// Computes algorithm
	jacobi_eigensolver(A, eps, eigval_jac, eigvec_jac, maxiter, iter, converged);

	// Stop if method did not converge
	if (!converged)
	{
		return 1;
	}

	// Compares numerical result to analytic solution
	bool approxval = arma::approx_equal(eigval_anal, eigval_jac, "absdiff", 1e-9);
	bool approxvec = arma::approx_equal(eigvec_anal, eigvec_jac, "absdiff", 1e-9);

	// Prints results to terminal
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


// Tests the efficiency of the algorithm, by running it for different values of N
int problem5(const double eps, const int maxiter, const int N_max)
{
	int N;
	double h, a, d;

	int iter;
	bool converged;

	arma::vec eigval;
	arma::mat eigvec;

	std::ofstream ofile;
	ofile.open("data/tri_iter.txt");

	int width = 10;

	// Tests algorithm for tridiagonal matrix
	for (int i = 1; i < N_max; i++)
	{
		N = 10*i;
		h = 1./(N+1.);

		a = -1./(h*h);
		d = 2./(h*h);

		arma::mat A = triMat(N, a, d, a);

		jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iter, converged);

		ofile << std::setw(width) << N;
		ofile << std::setw(width) << iter;
		ofile << std::endl;
	}

	ofile.close();

	// Tests algorithm for dense matrix
	ofile.open("data/dense_iter.txt");
	for (int i = 1; i < N_max; i++)
	{
		N = 10*i;
		h = 1./(N+1.);

		a = -1./(h*h);
		d = 2./(h*h);

		// Generate random symmetric matrix
		arma::mat A = arma::mat(N, N).randn();  
		A = arma::symmatu(A);

		jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iter, converged);

		ofile << std::setw(width) << N;
		ofile << std::setw(width) << iter;
		ofile << std::endl;
	}

	ofile.close();

	return 0;
}


int problem6(const double eps, const int maxiter)
{
	int N = 9;
	double h, a, d;

	int iter;
	bool converged;

	arma::vec eigval;
	arma::mat eigvec;

	h = 1./(N+1.);

	a = -1./(h*h);
	d = 2./(h*h);

	arma::mat A = triMat(N, a, d, a);

	// Finds numerical solution
	jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iter, converged);

	// Finds analytic solution
	arma::mat eigvec_anal = analEigvec(N, d, a);
	
	std::string filename;

	// Writes numerical solution to file
	filename = "x_v10.txt";
	printVec(filename, eigvec, h);

	// Writes analytic solution to file
	filename = "x_u10.txt";
	printVec(filename, eigvec_anal, h);

	// Same as above with different N (PS: Time was up...)
	N = 99;

	h = 1./(N+1.);

	a = -1./(h*h);
	d = 2./(h*h);

	A = triMat(N, a, d, a);

	jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iter, converged);

	eigvec_anal = analEigvec(N, d, a);
	
	filename = "x_v100.txt";
	printVec(filename, eigvec, h);

	filename = "x_u100.txt";
	printVec(filename, eigvec_anal, h);

	return 0;
}
