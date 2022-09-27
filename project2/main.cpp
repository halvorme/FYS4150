#include "utils.hpp"
#include "tasks.hpp"

int main()
{
	// Initialise 6x6 matrix A
	int N = 6;
	double h = 1./(N+1.);

	double a = -1./(h*h);
	double d = 2./(h*h);

	arma::vec eigval_anal;
	arma::mat eigvec_anal;

	arma::mat A = triMat(N, a, d, a);


	problem2(N, A, d, a, eigval_anal, eigvec_anal);


	problem3();


	double eps = 1e-8;
	int maxiter = pow(10,6);

	arma::vec eigval_jac(N);
	arma::mat eigvec_jac(N, N);

	problem4(A, eps, maxiter, eigval_jac, eigvec_jac, eigval_anal, eigvec_anal);

	// This number is multiplied by 10
	int N_max = 21; 

	problem5(eps, maxiter, N_max);

	problem6(eps, maxiter);

	return 0;
}
