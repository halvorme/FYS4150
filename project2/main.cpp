#include "utils.hpp"
#include "tasks.hpp"

int main()
{
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
	int maxiter = 1000;
	// int iter;

	arma::vec eigval_jac(N);
	arma::mat eigvec_jac(N, N);

	problem4(A, eps, maxiter, eigval_jac, eigvec_jac, eigval_anal, eigvec_anal);


	return 0;
}
