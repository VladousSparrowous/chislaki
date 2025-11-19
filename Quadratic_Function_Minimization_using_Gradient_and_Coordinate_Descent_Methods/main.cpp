#include "method.h" 


/*
 * Quadratic Function Minimization using Gradient and Coordinate Descent Methods
 * 
 * Solves: F(x) = 1/2*xᵀAx + bᵀx + c
 * 
 * Methods implemented:
 * - MNGS4: Gradient descent (x-change stopping)
 * - MNGS3: Gradient descent (gradient norm stopping)  
 * - MNGS2: Gradient descent with preconditioning
 * - MNPSvibor: Coordinate descent (max coordinate selection)
 * - MNPSperebor: Coordinate descent (cyclic coordinates)
 * 
 * Matrix A (3x3), vector b provided
 * Compares with exact Gauss solution
 * 
 * Output:
 *  The program prints tables showing iteration number, 
 *  current vector x, and function value f(x) for each method,
 *  along with the final error relative to the exact solution.
 */

using namespace std;

int main(){
	cout << "N = 2" << endl << endl;
	
	cout << "    " << A[0] << endl;
	cout << "A = " << A[1] << endl;
	cout << "    " << A[2] << endl;
	cout << endl;
	cout << "b = " << B << endl;
	vector<double> tochno = Gauss(A, B);
	vector<double> tochnotochno = tochno * (-1);
	cout <<endl<<"Gauss: "<< tochnotochno  << endl<<endl;
	vector<double> x0 = { 0,0,0};
	MNGS4(x0);
	MNGS3(x0);
	MNGS2(x0);
	//x0 = { 3,3,3};
	MNPSvibor(x0);
	MNPSperebor(x0);
	return 0;
}