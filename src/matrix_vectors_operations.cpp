#include "../include/matrix_vectors_operations.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

//vector operation
vector<double> operator*(vector<double>& b, double a) {
	vector<double> c;
	for (int i = 0;i < b.size();i++) {
		c.push_back(a*b[i]);
	}
	return c;
}

vector<double> operator/(vector<double>b, double a) {
	vector<double> c;
	for (int i = 0;i < b.size();i++) {
		c.push_back(b[i] / a);
	}
	return c;
}

vector<double> operator-(vector<double> a, vector<double>b) {
	vector<double> c;
	for (int i = 0;i < b.size();i++) {
		c.push_back(a[i] - b[i]);
	}
	return c;
}

vector<double> operator+(vector<double> a, vector<double>b) {
	vector<double> c;
	for (int i = 0;i < b.size();i++) {
		c.push_back(a[i] + b[i]);
	}
	return c;
}

//matrix operation

vector<vector<double>> operator*(vector<vector<double>>& b, double a) {
	vector<vector<double>> c;
	for (int i = 0;i < b.size();i++) {
		c.push_back(b[i] * a);
	}
	return c;
}

vector<vector<double>> operator*(vector<vector<double>>& a, vector<vector<double>>& b) {
	vector<vector<double>> c;
	for (int i = 0;i < a.size();i++) {
		vector<double> str;
		for (int j = 0; j < b[0].size();j++) {
			double f = 0;
			for (int k = 0; k < b.size();k++) {
				f += a[i][k] * b[k][j];
			}
			str.push_back(f);
		}
		c.push_back(str);
	}
	return c;
}

vector<vector<double>> operator+(vector<vector<double>> a, vector<vector<double>>b) {
	vector<vector<double>> c;
	for (int i = 0;i < b.size();i++) {
		c.push_back(a[i] + b[i]);
	}
	return c;
}

// vector and matrix operations

vector<double> operator*(vector<vector<double>>& a, vector<double>& b) {
	vector<double> c;
	for (int i = 0;i < a.size();i++) {
		double f = 0;
		for (int j = 0; j < b.size();j++) {
			f += a[i][j] * b[j];
		}
		c.push_back(f);
	}
	return c;
}

vector<double> operator*(vector<double>& b, vector<vector<double>>& a) {
	vector<double> c;
	for (int i = 0;i < b.size();i++) {
		double f = 0;
		for (int j = 0; j < b.size();j++) {
			f += a[j][i] * b[j];
		}
		c.push_back(f);
	}
	return c;
}

//output
ostream& operator << (ostream& stream,vector<double>&a) {
	stream << "( ";
	for (int i = 0; i < a.size();i++) {
		stream<<setw(8)<< a[i]<<"  ";
	}
	stream << ") ";
	return stream;
}


ostream& operator <<(ostream& stream, vector<vector<double>>& a) {
	int i = 0;
	while (i < a.size()) {
		stream << a[i];
		i += 1;
	}
	cout << endl;
	return stream;
}

//transpose

vector<vector<double>> Transp(vector<vector<double>>&A) {
	vector<vector<double>> C;
	for (int i = 0; i < A[0].size();i++) {
		vector<double> c;
		for (int j = 0; j < A.size();j++) {
			c.push_back(A[j][i]);
		}
		C.push_back(c);
	}
	return C;
}

vector<double> Gauss(vector<vector<double>>& P, vector<double>& b) {

	vector<vector<double>> A = P;
	vector<double> B = b;

	for (int i = 0; i < A.size(); i++) {
		int m = -1;
		if (A[i][i] == 0) {
			for (int l = i + 1; l < A.size();l++) {
				if (A[l][i] != 0) {
					m = l;
				}
				swap(m, i, A);
				swapB(m, i, B);
			}
		}

		B[i] = B[i] / A[i][i];
		A[i] = A[i] / A[i][i];


		for (int j = i + 1; j < A.size(); j++) {
			B[j] = B[j] - B[i] * A[j][i];
			A[j] = A[j] - A[i] * A[j][i];

		}
	}
	for (int p1 = A.size() - 1;p1 >= 0; p1--) {
		B[p1] = B[p1] / A[p1][p1];

		for (int p2 = p1 - 1;p2 >= 0; p2--) {
			B[p2] -= B[p1] * A[p2][p1];
			A[p2][p1] -= A[p2][p1];
		}
	}
	return B;
}

//geron's algorithm for sqrt

float geron(const double& c, double x) {
	float y = (x + c / x) / 2;
	if (fabs(x - y) < (10e-7 / 2.88)) {
		return y;
	}
	else {
		x = y;
		return geron(c, x);
	}
}

//euclidean norm with geron's sqrt

double En(vector<double> a) {
	double S = 0;
	for (int i = 0; i < a.size();i++) {
		S += a[i] * a[i];
	}
	return geron(S, 100);
}

//euclidean norm

double norm(const vector<double>& a) {
	double S = 0;
	for (int i = 0;i < a.size();i++) {
		S += a[i] * a[i];
	}
	return sqrt(S);
 }

//swop for Gauss

void swap(int i, int j, vector<vector<double>>& A) {
	vector<vector<double>> B = A;
	A[i] = B[j];
	A[j] = B[i];
}

void swapB(int i, int j, vector<double>& A) {
	vector<double> B = A;
	A[i] = B[j];
	A[j] = B[i];
}