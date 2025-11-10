#ifndef MATRIX_VECTORS_OPERATIONS_H
#define MATRIX_VECTORS_OPERATIONS_H

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

//vector operation
vector<double> operator*(vector<double>& b, double a);

vector<double> operator/(vector<double>b, double a);

vector<double> operator-(vector<double> a, vector<double>b);

vector<double> operator+(vector<double> a, vector<double>b);


//matrix operation
vector<vector<double>> operator*(vector<vector<double>>& b, double a);

vector<vector<double>> operator*(vector<vector<double>>& a, vector<vector<double>>& b);

vector<vector<double>> operator+(vector<vector<double>> a, vector<vector<double>>b);


// vector and matrix operations

vector<double> operator*(vector<vector<double>>& a, vector<double>& b);

vector<double> operator*(vector<double>& b, vector<vector<double>>& a);

//output
ostream& operator << (ostream& stream,vector<double>&a);

ostream& operator <<(ostream& stream, vector<vector<double>>& a);


//transpose
vector<vector<double>> Transp(vector<vector<double>>&A);


//Gauss
vector<double> Gauss(vector<vector<double>>& P, vector<double>& b);


//geron's algorithm for sqrt
float geron(const double& c, double x);


//euclidean norm with geron's sqrt
double En(vector<double> a);


//euclidean norm
double norm(const vector<double>& a);


//swop for Gauss
void swap(int i, int j, vector<vector<double>>& A);
void swapB(int i, int j, vector<double>& A);

#endif