#pragma once
#include "matrix_vectors_operations.h" 
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

extern vector<vector<double>> A;
extern vector<double> B;
extern vector<double> Xabs;

double operator*(vector<double>& a, vector<double>& b);

double Function(vector<double>& x);
double koef1(vector<double>& q);
double koef2(vector<double>& q, vector<double>& x);
void MNGS4(vector<double>& x);
void MNGS3(vector<double>& x);
double diagp();
void MNGS2(vector<double>& x);
int maxdx(vector<double>& x);
void MNPSvibor(vector<double>& x);
void MNPSperebor(vector<double>& x);