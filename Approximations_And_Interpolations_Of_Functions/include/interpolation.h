#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "polynom.h"
#include <vector>
using namespace std;

// constants
extern const double pi;
extern const double e;
extern const double a;
extern const double b;

// Functions
double F(double x);
double F2(double x);
double exp(double x);
double Max(vector<double>& X);

// Interpolation
Polynom Ln(int n, bool Kakoi);
Polynom Nn(int n, bool Kakoi);
Polynom Nn_drugoi(int n, bool Kakoi);
Polynom Pn(int n, bool Kakoi);

// least squares method
Polynom MNKnorm(int n, vector<double>& X, vector<double>& Y);
vector<Polynom> MNKq(int n, vector<double>& X, vector<double>& Y);

#endif