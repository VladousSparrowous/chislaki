#ifndef POLYNOM_H
#define POLYNOM_H

#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>


using namespace std;

//struct of variable
struct variable {
	double C;
	int power;
};

//variable operations
variable operator+(variable& a, variable& b);

variable operator*(variable a,variable b);

variable operator*(variable& a, double b);

variable operator*(double b,variable& a);

//polynom class 
class Polynom {
	
public:
	Polynom(int n);
	Polynom(vector<variable>& a);
	void in(variable& a);
	variable out(int n);

	int power();
	double calc(double x);
private:
	vector<variable> monom;
	int N;
};

//polynom operations
Polynom operator+(Polynom& a, Polynom& b) ;

Polynom operator*(Polynom& a, Polynom& b);

Polynom operator*(Polynom& a, double c);

//output
ostream& operator<<(ostream& stream, variable& A);

ostream& operator<<(ostream& stream,Polynom& A);

#endif