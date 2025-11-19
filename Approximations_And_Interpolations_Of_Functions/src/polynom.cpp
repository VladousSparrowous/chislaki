#include "polynom.h"

using namespace std;

//struct of variable
struct variable {
	double C;
	int power;
};

//variable operations
variable operator+(variable& a, variable& b) {
	variable c = { a.C+b.C, a.power};
	return c;
}

variable operator*(variable a,variable b) {
	variable c = {a.C*b.C, a.power+b.power};
	return c;
}

variable operator*(variable& a, double b) {
	variable c = { a.C*b, a.power};
	return c;
}

variable operator*(double b,variable& a) {
	variable c = { a.C*b, a.power };
	return c;
}

//polynom class 
class Polynom {
	
public:
	
	Polynom(int n) {
		for (int i = 0; i <= n; i++) {
			monom.push_back({0,i});
		}
		N = n;
	}
	Polynom(vector<variable>& a) {
		int n_max=-1;
		for (int i=0;i < a.size();i++) {
			if (a[i].power > n_max) {
				n_max = a[i].power;
			}
		}
		for (int i = 0; i <= n_max; i++) {
			monom.push_back({ 0,i });
		}
		N = n_max;
		for (int i=0;i < a.size();i++) {
			monom[a[i].power] = a[i];
		}
	}
	void in(variable& a) {
		if (a.power<= N) {
			monom[a.power] = monom[a.power]+a;
		}
		else {
			while (a.power>N+1) {
				monom.push_back({ 0, N + 1 });
				N += 1;
			}
			monom.push_back(a);
		}
	}
	variable out(int n) {
		if (n <= N) {
			return monom[n];
		}
		else return{ 0,n };
	}

	int power() {
		return N;
	}
	double calc(double x) {
		double S = monom[0].C;
		double e = x;
		for (int i = 1; i <= N;i++) {
			S += monom[i].C * e;
			e *= x;
		}
		
		return S;
	}
private:
	vector<variable> monom;
	int N;
};

//polynom operations
Polynom operator+(Polynom& a, Polynom& b) {
	int power = 0;
	if (a.power()>= b.power()) {
		power = a.power();
		Polynom c = a;
		for (int i = 0;i <= power;i++) {
		variable s = b.out(i);
		c.in(s);
		}
		return c;
	}
	else {
		power = b.power();
		Polynom c = b;
		for (int i = 0;i <= power;i++) {
			variable s = a.out(i);
			c.in(s);
		}
		return c;
	}
	
}

Polynom operator*(Polynom& a, Polynom& b) {
	
	int power = a.power()+b.power();
	
	Polynom answer(power);
	for (int i = 0;i <= a.power();i++) {
		for (int j = 0;j <= b.power();j++) {
			variable m = a.out(i) * b.out(j);
			answer.in(m);
		}
	}
	return answer;
}

Polynom operator*(Polynom& a, double c) {
	Polynom q(a.power());
	for (int i = 0; i < a.power()+1;i++) {
		variable m = a.out(i);
		variable g = c*m;
		q.in(g);
	}
	return q;
}

//output
ostream& operator<<(ostream& stream, variable& A) {
	if (A.C == 0) {
		return stream;
	}
	if (A.power != 0) {
		A.C > 0 ? stream << "+" << A.C : stream << A.C;
		stream << "x^" << A.power;
	}
	else stream << A.C;

	return stream;
}

ostream& operator<<(ostream& stream,Polynom& A) {
	for (int i = 0; i <= A.power();i++) {
		variable s = A.out(i);
		stream << s;
	}
	return stream;
}
