#include "../include/method.h"


using namespace std;

//constants
 vector<vector<double>> A = { {4,1,1},{1,2 * (3 + 0.1 * 2),-1},{1,-1,2 * (4 + 0.1 * 2)} };
 vector<double> B = { 1,-2, 3 };
 vector<double> Xabs = { -0.2541186161449, 0.30683690280065,-0.28912685337726 };

//dot product
double operator*(vector<double>& a, vector<double>& b) {
	double f = 0;
	for (int i = 0;i < b.size();i++) {
		f += a[i] * b[i];
	}
	return f;
}


double Function(vector<double>& x) {
    vector<double> div = (x * A);
    double z = div * x;
    double f = z / 2 + B * x + 2;
    return f;
}

double koef1(vector<double>& q) {
    vector<double> div = (q * A);
    double z = div * q;
    return norm(q) * norm(q) / z;
}

double koef2(vector<double>& q, vector<double>& x) {
    vector<double> div = (q * A);
    double z = div * q;
    vector<double> mno = A * x + B;
    double p = q * mno;
    return p / z;
}

void MNGS4(vector<double>& x) {
    cout << left;
    cout << "Gradient_po_norme_dx" << endl;
    cout << setw(3) << "N" << "|  " << setw(34) << "Xn" << "|  " << setw(20) << "F(Xn)" << "|  " << endl;
    cout << setw(3) << 0 << "|  " << x << "|  " << setw(20) << Function(x) << "|  " << endl;
    int counter = 0;
    bool t = true;
    while (t == true) {
        vector<double> q = A * x + B;
        vector<double> xold = x;
        x = x - q * koef1(q);
        counter += 1;
        cout << setw(3) << counter << "|  " << x << "|  " << setw(15) << Function(x) << "|  " << endl;
        if (norm(x - xold) < 10e-6) {
            t = false;
            cout << "pogreshnost = " << abs(Function(Xabs) - Function(x)) << endl;
        }
    }
}

void MNGS3(vector<double>& x) {
    cout << left;
    cout << "Gradient_po_grad" << endl;
    cout << setw(3) << "N" << "|  " << setw(34) << "Xn" << "|  " << setw(20) << "F(Xn)" << "|  " << endl;
    cout << setw(3) << 0 << "|  " << x << "|  " << setw(20) << Function(x) << "|  " << endl;
    int counter = 0;
    bool t = true;
    while (t == true) {
        vector<double> q = A * x + B;
        vector<double> xold = x;
        x = x - q * koef1(q);
        counter += 1;
        cout << setw(3) << counter << "|  " << x << "|  " << setw(15) << Function(x) << "|  " << endl;
        if (norm(A * x + B) < 10e-6) {
            t = false;
            cout << "pogreshnost = " << abs(Function(Xabs) - Function(x)) << endl;
        }
    }
}
//diagonally dominant
double diagp() {
    int c = 0;
    double preobl = -1;
    for (auto s : A) {
        double d = 0;
        for (int i = 0; i < s.size(); i++) {
            if (i == c) {
                d += abs(s[i]);
            }
            else d -= abs(s[i]);
        }
        c++;
        if (d >= preobl) preobl = d;
    }
    return preobl;
}

void MNGS2(vector<double>& x) {
    double Diag = diagp();
    cout << left;
    cout << "Gradient_s_preobl" << endl;
    cout << setw(3) << "N" << "|  " << setw(34) << "Xn" << "|  " << setw(20) << "F(Xn)" << "|  " << endl;
    cout << setw(3) << 0 << "|  " << x << "|  " << setw(20) << Function(x) << "|  " << endl;
    int counter = 0;
    bool t = true;
    while (t == true) {
        vector<double> q = A * x + B;
        vector<double> xold = x;
        x = x - q * koef1(q);
        counter += 1;
        cout << setw(3) << counter << "|  " << x << "|  " << setw(15) << Function(x) << "|  " << endl;
        if (norm(A * x + B) / Diag < 10e-6) {
            t = false;
            cout << "pogreshnost = " << abs(Function(Xabs) - Function(x)) << endl;
        }
    }
}

int maxdx(vector<double>& x) {
    int max_i = 0;
    for (int i = 1; i < x.size(); i++) {
        if (x[i] > x[max_i]) {
            max_i = i;
        }
    }
    return max_i;
}

void MNPSvibor(vector<double>& x) {
    cout << left;
    cout << "pokoordinat_vibor" << endl;
    cout << setw(3) << "N" << "|  " << setw(34) << "Xn" << "|  " << setw(20) << "F(Xn)" << "|  " << endl;
    cout << setw(3) << 0 << "|  " << x << "|  " << setw(20) << Function(x) << "|  " << endl;
    int counter = 0;
    bool t = true;
    vector<double> dx = { 0,0,0 };
    vector<double> q = { 0,0,0 };
    for (int i = 0; i < 3; i++) {
        q[i] = 1;
        vector<double> xold = x;
        x = x - q * koef2(q, x);
        counter += 1;
        cout << setw(3) << counter << "|  " << x << "|  " << setw(15) << Function(x) << "|  " << endl;
        dx[i] = abs(xold[i] - x[i]);
        q[i] = 0;
    }
    while (t == true) {
        int i = maxdx(dx);
        vector<double> xold = x;
        q[i] = 1;
        x = x - q * koef2(q, x);
        dx[i] = abs(xold[i] - x[i]);
        counter += 1;
        q[i] = 0;
        cout << setw(3) << counter << "|  " << x << "|  " << setw(15) << Function(x) << "|  " << endl;
        if (norm(x - xold) < 10e-6) {
            t = false;
            cout << "pogreshnost = " << abs(Function(Xabs) - Function(x)) << endl;
        }
    }
}

void MNPSperebor(vector<double>& x) {
    cout << left;
    cout << "pokoordinat_perebor" << endl;
    cout << setw(3) << "N" << "|  " << setw(34) << "Xn" << "|  " << setw(20) << "F(Xn)" << "|  " << endl;
    cout << setw(3) << 0 << "|  " << x << "|  " << setw(20) << Function(x) << "|  " << endl;
    int counter = 0;
    bool t = true;
    vector<double> dx = { 0,0,0 };
    vector<double> q = { 0,0,0 };
    for (int i = 0; i < 3; i++) {
        q[i] = 1;
        vector<double> xold = x;
        x = x - q * koef2(q, x);
        counter += 1;
        cout << setw(3) << counter << "|  " << x << "|  " << setw(15) << Function(x) << "|  " << endl;
        dx[i] = abs(xold[i] - x[i]);
        q[i] = 0;
    }

    while (t == true) {
        vector<double> xold = x;
        for (int i = 0; i < 3; i++) {
            q[i] = 1;
            x = x - q * koef2(q, x);
            dx[i] = abs(xold[i] - x[i]);
            counter += 1;
            q[i] = 0;
            cout << setw(3) << counter << "|  " << x << "|  " << setw(15) << Function(x) << "|  " << endl;
        }
        if (norm(x - xold) < 10e-6) {
            t = false;
            cout << "pogreshnost = " << abs(Function(Xabs) - Function(x)) << endl;
        }
    }
}