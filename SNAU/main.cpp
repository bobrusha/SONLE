#include "matrix.h"
using namespace std;

const double ACCURACY = 0.0001f;

double MethodNuton(double, double);
double func(double);
double deriv(double);

double MethodNuton(double , double , double , double , double& , double& );
double func2(double, double);
double func3(double, double);

int main(){
	double solution = MethodNuton(0.0, 1.0f);
	cout << "1) Solution of nonlinear equation:" << endl;
	cout << solution << endl<<endl;
	cout << "Check:" << endl;
	cout << "f("<<solution<<") = "<<func(solution) << endl<<endl;
	
	cout << "Solution of nonlinear system:" << endl;

	double xn = 0, yn = 0 ;

	MethodNuton(3.0, 4.0, 1.0, 2.0, xn, yn);

	std::cout << xn << " " << yn << std::endl;

	std::cout << "Check" << std::endl;
	std::cout << func2(xn, yn) << " " << func3(xn, yn) << std::endl;

	cout << "============" << endl << endl;

	system("pause");
	return 0;
}

double func(double x){
	// f(x) = sqrt(x) - cos(x)
	if ((sqrt(x) - cos(x)) < EPS) return 0.0;
	return (sqrt(x) - cos(x));
}

double deriv(double x){
	// f'(x) = - 1/(2*sqrt(x)) + sin(x)
	if ((1.0 / (2.0*sqrt(x)) + sin(x)) < EPS) return 0.0;
	return (1.0/(2.0*sqrt(x)) + sin(x));
}

double MethodNuton(const double a, const double b){
	double x = b - func(b) / deriv(b);
	if (fabs(x - b) > ACCURACY){
		if (fabs(x - a) > EPS) x = MethodNuton(a, x);
	}
	if (x < EPS){
		x = 0.0;
	}
	return x;
}

double func2(double x, double y){
	return (fabs(cos(x - 1) + y - 0.5f)<EPS ? 0.0 : (cos(x - 1) + y - 0.5f));
}

double func3(double x, double y){
	if (fabs(x - cos(y) - 3.0) < EPS)
		return 0.0;
	else
		return ( x - cos(y) - 3.0);
}

double MethodNuton(double a0, double a1, double c0, double c1, double& xn, double& yn){
	Matrix F(2, 1, 'F');
	Matrix J(2, 2, 'J');
	Matrix Tmp(2, 2, 'T');
	Matrix P(2, 2, 'T');
	Matrix L(2, 2, 'T');
	Matrix U(2, 2, 'T');

	F.M[0][0] = cos(a1 - 1) + c1 - 0.5f;
	F.M[1][0] = a1 - cos(c1) - 3.0;

	F.print();

	J.M[0][0] = -sin(a1 - 1);
	J.M[0][1] = 1;
	J.M[1][0] = 1;
	J.M[1][1] = sin(c1);

	std::cout << "Matrix Jacobi:" << std::endl;
	J.print();
	std::cout << std::endl;

	PLU(J, Tmp, P, L, U);

	Matrix InvJ(2, 2, 'I');

	getInverseMatrix(P, L, U, InvJ);
	cout << "Inverse Matrix:" << endl;
	InvJ.print();

	InvJ*F;
	xn = a1 - InvJ.M[0][0];
	yn = c1 - InvJ.M[1][0];

	if (fabs(xn) < EPS) xn = 0.0;
	if (fabs(yn) < EPS) yn = 0.0;

	std::cout << "============ " << std::endl;
	std::cout << xn << " " << yn << std::endl;

	if (fabs(a1 - xn) > ACCURACY && fabs(c1 - yn) > ACCURACY){
		MethodNuton(a0, xn, c0, yn, xn, yn);
	}

	return 0.0;
}
void MethodNuton(double& x1, double& x2, double& x3, double& x4, double& x5, double& x6, double& x7, double& x8, double& x9, double& x10){
	Matrix F(10, 1, 'F');
	F.M[0][0] = -(cos(x1 * x2) - exp(-3 * x3) + x4 * pow(x5, 2) - x6 - sinh(2 * x8) * x9 + 2 * x10 + 2.0004339741653854440);
	F.M[1][0] = -(sin(x1 * x2) + x3 * x9 * x7 - exp(-x10 + x6) + 3 * pow(x5, 2) - x6 * (x8 + 1) + 10.886272036407019994);
	F.M[2][0] = -(x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8 + x9 - x10 - 3.1361904761904761904);
	F.M[3][0] = -(2 * cos(-x9 + x4) + x5 / (x3 + x1) - sin(pow(x2, 2)) + pow(cos(x7 * x10), 2) - x8 - 0.1707472705022304757);
	F.M[4][0] = -(sin(x5) + 2 * x8 * (x3 + x1) - exp(-x7 * (-x10 + x6)) + 2 * cos(x2) + 1 / (x4 - x9) - 0.3685896273101277862);
	F.M[5][0] = -(exp(x1 - x4 - x9) + pow(x5, 2) / x8 + cos(3 * x10 * x2) / 2 - x6 * x3 + 2.0491086016771875115);
	F.M[6][0] = -(pow(x2, 3) * x7 - sin(x10 / x5 + x8) + (x1 - x6) * cos(x4) + x3 - 0.7380430076202798014);
	F.M[7][0] = -(x5 * pow(x1 - 2 * x6, 2) - 2 * sin(-x9 + x3) + 1.5 * x4 - exp(x2 * x7 + x10) + 3.5668321989693809040);
	F.M[8][0] = -(7 / x6 + exp(x5 + x4) - 2 * x2 * x8 * x10 * x7 + 3 * x9 - 3 * x1 - 8.4394734508383257499);
	F.M[9][0] = -(x10 * x1 + x9 * x2 - x8 * x3 + sin(x4 + x5 + x6) * x7 - 0.78238095238095238096);
	
	Matrix J(7, 10, 'J');

	J.M[0][0] = -sin(x1 * x2) * x2;
	J.M[0][1] = -sin(x1 * x2) * x1;
	J.M[0][2] = 3.0 * exp(-3.0 * x3);
	J.M[0][3] = x5 * x5;
	J.M[0][4] = 2.0 * x4 * x5;
	J.M[0][5] = -1.0;
	J.M[0][6] = 0.0;
	J.M[0][7] = -2.0 * cosh(2.0 * x8) * x9;
	J.M[0][8] = -sinh(2.0 * x8);
	J.M[0][9] = 2.0;
	J.M[1][0] = cos(x1 * x2) * x2;
	J.M[1][1] = cos(x1 * x2) * x1;
	J.M[1][2] = x9 * x7;
	J.M[1][3] = 0.0;
	J.M[1][4] = 6.0 * x5;
	J.M[1][5] = -exp(-x10 + x6) - x8 - 1.0;
	J.M[1][6] = x3 * x9;
	J.M[1][7] = -x6;
	J.M[1][8] = x3 * x7;
	J.M[1][9] = exp(-x10 + x6);
	J.M[2][0] = 1.0;
	J.M[2][1] = -1.0;
	J.M[2][2] = 1.0;
	J.M[2][3] = -1.0;
	J.M[2][4] = 1.0;
	J.M[2][5] = -1.0;
	J.M[2][6] = 1.0;
	J.M[2][7] = -1.0;
	J.M[2][8] = 1.0;
	J.M[2][9] = -1.0;
	J.M[3][0] = -x5 / pow(x3 + x1, 2.0);
	J.M[3][1] = -2.0 * cos(x2 * x2) * x2;
	J.M[3][2] = -x5 / pow(x3 + x1, 2.0);
	J.M[3][3] = -2.0 * sin(-x9 + x4);
	J.M[3][4] = pow(x3 + x1, -1.0);
	J.M[3][5] = 0.0;
	J.M[3][6] = -2.0 * cos(x7 * x10) * sin(x7 * x10) * x10;
	J.M[3][7] = -1.0;
	J.M[3][8] = 2.0 * sin(-x9 + x4);
	J.M[3][9] = -2.0 * cos(x7 * x10) * sin(x7 * x10) * x7;
	J.M[4][0] = 2.0 * x8;
	J.M[4][1] = -2.0 * sin(x2);
	J.M[4][2] = 2.0 * x8;
	J.M[4][3] = pow(-x9 + x4, -2.0);
	J.M[4][4] = cos(x5);
	J.M[4][5] = x7 * exp(-x7 * (-x10 + x6));
	J.M[4][6] = -(x10 - x6) * exp(-x7 * (-x10 + x6));
	J.M[4][7] = 2.0 * x3 + 2.0 * x1;
	J.M[4][8] = -pow(-x9 + x4, -2);
	J.M[4][9] = -x7 * exp(-x7 * (-x10 + x6));
	J.M[5][0] = exp(x1 - x4 - x9);
	J.M[5][1] = -3.0 / 2.0 * sin(3.0 * x10 * x2) * x10;
	J.M[5][2] = -x6;
	J.M[5][3] = -exp(x1 - x4 - x9);
	J.M[5][4] = 2.0 * x5 / x8;
	J.M[5][5] = -x3;
	J.M[5][6] = 0.0;
	J.M[5][7] = -pow(x5, 2.0) / pow(x8, 2.0);
	J.M[5][8] =   -exp(x1 - x4 - x9);
	J.M[5][9] =  -3.0 / 2.0 * sin(3.0 * x10 * x2) * x2;
	J.M[6][0] =  cos(x4);
	J.M[6][1] = 3.0 * x2 * x2 * x7;
	J.M[6][2] = 1.0;
	J.M[6][3] = -(x1 - x6) * sin(x4);
	J.M[6][4] = cos(x10 / x5 + x8) * x10 * pow(x5, -2.0);
	J.M[6][5] = -cos(x4);
	J.M[6][6] = pow(x2, 3.0);
	J.M[6][7] = -cos(x10 / x5 + x8);
	J.M[6][8] = 0.0;
	J.M[6][9] = -cos(x10 / x5 + x8) / x5;
	J.M[7][0] = 2.0 * x5 * (x1 - 2.0 * x6);
	J.M[7][1] = -x7 * exp(x2 * x7 + x10);
	J.M[7][2] = -2.0 * cos(-x9 + x3);
	J.M[7][3] = 1.5;
	J.M[7][4] = pow(x1 - 2.0 * x6, 2.0);
	J.M[7][5] = -4.0 * x5 * (x1 - 2.0 * x6);
	J.M[7][6] = -x2 * exp(x2 * x7 + x10);
	J.M[7][7] = 0.0;
	J.M[7][8] = 2.0 * cos(-x9 + x3);
	J.M[7][9] = -exp(x2 * x7 + x10);
	J.M[8][0] = -3.0;
	J.M[8][1] = -2.0 * x8 * x10 * x7;
	J.M[8][2] = 0.0;
	J.M[8][3] = exp(x5 + x4);
	J.M[8][4] = exp(x5 + x4);
	J.M[8][5] = -7.0 * pow(x6, -2.0);
	J.M[8][6] = -2.0 * x2 * x8 * x10;
	J.M[8][7] = -2.0 * x2 * x10 * x7;
	J.M[8][8] = 3.0;
	J.M[8][9] = -2.0 * x2 * (x8 * x7);
	J.M[9][0] = x10;
	J.M[9][1] = x9;
	J.M[9][2] = -x8;
	J.M[9][3] = cos(x4 + x5 + x6) * x7;
	J.M[9][4] = cos(x4 + x5 + x6) * x7;
	J.M[9][5] = cos(x4 + x5 + x6) * x7;
	J.M[9][6] = sin(x4 + x5 + x6);
	J.M[9][7] = -x3;
	J.M[9][8] = x2;
	J.M[9][9] = x1;

	return;
}