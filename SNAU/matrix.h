#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include <math.h>
#include <vector>

#define _MATRIX_
#ifdef _MATRIX_

double EPS = 10e-9;

class Matrix{
public:
	int _row, _column;
	double _det;
	double **M = new double*[_row];
	char _name;
	bool _issingular;
	int _kolvop;											//количество перестановок
	Matrix(const int n, const int m, char x)
	{
		_row = n;											// !!
		_column = m;
		_name = x;
		_det = 1.0;
		_issingular = false;
		_kolvop = 0;
		for (int i = 0; i < _row; i++){
			M[i] = new double[_column];
		}

		for (int i = 0; i < _row; i++)
		{
			for (int j = 0; j < _column; j++)
			{
				M[i][j] = 0.0;
			}
		}
	}
	Matrix(const int n, const int m, double* A, char x)
	{
		_name = x;
		_row = n; 
		_column = m;
		_kolvop = 0;
		for (int i = 0; i < _row; i++){
			M[i] = new double[_column];
		}

		for (int i = 0; i < n; i++){
			for (int j = 0; j < m; j++){
				M[i][j] = *(A + i*m + j);
			}
		}
	}
	Matrix ( const Matrix& X)									
	{
		if (_row == X._row && _column == X._column)
		{
			for (int i = 0; i < _row; i++)
			{
				for (int j = 0; j < _column; j++)
				{
					M[i][j] = X.M[i][j];
				}
			}
			_kolvop = X._kolvop;
		}
	}

	~Matrix(){
		//this->print();
		for (int i = 0; i < _row; i++)
			delete[] M[i];
		_kolvop = 0;
		//std::cout << "You killed me! My name is " << _name << std::endl;
	}
	//=============================================================================== перегрузка операторов

	Matrix& operator = (const Matrix& X) {
		//проверка на самоприсваивание
		if (this == &X) {
			return *this;
		}
		_row = X._row;
		_column = X._column;
		for (int count = 0; count < _row; count++)
			delete[] M[count];
		
		for (int i = 0; i < _row; i++){
			M[i] = new double[_column];
		}
		for (int i = 0; i < _row; i++){
			for (int j = 0; j < _column; j++){
				M[i][j] = X.M[i][j];
			}
		}
		return *this;
	}

	//============================================================
	void print(){
		for (int i = 0; i < _row; i++){
			for (int j = 0; j < _column; j++){
				std::cout<<M[i][j]<<" ";
			}
			std::cout << std::endl;
		}
		std::cout<<std::endl;
	}

	void full(const double k)
	{
		for (int i = 0; i < _row; i++)
		{
			for (int j = 0; j < _column; j++)
			{
				M[i][j] = k;
			}
		}
	}

	void IdentityMatrix()
	{
		if (_row != _column) return;
		this->full(0.0);
		for (int i = 0; i < _row; i++)
		{
			M[i][i] = 1.0;
		}
		return;
	}

	friend Matrix& operator + (Matrix&, const Matrix&);
	friend Matrix& operator - (Matrix&, const Matrix&);
	friend Matrix& operator * (const double, Matrix&);
	friend Matrix& operator * (Matrix&, const Matrix&);
	
	void det(){
		if (_row != _column) { 
			std::cout << "Matrix isn't square!"; 
			return; 
		}
		if (_row == 1) {
			_det = M[0][0];
			return;
		}
		else {
			_det = 1.0;
			for (int i = 0; i < _row; i++){
				_det *= M[i][i];
			}
			if(_kolvop%2 != 0){
				_det = -_det;
			}
			return;
		}
	}

	int rank()
	{
		int r = 0;
		
		for (int i = 0; i < _row; i++)
		{
			for (int j = 0; j < _column; j++)
			{
				if (M[i][j]) { 
					++r; 
					break; 
				}
			}
		}
		return r;
	}
	
	friend void PLU(Matrix&, Matrix&, Matrix&, Matrix&, Matrix&);

	void swapRows(int x, int y){
		if (x != y){
			for (int i = 0; i < _column; i++){
				double tmp = M[x][i];
				M[x][i] = M[y][i];
				M[y][i] = tmp;
				//this->print();
			}
		}
		++_kolvop;
		return;
	}
	void swapColumn(int x, int y)
	{
		double tmp;
		for (int i = 0; i < _row; i++){
			tmp = M[i][x];
			M[i][x] = M[i][y];
			M[i][y] = tmp;
			//this->print();
		}
		++_kolvop;
		return;
	}
	double getConditionNumber_1(Matrix& Inv){							//число обусловленности матрицы
		double mtr_n = 0.0;				// норма исходной матрицы
		double inv_n = 0.0;				// норма обратной матрицы
		
		for (int j = 0; j < _column; j++){
			for (int i = 0; i < _row; i++){
				mtr_n += M[i][j]*M[i][j];
			}
		}

		mtr_n = sqrt(mtr_n);

		for (int j = 0; j < _column; j++){
			for (int i = 0; i < _row; i++){
				inv_n += Inv.M[i][j]*Inv.M[i][j];
			}
		}
		inv_n = sqrt(inv_n);

		std::cout << std::endl<< "e) Condition number =" << mtr_n*inv_n<<std::endl;
		return (mtr_n*inv_n);
	}

	friend void solveLinerSystem(const Matrix& , const Matrix& , const Matrix& , Matrix&);
};

Matrix& operator + (Matrix& M1, const Matrix& M2){						
	for (int i = 0; i < M1._row; i++)
	{
		for (int j = 0; j < M1._column; j++)
		{
				M1.M[i][j] += M2.M[i][j];
		}
	}
	return M1;
}
Matrix& operator - (Matrix& M1, const Matrix& M2){						// make all const
	for (int i = 0; i < M1._row; i++)
	{
		for (int j = 0; j < M1._column; j++)
		{
			M1.M[i][j] -= M2.M[i][j];
		}
	}
	return M1;
}

Matrix& operator * (const double left, Matrix& right){			
	for (int i = 0; i < right._row; i++)
	{
		for (int j = 0; j < right._column; j++)
		{
			right.M[i][j] *= left;
		}
	}
	return right;
}

Matrix& operator * (Matrix& left, const Matrix& right)
{
	if (left._column == right._row){
		Matrix T (left._row, right._column, 'T');
		for (int i = 0; i < left._row; i++)
		{
			for (int j = 0; j < right._column; j++)
			{
				for (int r = 0; r < left._column; r++)
				{
					T.M[i][j] += left.M[i][r] * right.M[r][j];
				}
				if (fabs(T.M[i][j]) < EPS){ T.M[i][j] = 0.0; }
			}
		}
		left = T;
	}
	return left;
}

void PLU ( Matrix& A, Matrix&B, Matrix& P, Matrix& L, Matrix& U)			//работает правильно вроде
{
	P.IdentityMatrix();
	B = A;
	for (int i = 0; i < A._row; i++){
		double tmp = 0.0;
		int num = -1;

		for (int row = i; row < A._row; row++) {
			if (fabs(B.M[row][i]) > tmp) {
				tmp = fabs(A.M[row][i]);
				num = row;
			}
		}
		if (num != -1){
			B.swapRows(i, num);
			P.swapRows(i, num);
		}

		if (tmp == 0){
			std::cout << "Matrix is singular"<<std::endl;
			A._issingular = true;
			continue;
		}
		for (int j = i + 1; j < A._row; j++) {
			B.M[j][i] /= B.M[i][i];
			for (int k = i + 1; k < A._row; k++)
				B.M[j][k] -= B.M[j][i] * B.M[i][k];
		}
	}
	// B = L + U - E;
	// L - нижнетреугольная матрица с единицами на диагонали
	// U - верхнетреугольная матрица
	Matrix E(A._row, A._column, 'E');
	E.IdentityMatrix();

	B = B + E;
	L.full(0.0); U.full(0.0);

	for (int i = 0; i < A._row; i++){ L.M[i][i] = 1; }
	for (int i = 0; i < A._row; i++){
		for (int j = 0; j < i; j++){
			L.M[i][j] = B.M[i][j];
		}
	}

	U = B - L;
	std::cout << "PLU decomposition:" << std::endl;
	std::cout << "Matrix P:" << std::endl;
	P.print();
	std::cout << "Matrix L:" << std::endl;
	L.print();
	std::cout << "Matrix U:" << std::endl;
	U.print();
	return;
}

int GaussJordanElimination(Matrix& A, Matrix& b){
	Matrix Tmp_A(A._row, A._column + 1, 'T');

	for (int i = 0; i < A._row; i++){
		for (int j = 0; j < A._column; j++){
			Tmp_A.M[i][j] = A.M[i][j];
		}
	}

	for (int i = 0; i < A._row; i++){
		Tmp_A.M[i][A._column] = b.M[i][0];
	}

	std::vector<int> where(A._column, -1);
	
	const int n = A._row;
	const int m = A._column;

	for (int col = 0, row = 0; col < m && row < n; ++col){
		int sel = row;
		for (int i = row; i < n; i++){
			if (abs(Tmp_A.M[i][col]) > abs(Tmp_A.M[sel][col])){
				sel = i;
			}
		}
		if (abs(Tmp_A.M[sel][col]) < EPS)
			continue;
		for (int i = col; i <= m; i++){
			std::swap(Tmp_A.M[sel][i], Tmp_A.M[row][i]);
		}

		where[col] = row;

		for (int i = 0; i < n; i++){
			if (i != row){
				double tmp = Tmp_A.M[i][col] / Tmp_A.M[row][col];
				for (int j = col; j <= m; ++j){
					Tmp_A.M[i][j] -= Tmp_A.M[row][j] * tmp;
				}
			}
		}
		++row;
	}

	b.full(0.0);

	for (int i = 0; i < m; i++){
		if (where[i] != -1){
			b.M[i][0] = Tmp_A.M[where[i]][m] / Tmp_A.M[where[i]][i];
		}
	}

	for (int i = 0; i < A._row; i++){
		for (int j = 0; j < A._column; j++){
			A.M[i][j] = Tmp_A.M[i][j];
		}
	}
	//Tmp_A.print();

	for (int i = 0; i<n; ++i) {
		double sum = 0.0;
		for (int j = 0; j < m; ++j)
			sum += b.M[j][0] * Tmp_A.M[i][j];
		if (abs(sum - Tmp_A.M[i][m]) > EPS){
			return 0;
		}
	}

	for (int i = 0; i < m; i++){
		if (where[i] == -1){
			return 2;
		}
	}
	return 1;
}

void solveLinerSystem(const Matrix& P, const Matrix& L, const Matrix& U, Matrix& b, Matrix& res){
	const int N = L._row;
	Matrix x(L._row, 1, 'x');
	Matrix y(L._row, 1, 'y');

	double sum;
	double det = 1.0;

	for (int i = 0; i < U._row; i++){
		det *= U.M[i][i];
	}

	if (det){
		Matrix p(N, 1, 'p');

		for (int i = 0; i < N; ++i){
			for (int j = 0; j < N; ++j){
				if (P.M[i][j])
					p.M[i][0] = j;
				else
					p.M[i][0] = 0;
			}
		}

		for (int i = 0; i < N; ++i){
			sum = 0.0;
			for (int j = 0; j < i; ++j){
				sum += L.M[i][j] * y.M[j][0];
			}
			sum = b.M[(int)p.M[i][0]][0] - sum;
			y.M[i][0] = sum;
		}
		for (int i = N - 1; i >= 0; --i){
			sum = 0.0;
			for (int j = i + 1; j < N; ++j){
				sum += U.M[i][j] * x.M[j][0];
			}
			sum = (y.M[i][0] - sum) / U.M[i][i];
			x.M[i][0] = sum;
		}
		res = x;
	}

	return;
}

void getInverseMatrix( const Matrix& P, const Matrix& L, const Matrix& U, Matrix& Inv){
	double det = 1.0;
	for (int i = 0; i < U._row; i++){
		det *= U.M[i][i];
	}

	if (det == 0.0){
		std::cout << "Matrix is singular! Inverse matrix does not exist!" << std::endl;
		return;
	}

	Matrix b(L._row, 1, 'L');
	b.M[0][0] = 1;

	Matrix sol_of_lin_system(L._row, 1, 'S');

	for (int i = 0; i < L._row; i++){
		solveLinerSystem(P, L, U, b, sol_of_lin_system);
		for (int j = 0; j < L._row; j++){
			Inv.M[j][i] = sol_of_lin_system.M[j][0];
		}

		if (i + 1 < L._row){
			b.swapRows(i, i + 1);
		}
	}
	return;
}

#endif