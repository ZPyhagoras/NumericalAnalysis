#include <iostream>
#include <algorithm>
#include <iomanip>
using namespace std;

#define eps 1e-12

class Matrix {
public:
	Matrix(int _row, int _col) {
		row = _row;
		col = _col;
		init_value(row, col);
	}
	Matrix(int _n) {
		row = _n;
		col = _n;
		init_value(row, col);
	}
	void set(int i, int j, double val) {
		nums[i][j] = val;
	}
	double val(int i, int j) {
		return nums[i][j];
	}
	double t_val(int i, int j) {
		return nums[j][i];
	}
	int col_() {
		return col;
	}
	int row_() {
		return row;
	}
	void show() {
		for (int i = 1; i <= row; i++) {
			for (int j = 1; j <= col; j++) {
				cout << "\t\t" << nums[i][j];
			}
			cout << endl;
		}
	}
	Matrix T() {
		Matrix res(col, row);
		for (int i = 1; i <= col; i++) {
			for (int j = 1; j <= row; j++) {
				res.set(i, j, nums[j][i]);
			}
		}
		return res;
	}
	Matrix division(double _val) {
		Matrix res(row, col);
		for (int i = 1; i <= row; i++) {
			for (int j = 1; j <= col; j++) {
				res.set(i, j, nums[i][j] / _val);
			}
		}
		return res;
	}
private:
	void init_value(int _row, int _col) {
		int __row = _row + 1, __col = _col + 1;
		nums = (double**)malloc(sizeof(double*) * __row);
		for (int i = 0; i < __row; i++) {
			nums[i] = (double*)malloc(sizeof(double) * __col);
			for (int j = 0; j < __col; j++) {
				nums[i][j] = 0;
			}
		}
	}
	double** nums;
	int row, col;
};

Matrix operator+(Matrix& _A, Matrix& _B) {
	int ra = _A.row_(), ca = _A.col_();
	int rb = _B.row_(), cb = _B.col_();
	if (ra != rb || ca != cb) {
		return Matrix(0);
	}
	Matrix res(ra, cb);
	for (int i = 1; i <= ra; i++) {
		for (int j = 1; j <= cb; j++) {
			res.set(i, j, _A.val(i, j) + _B.val(i, j));
		}
	}
	return res;
}

Matrix operator-(Matrix& _A, Matrix& _B) {
	int ra = _A.row_(), ca = _A.col_();
	int rb = _B.row_(), cb = _B.col_();
	if (ra != rb || ca != cb) {
		return Matrix(0);
	}
	Matrix res(ra, cb);
	for (int i = 1; i <= ra; i++) {
		for (int j = 1; j <= cb; j++) {
			res.set(i, j, _A.val(i, j) - _B.val(i, j));
		}
	}
	return res;
}

Matrix operator*(Matrix& _A, Matrix& _B) {
	int ra = _A.row_(), ca = _A.col_();
	int rb = _B.row_(), cb = _B.col_();
	if (ca != rb) {
		cout << ca << " " << rb << endl;
		return Matrix(0);
	}
	Matrix res(ra, cb);
	for (int i = 1; i <= ra; i++) {
		for (int j = 1; j <= cb; j++) {
			for (int k = 1; k <= ca; k++) {
				res.set(i, j, res.val(i, j) + _A.val(i, k) * _B.val(k, j));
			}
		}
	}
	return res;
}

int sgn(double val) {
	return val > 0 ? 1 : -1;
}

Matrix Init(int _n) {
	Matrix res(_n);
	for (int i = 1; i <= _n; i++) {
		for (int j = 1; j <= _n; j++) {
			double new_val = 0;
			if (i == j) {
				new_val = (1.64 - 0.024) * sin(0.2 * i) - 0.64 * exp(0.1 / i);
			}
			else if (abs(i - j) == 1) {
				new_val = 0.16;
			}
			else if (abs(i - j) == 2) {
				new_val = -0.064;
			}
			res.set(i, j, new_val);
		}
	}
	return res;
}

Matrix IdentityMatrix(int _n) {
	Matrix res(_n);
	for (int i = 1; i <= _n; i++) {
		for (int j = 1; j <= _n; j++) {
			if (i == j) {
				res.set(i, j, 1);
			}
		}
	}
	return res;
}

Matrix QR(Matrix A) {
	int n = A.col_();
	Matrix Q = IdentityMatrix(n);
	for (int r = 1; r <= n - 1; r++) {
		// (1)
		double zero = 0;
		for (int i = r + 1; i <= n; i++) {
			zero += fabs(A.val(i, r));
		}
		if (zero <= eps) {
			continue;
		}
		// (2)
		double d = 0;
		for (int i = r; i <= n; i++) {
			d += A.val(i, r) * A.val(i, r);
		}
		d = sqrt(d);
		double c = -sgn(A.val(r, r)) * d;
		double h = c * c - c * A.val(r, r);
		// (3)
		Matrix u(n, 1);
		u.set(r, 1, A.val(r, r) - c);
		for (int i = r + 1; i <= n; i++) {
			u.set(i, 1, A.val(i, r));
		}
		Matrix ut = u.T();
		// (4)
		Matrix w = Q * u;
		Matrix tmp = (w * ut).division(h);
		Q = Q - tmp;
		Matrix At = A.T();
		Matrix p = (At * u).division(h), pt = p.T();
		tmp = u * pt;
		A = A - tmp;
	}
	return Q;
}

int main() {
	Matrix A = Init(9);
	for (int i = 0; i < 10000; i++) {
		cout << i << endl;
		Matrix Q = QR(A), Qt = Q.T();
		Matrix oldA = A;
		Matrix tmp = Qt * A;
		A = tmp * Q;
		double zero = 0;
		for (int i = 1; i <= 501; i++) {
			zero += fabs(oldA.val(i, i) - A.val(i, i));
		}
		if (zero <= eps) {
			break;
		}
	}
	double val = 0;
	cout << "???" << endl;
	for (int i = 1; i <= 501; i++) {
		cout << setiosflags(ios::fixed) << scientific
			<< setprecision(12) << A.val(i, i) << endl;
	}
	return 0;
}