#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

typedef double** Matrix;     // 定义数据类型Matrix, 二维浮点型指针
typedef double* Vector;      // 定义数据类型Vector, 一维浮点型指针
const int max_iter = 10000;  // 定义迭代最大次数
const double eps = 1e-12;    // 定义代码中的浮点数精度

/*
函数名：CreateMatrix
功能：根据输入行列的参数，产生对应大小的矩阵，返回矩阵指针
返回类型：Matrix(double**)
*/
Matrix CreateMatrix(int row, int col) {
	Matrix mat = (double**)malloc(sizeof(double*) * row);  // 根据矩阵行数分配对应数量的空间
	for (int i = 0; i < row; i++) {
		mat[i] = (double*)malloc(sizeof(double) * col);  // 根据矩阵列行数分配对应数量的空间
		for (int j = 0; j < col; j++) {
			mat[i][j] = 0;  // 将所有元素初始化为0
		}
	}
	return mat;
}

/*
函数名：CopyMatrix
功能：复制输入的矩阵，返回新的矩阵指针，新旧矩阵不在同一地址空间
返回类型：Matrix(double**)
*/
Matrix CopyMatrix(Matrix mat, int row, int col) {
	Matrix new_mat = (double**)malloc(sizeof(double*) * row);  // 根据矩阵行数分配对应数量的空间
	for (int i = 0; i < row; i++) {
		new_mat[i] = (double*)malloc(sizeof(double) * col);  // 根据矩阵列行数分配对应数量的空间
		for (int j = 0; j < col; j++) {
			new_mat[i][j] = mat[i][j];  // 新矩阵元素与旧矩阵元素相同
		}
	}
	return new_mat;
}

/*
函数名：Init
功能：根据题目要求对矩阵进行初始化，采用压缩存储的方式，不存储0元素
返回类型：void
*/
void Init(Matrix mat, int n, int r, int s) {
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			if (i == j) {  // 根据要求初始化主对角线元素a
				mat[i - j + s + 1][j] = (1.64 - 0.024 * i) * sin(0.2 * i) - 0.64 * exp(0.1 / i);
			}
			else if (abs(i - j) == 1) {  // 初始化b
				mat[i - j + s + 1][j] = 0.16;
			}
			else if (abs(i - j) == 2) {  // 初始化c
				mat[i - j + s + 1][j] = -0.064;
			}
		}
	}
}

/*
函数名：Panning
功能：对矩阵进行平移，将输入矩阵减去相应倍数的单位矩阵，直接修改原矩阵
返回类型：void
*/
void Panning(Matrix mat, double val, int n, int r, int s) {
	for (int j = 1; j <= n; j++) {
		mat[s + 1][j] -= val;  // 主对角线元素减去相应值
	}
}

/*
函数名：PowerMethod
功能：根据幂法计算输入矩阵模最大的特征值，返回模最大的特征值
返回类型：double
*/
double PowerMethod(Matrix mat, int n, int r, int s) {
	Vector u = (double*)malloc(sizeof(double) * n + 1);  // 分配u和y的地址空间，并初始化
	Vector y = (double*)malloc(sizeof(double) * n + 1);
	for (int i = 0; i <= n; i++) {
		u[i] = y[i] = 1;
	}
	double error = 1, beta = 0;  // 初始化误差为1，特征值为0
	for (int k = 1; k <= max_iter && error > eps; k++) {  // 根据公式进行迭代，直到达到最大迭代次数或误差小于精度要求时停止
		double eta = 0;  // 根据公式计算eta
		for (int i = 1; i <= n; i++) {
			eta += u[i] * u[i];
		}
		eta = sqrt(eta);
		for (int i = 1; i <= n; i++) {  // 迭代更新特征向量
			y[i] = u[i] / eta;
		}
		for (int i = 1; i <= n; i++) {
			double tmp = 0;
			for (int j = 1; j <= n; j++) {
				double val = abs(i - j) > 2 ? 0 : mat[i - j + s + 1][j];
				tmp += val * y[j];
			}
			u[i] = tmp;
		}
		double new_beta = 0;  // 计算新的特征值
		for (int i = 1; i <= n; i++) {
			new_beta += y[i] * u[i];
		}
		error = fabs(new_beta - beta) / fabs(new_beta);  // 计算当前特征值迭代误差
		beta = new_beta;
	}
	return beta;
}

/*
函数名：BandedLU
功能：对带状矩阵进行LU分解，直接修改原矩阵
返回类型：void
*/
void BandedLU(Matrix mat, int n, int r, int s) {
	for (int k = 1; k <= n; k++) {
		int col = k + s < n ? k + s : n;
		for (int j = k; j <= col; j++) {  // 根据公式更新每一行的值
			int t = k - r > j - s ? k - r : j - s;
			double tmp = 0;
			for (t = t > 1 ? t : 1; t <= k - 1; t++) {
				tmp += mat[k - t + s + 1][t] * mat[t - j + s + 1][j];
			}
			mat[k - j + s + 1][j] = mat[k - j + s + 1][j] - tmp;
		}
		int row = k + r < n ? k + r : n;
		for (int i = k + 1; i <= row && k < n; i++) {  // 根据公式更新每一列的值
			int t = i - r > k - s ? i - r : k - s;
			double tmp = 0;
			for (t = t > 1 ? t : 1; t <= k - 1; t++) {
				tmp += mat[i - t + s + 1][t] * mat[t - k + s + 1][k];
			}
			mat[i - k + s + 1][k] = (mat[i - k + s + 1][k] - tmp) / mat[s + 1][k];
		}
	}
}

/*
函数名：InversePowerMethod
功能：根据反幂法计算输入矩阵模最小的特征值，返回模最小的特征值
返回类型：double
*/
double InversePowerMethod(Matrix mat, int n, int r, int s) {
	Matrix lu_mat = CopyMatrix(mat, r + s + 1 + 1, n + 1);  // 复制矩阵并进行LU分解
	BandedLU(lu_mat, n, r, s);
	Vector u = (double*)malloc(sizeof(double) * n + 1);  // 分配u和y的地址空间，并初始化
	Vector y = (double*)malloc(sizeof(double) * n + 1);

	for (int i = 0; i <= n; i++) {
		u[i] = y[i] = 1;
	}
	double error = 1, beta = 0;  // 初始化误差为1，特征值为0
	for (int k = 1; k <= max_iter && error > eps; k++) {// 根据公式进行迭代，直到达到最大迭代次数或误差小于精度要求时停止
		double eta = 0;  // 根据公式计算eta
		for (int i = 1; i <= n; i++) {
			eta += u[i] * u[i];
		}
		eta = sqrt(eta);
		for (int i = 1; i <= n; i++) {  // 迭代更新特征向量
			y[i] = u[i] / eta;
		}
		Vector b = (double*)malloc(sizeof(double) * n + 1);
		for (int i = 0; i <= n; i++) {
			b[i] = y[i];
		}
		for (int i = 2; i <= n; i++) {
			double tmp = 0;
			for (int t = i - r > 1 ? i - r : 1; t <= i - 1; t++) {
				tmp += lu_mat[i - t + s + 1][t] * b[t];
			}
			b[i] = b[i] - tmp;
		}
		u[n] = b[n] / lu_mat[s + 1][n];  // 根据LU分解计算向量u
		for (int i = n - 1; i >= 1; i--) {
			double tmp = 0;
			for (int t = i + 1; t <= (i + s < n ? i + s : n); t++) {
				tmp += lu_mat[i - t + s + 1][t] * u[t];
			}
			u[i] = (b[i] - tmp) / lu_mat[s + 1][i];
		}

		double new_beta = 0;  // 计算新的特征值
		for (int i = 1; i <= n; i++) {
			new_beta += y[i] * u[i];
		}
		error = fabs(new_beta - beta) / fabs(new_beta);  // 计算当前特征值迭代误差
		beta = new_beta;
	}
	return beta;
}

/*
函数名：test
功能：生成测试矩阵，验证计算是否正确
返回类型：void
*/
void test() {
	Matrix test = (double**)malloc(sizeof(double*) * 4);
	for (int i = 0; i < 4; i++) {
		test[i] = (double*)malloc(sizeof(double) * 4);
		for (int j = 0; j < 4; j++) {
			test[i][j] = 0;
		}
	}
	test[1][1] = 5; test[1][2] = 30; test[1][3] = -48;
	test[2][1] = 3; test[2][2] = 14; test[2][3] = -24;
	test[3][1] = 3; test[3][2] = 15; test[3][3] = -25;

	//Matrix test = (double**)malloc(sizeof(double*) * 11);
	//for (int i = 0; i < 11; i++) {
	//	test[i] = (double*)malloc(sizeof(double) * 11);
	//	for (int j = 0; j < 11; j++) {
	//		test[i][j] = 0;
	//	}
	//}
	//for (int i = 1; i <= 10; i++) {
	//	for (int j = 1; j <= 10; j++) {
	//		test[i][j] = i == j ? 1.52 * cos(i + 1.2 * j) : sin(0.5 * i + 0.2 * j);
	//	}
	//}
}

int main() {
	int n = 501, r = 2, s = 2;
	Matrix A = CreateMatrix(r + s + 1 + 1, n + 1);
	Init(A, n, r, s);

	cout << setiosflags(ios::fixed) << scientific << setprecision(12);

	double lambda_1 = PowerMethod(A, n, r, s);  // 使用幂法，计算模最大的特征值，即lambda_1
	cout << "lambda_1\t= " << lambda_1 << endl;

	Panning(A, lambda_1, n, r, s);  // 平移矩阵A，减去模最大的特征值
	double lambda_501 = PowerMethod(A, n, r, s) + lambda_1;  // 计算当前模最大的特征值，即lambda_501
	Panning(A, -lambda_1, n, r, s);  // 加上模最大的特征值，恢复原始矩阵A
	cout << "lambda_501\t= " << lambda_501 << endl;

	double lambda_s = 1.0 / InversePowerMethod(A, n, r, s);  // 使用反幂法，计算模最小的特征值
	cout << "lambda_s\t= " << lambda_s << endl;

	double tmp = (lambda_501 - lambda_1) / 40;
	for (int k = 1; k <= 39; k++) {
		double mu_k = lambda_1 + k * tmp;  // 计算全部mu_k
		Panning(A, mu_k, n, r, s);  // 平移矩阵A，减去mu_k
		double lambda_ik = 1.0 / InversePowerMethod(A, n, r, s) + mu_k;  // 使用反幂法计算平移后矩阵模最小的特征值，即与mu_k最接近的特征值
		cout << "lambda_i_" << k << "\t= " << lambda_ik << endl;
		Panning(A, -mu_k, n, r, s);  // 加上mu_k，恢复原始矩阵A
	}

	double cond_2 = fabs(lambda_1 / lambda_s);  // 计算矩阵A的谱范数
	cout << "cond_2\t= " << cond_2 << endl;

	Matrix LU = CopyMatrix(A, r + s + 1 + 1, n + 1);  // 对矩阵A进行LU分解，此时分解后的矩阵主对角线元素乘积与矩阵A的行列式相等
	BandedLU(LU, n, r, s);
	double detA = 1.0;
	for (int j = 1; j <= n; j++) {  // 主对角线元素相乘得到A的行列式
		detA *= LU[s + 1][j];
	}
	cout << "detA\t= " << detA << endl;

	return 0;
}