#include "Calculation.h"

int main() {
	int n = 0;
	FILE* stream;
	freopen_s(&stream, "input2.txt", "r", stdin);
	// freopen_s("output.txt", "w", stdout);
	scanf_s("%d", &n);
	matrix A, b, x;
	matrix_init(n, n, &A);
	matrix_init(n, 1, &b);
	matrix_init(n, 0, &x);
	// gauss(n, &A, &b, &x);
	// max_colmn(n, &A, &b, &x);
	doolittle(n, &A, &b, &x);
	return 0;
}