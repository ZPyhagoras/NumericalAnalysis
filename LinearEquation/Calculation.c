#include "Calculation.h"

void equation_show(int n, matrix_p A, matrix_p b)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("\t%5f", A->nums[i][j]);
		}
		printf("\t%5f\n", b->nums[i][0]);
	}
	printf("\n");
}

void backtrace(int n, matrix_p A, matrix_p b, matrix_p x) 
{
	x->nums[n - 1][0] = b->nums[n - 1][0] / A->nums[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < n; j++)
		{
			sum += A->nums[i][j] * x->nums[j][0];
		}
		x->nums[i][0] = (b->nums[i][0] - sum) / A->nums[i][i];
	}
	for (int i = 0; i < n; i++)
	{
		printf("%f\n", x->nums[i][0]);
	}
}

void gauss(int n, matrix_p A, matrix_p b, matrix_p x)
{
	equation_show(n, A, b);
	double m_ik = 0;
	// 顺序消去
	for (int k = 0; k < n - 1; k++)
	{
		for (int i = k + 1; i < n; i++)
		{
			m_ik = A->nums[i][k] / A->nums[k][k];
			for (int j = 0; j < k + 1; j++) 
			{
				A->nums[i][j] = 0;
			}
			for (int j = k + 1; j < n; j++) 
			{
				A->nums[i][j] -= m_ik * A->nums[k][j];
			}
			b->nums[i][0] -= m_ik * b->nums[k][0];
		}
		equation_show(n, A, b);
	}
	// 回代
	backtrace(n, A, b, x);
}

void max_colmn(int n, matrix_p A, matrix_p b, matrix_p x)
{
	equation_show(n, A, b);
	double m_ik = 0;
	// 列主元素消去
	for (int k = 0; k < n - 1; k++)
	{
		// 寻找主元素
		int max_idx = k;
		for (int i = k + 1; i < n; i++)
		{
			if (fabs(A->nums[i][k]) > fabs(A->nums[max_idx][k]))
			{
				max_idx = i;
			}
		}
		// 交换主元素所在行
		double tmp = 0;
		for (int j = 0; j < n; j++) 
		{
			tmp = A->nums[max_idx][j];
			A->nums[max_idx][j] = A->nums[k][j];
			A->nums[k][j] = tmp;
		}
		tmp = b->nums[max_idx][0];
		b->nums[max_idx][0] = b->nums[k][0];
		b->nums[k][0] = tmp;
		for (int i = k + 1; i < n; i++)
		{
			m_ik = A->nums[i][k] / A->nums[k][k];
			for (int j = 0; j < k + 1; j++)
			{
				A->nums[i][j] = 0;
			}
			for (int j = k + 1; j < n; j++)
			{
				A->nums[i][j] -= m_ik * A->nums[k][j];
			}
			b->nums[i][0] -= m_ik * b->nums[k][0];
		}
		equation_show(n, A, b);
	}
	// 回代
	backtrace(n, A, b, x);
}

void doolittle(int n, matrix_p A, matrix_p b, matrix_p x)
{
	for (int k = 0; k < n; k++)
	{
		for (int j = k; j < n; j++) 
		{
			double sum_lu = 0;
			for (int t = 0; t < k; t++) 
			{
				sum_lu += A->nums[k][t] * A->nums[t][j];
			}
			A->nums[k][j] = A->nums[k][j] - sum_lu;
		}
		for (int i = k + 1; i < n; i++) 
		{
			double sum_lu = 0;
			for (int t = 0; t < k; t++)
			{
				sum_lu += A->nums[i][t] * A->nums[t][k];
			}
			A->nums[i][k] = (A->nums[i][k] - sum_lu) / A->nums[k][k];
		}
	}
	equation_show(n, A, b);

	matrix Y;
	matrix_init(n, 0, &Y);
	matrix_p y = &Y;
	y->nums[0][0] = b->nums[0][0];
	for (int i = 1; i < n; i++) 
	{
		double sum_ly = 0;
		for (int t = 0; t < i; t++)
		{
			sum_ly += A->nums[i][t] * y->nums[t][0];
		}
		y->nums[i][0] = b->nums[i][0] - sum_ly;
	}
	x->nums[n - 1][0] = y->nums[n - 1][0] / A->nums[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) 
	{
		double sum_ux = 0;
		for (int t = i + 1; t < n; t++)
		{
			sum_ux += A->nums[i][t] * x->nums[t][0];
		}
		x->nums[i][0] = (y->nums[i][0] - sum_ux) / A->nums[i][i];
	}
	for (int i = 0; i < n; i++)
	{
		printf("%f\n", y->nums[i][0]);
	}
	for (int i = 0; i < n; i++)
	{
		printf("%f\n", x->nums[i][0]);
	}
}

void crout(int n, matrix_p A, matrix_p b, matrix_p x)
{
}
