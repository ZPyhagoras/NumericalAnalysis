#include "Matrix.h"

void matrix_init(int row, int col, matrix_p mat)
{
	mat->row = row;
	mat->col = col;

	printf("Please enter the %d*%d matrix\n", mat->row, mat->col);
	mat->nums = (double**)malloc(sizeof(double*) * row);
	for (int i = 0; i < row; i++)
	{
		mat->nums[i] = (double*)malloc(sizeof(double) * col);
		for (int j = 0; j < col; j++)
		{
			scanf_s("%lf", &mat->nums[i][j]);
		}
	}
}

void matrix_show(matrix_p mat)
{
	for (int i = 0; i < mat->row; i++)
	{
		for (int j = 0; j < mat->col; j++)
		{
			printf("\t%f", mat->nums[i][j]);
		}
		printf("\n");
	}
}
