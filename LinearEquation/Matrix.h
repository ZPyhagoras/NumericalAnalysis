#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct Matrix* matrix_p;
typedef struct Matrix matrix;

struct Matrix
{
	int row, col;
	double** nums;
};

void matrix_init(int row, int col, matrix_p mat);

void matrix_show(matrix_p mat);