#pragma once

#include "Matrix.h"

void equation_show(int n, matrix_p A, matrix_p b);

void backtrace(int n, matrix_p A, matrix_p b, matrix_p x);

void gauss(int n, matrix_p A, matrix_p b, matrix_p x);

void max_colmn(int n, matrix_p A, matrix_p b, matrix_p x);

void doolittle(int n, matrix_p A, matrix_p b, matrix_p x);

void crout(int n, matrix_p A, matrix_p b, matrix_p x);