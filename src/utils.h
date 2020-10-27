/******************************************************************************
 * utils.h
 ******************************************************************************/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <cblas.h>
#include <lapacke.h>

#define MAX_LINE_SIZE 50000

/******************************************************************************
 * matrix structure
 ******************************************************************************/

typedef struct {
	int row;
	int col;
	double *data;
} Matrix;

/******************************************************************************
 * functions
 ******************************************************************************/

Matrix *create_matrix(int row_s, int col_s);
void copy_matrix(Matrix *c_mt, Matrix *mt);
void print_matrix(Matrix *mt);
void mul_matrix(Matrix *a, Matrix *b, Matrix *ret);
void inner_product(double *a, double *b, double *ret, int size);
void free_matrix(Matrix *mt);
Matrix *load_matrix(char *file_name, double alpha, double eps);
void save_matrix(char *file_name, Matrix *mt);
void init(Matrix *W, Matrix *H, double eps, int seed, int init_max);
double alpha_div(Matrix *X, Matrix *W, Matrix *H, double alpha);
double alpha_div_diff1_W(Matrix *X,Matrix *W,Matrix *H,int i,int k, double alpha);
double alpha_div_diff1_H(Matrix *X, Matrix *W, Matrix *H, int j, int k, double alpha);
double alpha_div_diff2_W(Matrix *X, Matrix *W, Matrix *H, int i, int k, double alpha);
double alpha_div_diff2_H(Matrix *X, Matrix *W, Matrix *H, int j, int k, double alpha);
int stop(Matrix *X, Matrix *W, Matrix *H,double alpha, double eps, double delta1,double delta2);