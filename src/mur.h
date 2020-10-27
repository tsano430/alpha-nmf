/******************************************************************************
 * mur.h
 ******************************************************************************/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "utils.h"

double _mur_update_W_ik(Matrix *X, Matrix *W, Matrix *H, double *WHT, int i, int k,double alpha);
double _mur_update_H_jk(Matrix *X, Matrix *W, Matrix *H, double *WHT, int j, int k, double alpha);
void _mur_update_W(Matrix *X, Matrix *W, Matrix *H, Matrix *W_new,double alpha, double eps);
void _mur_update_H(Matrix *X, Matrix *W, Matrix *H, Matrix *H_new, double alpha, double eps);
void _mur_update(Matrix *X, Matrix *W, Matrix *H, Matrix *upd_W, Matrix *upd_H, double alpha,double eps);
int mur_solver(Matrix *X, Matrix *W_ret, Matrix *H_ret, int *flag,double alpha, double eps, double delta1,double delta2,int seed, int init_max);
