/******************************************************************************
 * proposed.h
 ******************************************************************************/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "utils.h"

double _proposed_update_W_ik(Matrix *X,Matrix *W,Matrix *H,int i,int k, double alpha,double eps,double delta1,double delta2);
double _proposed_update_H_jk(Matrix *X,Matrix *W,Matrix *H,int j,int k, double alpha,double eps,double delta1,double delta2);
void _proposed_update_W(Matrix *X, Matrix *W, Matrix *H, Matrix *W_new, int *flag,double alpha,double eps,double delta1,double delta2);
void _proposed_update_H(Matrix *X, Matrix *W, Matrix *H, Matrix *H_new, int *flag,double alpha,double eps,double delta1,double delta2);
void _proposed_update(Matrix *X, Matrix *W, Matrix *H, Matrix *W_new, Matrix *H_new, int *flag,double alpha,double eps,double delta1,double delta2);
int proposed_solver(Matrix *X, Matrix *W_ret, Matrix *H_ret,int *flag,double alpha, double eps,double delta1,double delta2,int seed,int init_max);

