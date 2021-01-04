/******************************************************************************
 * utils.c: Utility functions
 ******************************************************************************/

#include "utils.h"

/******************************************************************************
 * matrix structure
 *
 * typedef struct {
 *	int row;
 *	int col;
 *	double *data;
 * } Matrix;
 ******************************************************************************/

Matrix *create_matrix(int row_s, int col_s){
	Matrix *mt = (Matrix*)malloc(sizeof(Matrix));
	mt -> row = row_s;
	mt -> col = col_s;
	mt -> data = (double*)malloc(row_s * col_s * sizeof(double));
	return mt;
}

void copy_matrix(Matrix *c_mt, Matrix *mt){
	int i, j;
	
	if (c_mt->row != mt->row || c_mt->col != mt->col){
		printf("copy_matrix: matrix size error\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < c_mt->row; ++i)
		for (j = 0; j < c_mt->col; ++j)
			c_mt->data[i*c_mt->col + j] = mt->data[i*c_mt->col + j];
}

void print_matrix(Matrix *mt){
	int i, j;

	for (i = 0; i < mt->row; ++i){
		for (j = 0; j < mt->col; ++j){
			printf("%.4f ", mt->data[i*mt->col + j]);
		}
		printf("\n");
	}
}

void mul_matrix(Matrix *a, Matrix *b, Matrix *ret){
	int M = a->row, N = b->row, K = a->col;
	double c1 = 1.0, c2 = 0.0;

	// CBLAS: mc = c1 * ma * mb^T + c2 * mc
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 
			c1, a->data, K, b->data, K, c2, ret->data, N);
}

void free_matrix(Matrix *mt){
	free(mt->data);
	free(mt);
}

Matrix *load_matrix(char *file_name, double alpha, double eps){
	char line[MAX_LINE_SIZE+1];
	FILE *fp;
	int row_s, col_s;
	char *row_p, *col_p;
	int i, j;
	Matrix *mt;

	fp = fopen(file_name, "r");

	if (fgets(line,MAX_LINE_SIZE+1,fp) == NULL){
		printf("fgets: error");
		exit(EXIT_FAILURE);
	}
	
	strtok(line, " "); 
	row_p = strtok(NULL, " "); col_p = strtok(NULL, " ");
	row_s = atoi(row_p); col_s = atoi(col_p);

	mt = create_matrix(row_s, col_s);
	
	// assign data
	for(i = 0; i < row_s; ++i){
		// error handling
		if (fgets(line,MAX_LINE_SIZE+1,fp) == NULL){
			printf("fgets: error");
			exit(EXIT_FAILURE);
		}

		double tmp = atof(strtok(line,","));
		if (alpha < 0.0 && tmp < DBL_EPSILON)
			tmp = eps;
		mt->data[col_s*i + 0] = tmp;
		for (j = 1; j < col_s; ++j){
			double tmp = atof(strtok(NULL,","));
			if (alpha < 0.0 && tmp < DBL_EPSILON)
				tmp = eps;
			mt->data[col_s*i + j] = tmp;
		}
	}

	// preprocessing for WDBC dataset
	if (strcmp(file_name, "wdbc.csv") == 0){
		for (i = 0; i < row_s; ++i){
			double col_max = -1;
			// find max value in each column
			for (j = 0; j < col_s; ++j){
				if (mt->data[col_s*i + j] > col_max)
					col_max = mt->data[col_s*i+j];
			}
			// devide each entry by the max value found above
			for (j = 0; j < col_s; ++j)
				mt->data[col_s*i+j] /= col_max;
		}
	}

	fclose(fp);
	return mt;
}

void save_matrix(char *file_name, Matrix *mt){
	int i, j;
	FILE *fp = fopen(file_name, "w");

	for (i = 0; i < mt->row; ++i){
		for (j = 0; j < mt->col; ++j){
			if (j == mt->col-1)
				fprintf(fp,"%.6f\n",mt->data[i*mt->col+j]);
			else
				fprintf(fp, "%.6f,", mt->data[i*mt->col+j]);
		}
	}
	fclose(fp);
}

void init(Matrix *W, Matrix *H, double eps, int seed, int init_max){
	int i,j,k;

	srand48(seed);
	
	// value range: [eps, init_max)
	for (i = 0; i < W->row; i++)
		for (k = 0; k < W->col; ++k)
			W->data[i*W->col + k] = (init_max-eps) * drand48() + eps;

	for (j = 0; j < H->row; ++j)
		for (k = 0; k < H->col; ++k)
			H->data[j*H->col + k] = (init_max-eps) * drand48() + eps;
}

double alpha_div(Matrix *X, Matrix *W, Matrix *H, double alpha){
	int i, j, M = X->row, N = X->col;
	double fi_term = 0.0, se_term = 0.0, th_term = 0.0;
	double ret = 0.0;
	Matrix *mul = create_matrix(M,N);

	mul_matrix(W,H,mul);

	for (i = 0; i < M; ++i)
		for (j = 0; j < N; ++j){
			fi_term += X->data[i*N+j];
			se_term += mul->data[i*N+j];
			th_term += pow(X->data[i*N+j],alpha) * 
				pow(mul->data[i*N+j],1-alpha);
		}
		
	ret = fi_term/(1-alpha) + se_term/alpha - th_term/(alpha*(1-alpha));
	
	free_matrix(mul);

	return ret;
}

double alpha_div_diff1_W(Matrix *X,Matrix *W,Matrix *H,int i,int k, double alpha){
	int j, m, M = X->row, N = X->col, K = W->col;
	double fi_term = 0.0, se_term = 0.0, mul;
	double ret = 0.0;

	for (j = 0; j < N; ++j){
		mul = 0.0;
		for (m = 0; m < K; ++m)
			mul += W->data[i*K+m] * H->data[j*K+m];
		fi_term += H->data[j*K+k];
		se_term += pow(X->data[i*N+j],alpha) * 
			H->data[j*K+k] * pow(mul,(-alpha));
	}
	
	ret = (fi_term - se_term) / alpha;

	return ret;

}

double alpha_div_diff1_H(Matrix *X, Matrix *W, Matrix *H, int j, int k, double alpha){
	int i,m, M = X->row, N = X->col, K = W->col;
	double fi_term = 0.0, se_term = 0.0, mul;
	double ret = 0.0;

	for (i = 0; i < M; ++i){
		mul = 0.0;
		for (m = 0; m < K; ++m)
			mul += W->data[i*K+m] * H->data[j*K+m];
		fi_term += W->data[i*K+k];
		se_term += pow(X->data[i*N+j],alpha) * 
			W->data[i*K+k] * pow(mul, -alpha);
	}
	ret = (fi_term - se_term) / alpha; 

	return ret;
}

double alpha_div_diff2_W(Matrix *X, Matrix *W, Matrix *H, int i, int k, double alpha){
    int j, m, M = X->row, N = X->col, K = W->col;
    double ret = 0.0, mul;
        
    for (j = 0; j < N; ++j){
	    mul = 0.0;
	    for (m = 0; m < K; ++m)
		    mul += W->data[i*K+m] * H->data[j*K+m];

	    ret += pow(X->data[i*N+j],alpha) * H->data[j*K+k] * H->data[j*K+k] * pow(mul,-alpha-1);
    }

    return ret;
}

double alpha_div_diff2_H(Matrix *X, Matrix *W, Matrix *H, int j, int k, double alpha){
    int i,m, M = X->row, N = X->col, K = W->col;
    double ret = 0.0, mul;
    
    for (i = 0; i < M; ++i){
	    mul = 0.0;
	    for (m = 0; m < K; ++m)
		    mul += W->data[i*K+m] * H->data[j*K+m];
	    
	    ret += pow(X->data[i*N+j],alpha) * W->data[i*K+k] * W->data[i*K+k] * pow(mul,-alpha-1);
    }
    
    return ret;
}

int stop(Matrix *X, Matrix *W, Matrix *H,double alpha, double eps, double delta1, double delta2){
	int M = X->row, N = X->col, K = W->col, i, j, k;
	double diff1_W, diff1_H;

	for (i = 0; i < M; ++i)
		for (k = 0; k < K; ++k){
			diff1_W = alpha_div_diff1_W(X,W,H,i,k,alpha);
			if (diff1_W > delta1 && W->data[i*K + k] <= eps+delta2)
				;
			else if (fabs(diff1_W) <= delta1)
				;
			else
				return 0; // False
		}

	for (j = 0; j < N; ++j)
		for (k = 0; k < K; ++k){
			diff1_H = alpha_div_diff1_H(X,W,H,j,k,alpha);
			if (diff1_H > delta1 && H->data[j*K + k] <= eps+delta2)
				;
			else if (fabs(diff1_H) <= delta1)
				;
			else
				return 0; // False
		}

	return 1; // True
}

