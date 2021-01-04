/******************************************************************************
 * mur.c: MUR algorithm
 ******************************************************************************/

#include "mur.h"

#define LOG_FLAG 0

double _mur_update_W_ik(Matrix *X, Matrix *W, Matrix *H, double *WHT, int i, int k,double alpha){
	double nume = 0.0, deno = 0.0, Z, W_ik_new;
	int j, M = X->row, N = X->col, K = W->col;

	for (j = 0; j < N; ++j){
		Z = X->data[i*N + j] / WHT[i*N + j];
		nume += pow(Z,alpha) * H->data[j*K + k];
		deno += H->data[j*K + k];
	}
	W_ik_new = W->data[i*K + k] * pow(nume/deno, 1/alpha);
	
	return W_ik_new;
}

double _mur_update_H_jk(Matrix *X, Matrix *W, Matrix *H, double *WHT, int j, int k, double alpha){
	double nume = 0.0, deno = 0.0, Z, H_jk_new;
	int i, M = X->row, N = X->col, K = W->col;

	for (i = 0; i < M; ++i){
		Z = X->data[i*N + j] / WHT[i*N + j];
		nume += pow(Z,alpha) * W->data[i*K + k];
		deno += W->data[i*K + k];
	}
	H_jk_new = H->data[j*K + k] * pow(nume/deno, 1/alpha);

	return H_jk_new;
}

void _mur_update_W(Matrix *X, Matrix *W, Matrix *H, Matrix *W_new,double alpha, double eps){
	int M = X->row, N = X->col, K = W->col, i, k;
	double w_tmp, diff1;
	double *WHT = (double *)malloc(M*N*sizeof(double));

	//mul_matrix(W,H,mul);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, W->data, K, H->data, K, 0.0, WHT, N);

	for (i = 0; i < M; ++i){
		for (k = 0; k < K; ++k){
			w_tmp = _mur_update_W_ik(X,W,H,WHT,i,k,alpha);
			W_new->data[i*K + k] = fmax(eps, w_tmp);
		}
	}

	free(WHT);
}


void _mur_update_H(Matrix *X, Matrix *W, Matrix *H, Matrix *H_new, double alpha, double eps){
	int M = X->row, N = X->col, K= W->col, j, k;
	double h_tmp, diff1;
	double *WHT = (double *)malloc(M*N*sizeof(double));

	//mul_matrix(W,H,mul);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, W->data, K, H->data, K, 0.0, WHT, N);

	for (j = 0; j < N; ++j){
		for (k = 0; k < K; ++k){
			h_tmp = _mur_update_H_jk(X,W,H,WHT,j,k,alpha);
			H_new->data[j*K + k] = fmax(eps, h_tmp);
		}
	}

	free(WHT);
}

void _mur_update(Matrix *X, Matrix *W, Matrix *H, Matrix *upd_W, Matrix *upd_H,double alpha,double eps){
	_mur_update_W(X,W,H,upd_W,alpha,eps);
	_mur_update_H(X,upd_W,H,upd_H,alpha,eps);
}

int mur_solver(Matrix *X, Matrix *W_ret, Matrix *H_ret, int *flag,double alpha, 
	double eps, double delta1,double delta2,int seed, int init_max){
	int rounds, M = X->row, N = X->col, K = W_ret->col;
	double alpha_div_val, pre_alpha_div_val, alpha_div_init, alpha_div_ret;
	Matrix *W, *H, *W_new, *H_new;
	FILE *fp;
	clock_t t,base_t;

	// initial processing
	W = create_matrix(M,K); H = create_matrix(N,K);
	W_new = create_matrix(M,K); H_new = create_matrix(N,K);
	init(W,H,eps,seed,init_max);
	
	// start message
	alpha_div_init = alpha_div(X,W,H,alpha);
	if (flag[0]) // MSG_OUTPUT
		printf("alpha div init: %f\n", alpha_div_init);
	
	pre_alpha_div_val = alpha_div_init;

	// log setup
	if (LOG_FLAG){
		fp = fopen("log_mur.csv","w");
		base_t = clock();
	}

	for (rounds = 1; rounds <= IT_LIMITS; ++rounds){
		_mur_update(X,W,H,W_new,H_new,alpha,eps);
		
		cblas_dcopy(W->row * W->col, W_new->data, 1, W->data, 1);
		cblas_dcopy(H->row * H->col, H_new->data, 1, H->data, 1);
		

		if (flag[0] || flag[1] || LOG_FLAG)
			alpha_div_val = alpha_div(X,W,H,alpha);

		// log
		if (LOG_FLAG){
			t = clock();
			fprintf(fp,"%d,%.5f,%.5f\n",rounds,
					(double)(t-base_t)/CLOCKS_PER_SEC,alpha_div_val);
		}

		// message
		if (flag[0])
			printf("rounds %d -- alpha div value: %f\n", rounds, alpha_div_val);

		// error handling
		if (flag[1]){
			if (pre_alpha_div_val < alpha_div_val){
				printf("MUsolver: global convergence error");
				exit(EXIT_FAILURE);
			}

			pre_alpha_div_val = alpha_div_val;
		}
		
		if (stop(X,W,H,alpha,eps,delta1,delta2))
			break;
	}
	
	// end message
	alpha_div_ret = alpha_div(X,W,H,alpha);
	if (flag[0])
		printf("alpha div result: %f\n", alpha_div_ret);

	copy_matrix(W_ret, W); 
	copy_matrix(H_ret, H);

	free_matrix(W); free_matrix(H);
	free_matrix(W_new); free_matrix(H_new);

	if (LOG_FLAG)
		fclose(fp);

	return rounds;
}

