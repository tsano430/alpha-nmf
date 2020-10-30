/******************************************************************************
 * mono.c -- Proposed Algorithm
 ******************************************************************************/

#include "proposed.h"

#define LOG_FLAG 0

const double beta = 0.95;

double _proposed_update_W_ik(Matrix *X,Matrix *W,Matrix *H,int i,int k,
		double alpha,double eps,double delta1,double delta2){
	int M = X->row, K = W->col, m = 1;
	double diff1,diff2,diff1_n,dy;
	double W_n_ik, W_ik, W_i_ik, d, y;
	Matrix *W_n;

	// init
	W_ik = W->data[i*K+k];
	diff1 = alpha_div_diff1_W(X,W,H,i,k,alpha);
	W_n = create_matrix(M,K);
	cblas_dcopy(W->row * W->col,W->data,1,W_n->data,1);	

	// stop condition for (i,k) element
	if (diff1 > delta1 && W_ik < eps+delta2){
		free_matrix(W_n);
		return W_ik;
	}
	else if (fabs(diff1) <= delta1){
		free_matrix(W_n);
		return W_ik;
	}

	//printf("diff1: %f\n",diff1);

	diff2 = alpha_div_diff2_W(X, W_n, H, i, k, alpha);

	d = -diff1 / diff2; W_n_ik = W_ik + d;

	// linesearch to construct separating hyperplane
	W_n->data[i * K + k] = W_n_ik;
	diff1_n = alpha_div_diff1_W(X, W_n, H, i, k, alpha);

	if (diff1_n * d <= 0.0){
		y = W_n_ik;
	}else{
		while (1){
			y = W_ik + pow(beta,m) * d;
			W_n->data[i * K + k] = y;
			dy = alpha_div_diff1_W(X, W_n, H, i, k, alpha);

			//printf("dy: %f\n",dy);

			if (dy * d <= 0.0)
				break;
			m++;
		}
	}

	W_i_ik = W_ik - (W_n_ik - W_ik) / (diff1_n - diff1) * diff1;

	if (d > 0.0 && m >= 1 && y < W_i_ik)
		W_ik = W_i_ik;
	else if (d < 0.0 && m >= 1 && y > W_i_ik)
		W_ik = W_i_ik;
	else
		W_ik = y;
	
	// update W_ik
	W_ik = fmax(eps,W_ik);

	free_matrix(W_n);
	return W_ik;
}

double _proposed_update_H_jk(Matrix *X,Matrix *W,Matrix *H,int j,int k,
		double alpha,double eps,double delta1,double delta2){
	int N = X->col, K = W->col, m = 1;
	double diff1,diff2,diff1_n,dy;
	double H_n_jk,H_jk,H_i_jk,d,y;
	Matrix *H_n;

	// init
	H_jk = H->data[j*K+k];
	diff1 = alpha_div_diff1_H(X,W,H,j,k,alpha);
	H_n = create_matrix(N, K);
	cblas_dcopy(H->row * H->col,H->data,1,H_n->data,1);

	if (diff1 > delta1 && H_jk < eps+delta2){
		free_matrix(H_n);
		return H_jk;
	}
	else if (fabs(diff1) <= delta1){
		free_matrix(H_n);
		return H_jk;
	}

	diff2 = alpha_div_diff2_H(X, W, H_n, j, k, alpha);

	d = -diff1/diff2; H_n_jk = H_jk + d;

	// linesearch
	H_n->data[j * K + k] = H_n_jk;
	diff1_n = alpha_div_diff1_H(X, W, H_n, j, k, alpha);
	
	if (diff1_n * d <= 0.0){
		y = H_n_jk;
	}else{
		while (1){
			y = H_jk + pow(beta,m) * d;
			H_n->data[j * K + k] = y;
			dy = alpha_div_diff1_H(X, W, H_n, j, k, alpha);

			if (dy * d <= 0.0)
				break;
			m++;
		}
	}

	H_i_jk = H_jk - (H_n_jk - H_jk) / (diff1_n - diff1) * diff1;

	if (d > 0.0 && m >= 1 && y < H_i_jk)
		H_jk = H_i_jk;
	else if (d < 0.0 && m >= 1 && y > H_i_jk)
		H_jk = H_i_jk;
	else
		H_jk = y;
	
	// update H_jk
	H_jk = fmax(eps,H_jk);

	free_matrix(H_n);

	return H_jk;
}

void _proposed_update_W(Matrix *X, Matrix *W, Matrix *H, Matrix *W_new,int *flag,
		double alpha,double eps,double delta1,double delta2){
	int M = X->row, N = X->col, K = W->col, i, k;
	double diff1;
	double bef,aft;

	cblas_dcopy(W->row * W->col, W->data, 1, W_new->data, 1);

	for (i = 0; i < M; ++i){
		for (k = 0; k < K; ++k){
			if (flag[1])
				bef = alpha_div(X,W_new,H,alpha);

			// update_W_ik
			W_new->data[i*K + k] = _proposed_update_W_ik(X,W_new,H,i,k,alpha,eps,delta1,delta2);

			if (flag[1])
				aft = alpha_div(X,W_new,H,alpha);

			// global convergence debug
			if (flag[1]){
				if (fabs(bef-aft) < 1e-12 || bef > aft)
					;
				else if (bef < aft){
					printf("update_W: W[%d,%d] global convergence error\n",i,k);
					exit(EXIT_FAILURE);
				}
			}
		}
	}
}

void _proposed_update_H(Matrix *X, Matrix *W, Matrix *H, Matrix *H_new,int *flag,
		double alpha,double eps,double delta1,double delta2){
	int M = X->row, N = X->col, K = W->col, j, k;
	double diff1;
	double bef,aft;

	cblas_dcopy(H->row * H->col, H->data, 1, H_new->data, 1);

	for (j = 0; j < N; ++j){
		for (k = 0; k < K; ++k){
			if (flag[1])
				bef = alpha_div(X,W,H_new,alpha);

			// update_H_jk
			H_new->data[j*K + k] = _proposed_update_H_jk(X,W,H_new,j,k,alpha,eps,delta1,delta2);

			if (flag[1])
				aft = alpha_div(X,W,H_new,alpha);

			// global convergence debug
			if (flag[1]){
				if (fabs(bef-aft) < 1e-12 || aft < bef)
					;
				else if (aft > bef){
					printf("update_H: H[%d,%d] global convergence error\n",j,k);
					exit(EXIT_FAILURE);
				}
			}
		}
	}
}

void _proposed_update(Matrix *X, Matrix *W, Matrix *H, Matrix *W_new, Matrix *H_new,
		int *flag,double alpha,double eps,double delta1,double delta2){
	_proposed_update_W(X,W,H,W_new,flag,alpha,eps,delta1,delta2);
	_proposed_update_H(X,W_new,H,H_new,flag,alpha,eps,delta1,delta2);
}

int proposed_solver(Matrix *X, Matrix *W_ret, Matrix *H_ret,int *flag,double alpha,
		double eps,double delta1,double delta2,int seed,int init_max){
	int rounds, M = X->row, N = X->col, K = W_ret->col;
	double alpha_div_val, pre_alpha_div_val, alpha_div_init, alpha_div_ret;
	Matrix *W, *H, *W_new, *H_new;
	FILE *fp;

	// initial processing
	W = create_matrix(M,K); H = create_matrix(N,K);
	W_new = create_matrix(M,K); H_new = create_matrix(N,K);
	init(W,H,eps,seed,init_max);
	
	// start message
	alpha_div_init = alpha_div(X,W,H,alpha);
	if (flag[0]) 
		printf("alpha div init: %f\n", alpha_div_init);
	pre_alpha_div_val = alpha_div_init;

	// log setup
	if (LOG_FLAG)
		fp = fopen("logdata.csv","w");

	// solve
	for (rounds = 1; rounds <= IT_LIMITS; ++rounds){
		_proposed_update(X,W,H,W_new,H_new,flag,alpha,eps,delta1,delta2);
		
		cblas_dcopy(W->row * W->col, W_new->data, 1, W->data, 1); 
		cblas_dcopy(H->row * H->col, H_new->data, 1, H->data, 1);

		if (flag[0] || flag[1] || LOG_FLAG)
			alpha_div_val = alpha_div(X,W,H,alpha);
		
		// log
		if (LOG_FLAG)
			fprintf(fp,"%d,%.5f\n",rounds, alpha_div_val);

		// message
		if (flag[0])
			printf("rounds %d -- alpha div value: %f\n", rounds, alpha_div_val);

		// global convergence debug
		if (flag[1]){
			if (pre_alpha_div_val < alpha_div_val){
				printf("mono_solver: entire global convergence error");
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