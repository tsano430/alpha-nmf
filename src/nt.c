/******************************************************************************
 * nt.c: nakatsu and takahashi algorithm
 ******************************************************************************/

#include "nt.h"

#define LOG_FLAG 0

double _nt_update_W_ik(Matrix *X,Matrix *W,Matrix *H,int i,int k,
		double alpha,double eps,double delta1,double delta2){
	int M = X->row, K = W->col;
	double diff1,diff2,diff1_n,W_n_ik,W_ik = W->data[i*K+k];
	Matrix *W_n;

	// step 1
	diff1 = alpha_div_diff1_W(X,W,H,i,k,alpha);

	if (W_ik <= eps+delta2 && W_ik >= eps && diff1 >= -delta1)
		return W_ik;
	else if (W_ik > eps+delta2 && fabs(diff1) <= delta1)
		return W_ik;

	diff2 = alpha_div_diff2_W(X,W,H,i,k,alpha);

	W_n_ik = W_ik - diff1/diff2;

	// step 2
	W_n_ik = fmax(eps,W_n_ik);

	// step 3
	W_n = create_matrix(M,K); copy_matrix(W_n,W);
	W_n->data[i*K+k] = W_n_ik;
	diff1_n = alpha_div_diff1_W(X,W_n,H,i,k,alpha);
	free_matrix(W_n);

	if (diff1_n * diff1 >= 0.0)
		return W_n_ik;
	else
		return W_ik - (W_n_ik - W_ik) / (diff1_n - diff1) * diff1;
}

double _nt_update_H_jk(Matrix *X,Matrix *W,Matrix *H,int j,int k,
		double alpha,double eps,double delta1,double delta2){
	int N = X->col, K = W->col;
	double diff1,diff2,diff1_n,H_n_jk,H_jk = H->data[j*K+k];
	Matrix *H_n;

	// step 1
	diff1 = alpha_div_diff1_H(X,W,H,j,k,alpha);

	if (H_jk <= eps+delta2 && H_jk >= eps && diff1 >= -delta1)
		return H_jk;
	else if (H_jk > eps+delta2 && fabs(diff1) <= delta1)
		return H_jk;

	diff2 = alpha_div_diff2_H(X,W,H,j,k,alpha);

	H_n_jk = H_jk - diff1/diff2;

	// step 2
	H_n_jk = fmax(eps,H_n_jk);

	// step 3
	H_n = create_matrix(N,K); copy_matrix(H_n,H);
	H_n->data[j*K+k] = H_n_jk;
	diff1_n = alpha_div_diff1_H(X,W,H_n,j,k,alpha);
	free_matrix(H_n);

	if (diff1_n * diff1 >= 0.0)
		return H_n_jk;
	else
		return H_jk - (H_n_jk - H_jk) / (diff1_n - diff1) * diff1;
}

void _nt_update_W(Matrix *X, Matrix *W, Matrix *H, Matrix *W_new,
		int *flag,double alpha,double eps,double delta1,double delta2){
	int M = X->row, N = X->col, K = W->col, i, k;
	double w_tmp, diff1;
	double bef,aft;

	copy_matrix(W_new, W);

	for (i = 0; i < M; ++i){
		for (k = 0; k < K; ++k){
			if (flag[1])
				bef = alpha_div(X,W_new,H,alpha);

			// update_W_ik
			w_tmp = _nt_update_W_ik(X,W_new,H,i,k,alpha,eps,delta1,delta2);
			W_new->data[i*K + k] = fmax(eps,w_tmp);

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

void _nt_update_H(Matrix *X, Matrix *W, Matrix *H, Matrix *H_new,
		int *flag,double alpha,double eps,double delta1,double delta2){
	int M = X->row, N = X->col, K = W->col, j, k;
	double h_tmp, diff1;
	double bef,aft;

	copy_matrix(H_new, H);

	for (j = 0; j < N; ++j){
		for (k = 0; k < K; ++k){
			if (flag[1])
				bef = alpha_div(X,W,H_new,alpha);

			// update_H_jk
			h_tmp = _nt_update_H_jk(X,W,H_new,j,k,alpha,eps,delta1,delta2);
			H_new->data[j*K + k] = fmax(eps,h_tmp);

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

//void update_W(Matrix *X, Matrix *W, Matrix *H, Matrix *W_new){
//	int M = X->row, N = X->col, i, k;
//	double w_tmp, diff1;
//	double bef,aft;
//
//       	copy_matrix(W_new, W);
//
//	for (i = 0; i < M; ++i){
//		for (k = 0; k < K; ++k){
//			diff1 = alpha_div_diff1_W(X,W_new,H,i,k);
//
//			// unstable value skipped
//			if (fabs(W_new->data[i*K + k] - eps) < TOLERANT_ERROR &&
//				diff1 >= -delta3)
//				;
//			else if (W_new->data[i*K + k] > eps && fabs(diff1) <= delta3)
//				;
//			else{
//				if (CONVERGENCE_DEBUG)
//					bef = alpha_div(X,W_new,H);
//
//				// update_W_ik
//				double before_x = W_new->data[i*K+k];
//				w_tmp = fmax(eps,update_W_ik(X,W_new,H,i,k));
//				W_new->data[i*K + k] = w_tmp;
//
//				if (CONVERGENCE_DEBUG)
//					aft = alpha_div(X,W_new,H);
//
//				// global convergence debug
//				if (CONVERGENCE_DEBUG){
//					double d_aft = alpha_div_diff1_W(X,W_new,H,i,k);
//					if (fabs(bef-aft) < TOLERANT_ERROR || bef > aft)
//						;
//					else if (bef < aft){
//						printf("update_W: W[%d,%d] global convergence error\n",i,k);
//						printf("before_x, after_x: %.25f,%.25f\n", before_x,w_tmp);
//						printf("before_val, after_val: %.25f,%.25f\n",bef,aft);
//						printf("before_diff, after_diff: %.25f,%.25f\n",diff1,d_aft);
//						save_matrix("W_debug.csv", W_new);
//						save_matrix("H_debug.csv", H);
//						exit(EXIT_FAILURE);
//					}
//				}
//			}
//		}
//	}
//}
//
//
//void update_H(Matrix *X, Matrix *W, Matrix *H, Matrix *H_new){
//	int M = X->row, N = X->col, j, k;
//	double h_tmp, diff1;
//	double bef,aft;
//
//	copy_matrix(H_new, H);
//
//	for (j = 0; j < N; ++j){
//		for (k = 0; k < K; ++k){
//			diff1 = alpha_div_diff1_H(X,W,H_new,j,k);
//			
//			// unstable value skipped
//			if (fabs(H_new->data[j*K + k] - eps) < TOLERANT_ERROR &&
//				diff1 >= -delta3)
//				;
//			else if (H_new->data[j*K + k] > eps && fabs(diff1) <= delta3)
//				;
//			else{
//				if (CONVERGENCE_DEBUG)
//					bef = alpha_div(X,W,H_new);
//
//				// update_H_jk
//				h_tmp = fmax(eps,update_H_jk(X,W,H_new,j,k));
//				H_new->data[j*K + k] = h_tmp;
//
//				if (CONVERGENCE_DEBUG)
//					aft = alpha_div(X,W,H_new);
//
//				// global convergence debug
//				if (CONVERGENCE_DEBUG){
//					double d_aft = alpha_div_diff1_H(X,W,H_new,j,k);
//
//					if (fabs(bef-aft) < TOLERANT_ERROR || aft < bef)
//						;
//					else if (aft > bef){
//						printf("update_H: H[%d,%d] global convergence error\n",j,k);
//						printf("before_val,after_val: %.25f,%.25f\n",bef,aft);
//						printf("before_diff, after_diff: %.25f,%.25f\n",diff1,d_aft);
//						exit(EXIT_FAILURE);
//					}
//				}
//			}
//		}
//	}
//}

void _nt_update(Matrix *X, Matrix *W, Matrix *H, Matrix *W_new, Matrix *H_new,
		int *flag,double alpha,double eps,double delta1,double delta2){
	_nt_update_W(X,W,H,W_new,flag,alpha,eps,delta1,delta2);
	_nt_update_H(X,W_new,H,H_new,flag,alpha,eps,delta1,delta2);
}

int nt_solver(Matrix *X, Matrix *W_ret, Matrix *H_ret,int *flag,double alpha,
		double eps,double delta1,double delta2,int seed,int init_max){
	int rounds, M = X->row, N = X->col, K = W_ret->col;
	double alpha_div_val, pre_alpha_div_val, alpha_div_init, alpha_div_ret;
	Matrix *W, *H, *W_new, *H_new;
	FILE *fp;
	clock_t base_t,t;

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
	if (LOG_FLAG){
		fp = fopen("log_n1.csv","w");
		base_t = clock();
	}

	// solve
	for (rounds = 1; rounds <= IT_LIMITS; ++rounds){
		_nt_update(X,W,H,W_new,H_new,flag,alpha,eps,delta1,delta2);
		
		copy_matrix(W, W_new); copy_matrix(H, H_new);

		if (flag[0] || flag[1] || LOG_FLAG)
			alpha_div_val = alpha_div(X,W,H,alpha);
		
		// log
		if (LOG_FLAG){
			t = clock();
			fprintf(fp,"%d,%.5f,%.5f\n",rounds,
			(double)(t-base_t)/CLOCKS_PER_SEC, alpha_div_val);
		}

		// message
		if (flag[0])
			printf("rounds %d -- alpha div value: %f\n", rounds, alpha_div_val);

		// global convergence debug
		if (flag[1]){
			if (pre_alpha_div_val < alpha_div_val){
				printf("n1_solver: entire global convergence error");
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


