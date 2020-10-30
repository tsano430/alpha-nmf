/******************************************************************************
 * simulation.c
 ******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <time.h>
#include "utils.h"
#include "mur.h"
#include "nt.h"
#include "proposed.h"

int main(void){
	char *fname[2] = {"../input/orl.csv", "../input/cluto_tr23.csv"};
	int flag[2] = {1,0}, init_max = 10,rounds, M, N, K, i, j;
	double eps, delta1, delta2, alpha;
	clock_t bef_time, aft_time;
	Matrix *X,*W,*H;
    double ret_orl[3][26], ret_cluto[3][26];
	int seed, trial = 5;
    FILE *fp;

    for (i = 0; i < 3; ++i){
        for (j = 0; j <= 25; ++j){
            ret_orl[i][j] = 0.0;
            ret_cluto[i][j] = 0.0;
        }
    }

    /* ORL dataset */
    M = 2576; 
    N = 400; 
    K = 5;
    eps = 1.0; 
    delta1 = 10.0; 
    delta2 = 1.0;

    for (i = 0; i <= 25; ++i){
        alpha = i * 0.2 - 2.5;

        X = create_matrix(M,N);
        W = create_matrix(M,K); 
        H = create_matrix(N,K);
        X = load_matrix(fname[0],alpha,eps);

        for (seed = 0; seed < trial; ++seed){
            // MUR
            bef_time = clock();
            rounds = mur_solver(X,W,H,flag,alpha,eps,delta1,delta2,seed,init_max);
            aft_time = clock();
            ret_orl[0][i] += (aft_time - bef_time) / CLOCKS_PER_SEC;
            // NT
            bef_time = clock();
            rounds = nt_solver(X,W,H,flag,alpha,eps,delta1,delta2,seed,init_max);
            aft_time = clock();
            ret_orl[1][i] += (aft_time - bef_time) / CLOCKS_PER_SEC;
            // Proposed
            bef_time = clock();
            rounds = proposed_solver(X,W,H,flag,alpha,eps,delta1,delta2,seed,init_max);
            aft_time = clock();
            ret_orl[2][i] += (aft_time - bef_time) / CLOCKS_PER_SEC;
        }
        ret_orl[0][i] /= trial;
        ret_orl[1][i] /= trial;
        ret_orl[2][i] /= trial;
        printf("%f\n", ret_orl[2][i]);
    }
    free_matrix(X); free_matrix(W); free_matrix(H);

    fp = fopen("./log/ret_orl.log", "w");
    for (i = 0; i <= 25; ++i){
        alpha = i * 0.2 - 2.5;
        if (i == 25)
                fprintf(fp, "%.1f\n", alpha);
            else
                fprintf(fp, "%.1f,", alpha);
    }
    for (i = 0; i < 3; ++i){
        for (j = 0; j <= 25; ++j){
            if (j == 25)
                fprintf(fp, "%.6f\n", ret_orl[i][j]);
            else
                fprintf(fp, "%.6f,", ret_orl[i][j]);
        }
    }
    fclose(fp);

    /* CLUTO (tr23) dataset */
    M = 5832; 
    N = 204; 
    K = 6;
    eps = 0.1; 
    delta1 = 1.0; 
    delta2 = 0.1;

    for (i = 0; i <= 25; ++i){
        alpha = i * 0.2 - 2.5;

        X = create_matrix(M,N);
        W = create_matrix(M,K); 
        H = create_matrix(N,K);
        X = load_matrix(fname[1],alpha,eps);

        for (seed = 0; seed < trial; ++seed){
            // MUR
            bef_time = clock();
            rounds = mur_solver(X,W,H,flag,alpha,eps,delta1,delta2,seed,init_max);
            aft_time = clock();
            ret_cluto[0][i] += (aft_time - bef_time) / CLOCKS_PER_SEC;
            // NT
            bef_time = clock();
            rounds = nt_solver(X,W,H,flag,alpha,eps,delta1,delta2,seed,init_max);
            aft_time = clock();
            ret_cluto[1][i] += (aft_time - bef_time) / CLOCKS_PER_SEC;
            // Proposed
            bef_time = clock();
            rounds = proposed_solver(X,W,H,flag,alpha,eps,delta1,delta2,seed,init_max);
            aft_time = clock();
            ret_cluto[2][i] += (aft_time - bef_time) / CLOCKS_PER_SEC;
        }
        ret_cluto[0][i] /= trial;
        ret_cluto[1][i] /= trial;
        ret_cluto[2][i] /= trial;
	}
    free_matrix(X); free_matrix(W); free_matrix(H);

    fp = fopen("./log/ret_cluto.log", "w");
    for (i = 0; i <= 25; ++i){
        alpha = i * 0.2 - 2.5;
        if (i == 25)
                fprintf(fp, "%.1f\n", alpha);
            else
                fprintf(fp, "%.1f,", alpha);
    }
    for (i = 0; i < 3; ++i){
        for (j = 0; j <= 25; ++j){
            if (j == 25)
                fprintf(fp, "%.6f\n", ret_cluto[i][j]);
            else
                fprintf(fp, "%.6f,", ret_cluto[i][j]);
        }
    }
    fclose(fp);

	return 0;
}
