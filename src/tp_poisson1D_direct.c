/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "../include/lib_poisson1D.h"
#include "../include/atlas_headers.h"

int main(int argc,char *argv[])
    /* ** argc: Nombre d'arguments */
    /* ** argv: Valeur des arguments */
{
    int ierr;
    int jj;
    int nbpoints, la;
    int ku, kl, kv, lab;
    int *ipiv;
    int info;
    int NRHS;
    double T0, T1;
    double *RHS, *EX_SOL, *X;
    double **AAB;
    double *AB;

    double temp, relres;

    double norm_exsol;
    double norm_sol;
    double norm_res;

    struct timespec start, end;
    double cpu_time_used;

    NRHS=1;
    nbpoints=10;
    la=nbpoints-2;
    T0=-5.0;
    T1=5.0;

    printf("--------- Poisson 1D ---------\n\n");
    RHS=(double *) malloc(sizeof(double)*la);
    EX_SOL=(double *) malloc(sizeof(double)*la);
    X=(double *) malloc(sizeof(double)*la);

    // TODO : you have to implement those functions
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    write_vec(RHS, &la, "RHS.dat");
    write_vec(EX_SOL, &la, "EX_SOL.dat");
    write_vec(X, &la, "X_grid.dat");

    kv=1;
    ku=1;
    kl=1;
    lab=kv+kl+ku+1;

    AB = (double *) malloc(sizeof(double)*lab*la);

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);


    /* write_vec(RHS, &la, "RHS.dat"); */

    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
    

    printf("Solution with LAPACK\n");
    //----------------------TEST ON DGBMV----------------------------
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,kl,ku,1.0,AB+1,lab,EX_SOL,1,0.0,RHS,1);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    cpu_time_used = ( (double) end.tv_sec + (double) end.tv_nsec/1000000000) - ((double) start.tv_sec + (double) start.tv_nsec/1000000000);

    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    norm_sol = cblas_dnrm2(la, RHS, 1);
    norm_res = cblas_dnrm2(la, EX_SOL, 1);
    relres = norm_res / (norm_exsol * norm_sol);
    printf("DGBMV :\n\tTemps: %f\n\tRelres:%e\n",cpu_time_used,relres);
    /* LU Factorization */
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    info=0;
    ipiv = (int *) calloc(la, sizeof(int));
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    /* LU for tridiagonal matrix  (can replace dgbtrf_) */
    /* ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info); */

    /* write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat"); */

    /* Solution (Triangular) */
    if (info==0){
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info,1);
        if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
        printf("\n INFO = %d\n",info);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    cpu_time_used = ( (double) end.tv_sec + (double) end.tv_nsec/1000000000) - ((double) start.tv_sec + (double)start.tv_nsec/1000000000);

    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    norm_sol = cblas_dnrm2(la, RHS, 1);
    norm_res = cblas_dnrm2(la, EX_SOL, 1);
    relres = norm_res / (norm_exsol * norm_sol);
    printf("DGBTRF&DGBTRS :\n\tTemps: %f\n\tRelres:%e\n",cpu_time_used,relres);

    //----------------DGBTRFTRIDIAG---------------------
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    info=0;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    /* LU for tridiagonal matrix  (can replace dgbtrf_) */
    ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    /* write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat"); */

    /* Solution (Triangular) */
    if (info==0){
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info,1);
        if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
        printf("\n INFO = %d\n",info);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    cpu_time_used = ( (double) end.tv_sec + (double) end.tv_nsec/1000000000) - ((double) start.tv_sec + (double)start.tv_nsec/1000000000);

    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    norm_sol = cblas_dnrm2(la, RHS, 1);
    norm_res = cblas_dnrm2(la, EX_SOL, 1);
    relres = norm_res / (norm_exsol * norm_sol);
    printf("Tridiag&DGBTRS :\n\tTemps: %f\n\tRelres:%e\n",cpu_time_used,relres);
    /* It can also be solved with dgbsv */
    // TODO : use dgbsv
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    cpu_time_used = ( (double) end.tv_sec + (double) end.tv_nsec/1000000000) - ((double) start.tv_sec + (double)start.tv_nsec/1000000000);

    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    norm_sol = cblas_dnrm2(la, RHS, 1);
    norm_res = cblas_dnrm2(la, EX_SOL, 1);
    relres = norm_res / (norm_exsol * norm_sol);
    printf("DGBSV :\n\tTemps: %f\n\tRelres:%e\n",cpu_time_used,relres);
    if (info!=0){printf("\n INFO DGBSV = %d\n",info);}

    write_xy(RHS, X, &la, "SOL.dat");

    /* Relative forward error */
    // TODO : Compute relative norm of the residual

    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    norm_sol = cblas_dnrm2(la, RHS, 1);
    norm_res = cblas_dnrm2(la, EX_SOL, 1);
    relres = norm_res / (norm_exsol * norm_sol);

    printf("\nThe relative forward error is relres = %e\n",relres);

    free(RHS);
    free(EX_SOL);
    free(X);
    free(AB);
    printf("\n\n--------- End -----------\n");
}
