/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
/* #include "lib_poisson1D.h" */
#include "../include/lib_poisson1D.h"
#include "../include/atlas_headers.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
    for(int i = 0; i < (*la); i++){
        for(int j = 0; j < (*lab); j++){
            if (j < (*kv))AB[i*(*lab)+j]=0;
            else if(j==(*kv)+1)AB[i*(*lab)+j]=2; // D
            else AB[i*(*lab)+j]=-1; // upper & lower
        }
    }
    AB[(*kv)]=0;
    AB[(*lab)*(*la)-1]=0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1) {
    RHS[0] = (*BC0);
    for (int i = 1; i < (*la); i++) {
        RHS[i] = 0.0;
    }
    RHS[(*la)-1] = (*BC1);
}


void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
    for (int i = 0; i < (*la); i++) {
        EX_SOL[i]=(*BC0)+(X[i]*((*BC1)-(*BC0)));
    }
}

void set_grid_points_1D(double* x, int* la) {
    double h = 1.0 / (*la + 1);
    x[0] = h;
    for (int i = 1; i < (*la); i++) {
        x[i] = x[i-1] + h;
    }
}


void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
    FILE * file;
    int ii,jj;
    file = fopen(filename, "w");
    //Numbering from 1 to la
    if (file != NULL){
        for (ii=0;ii<(*lab);ii++){
            for (jj=0;jj<(*la);jj++){
                fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
            }
            fprintf(file,"\n");
        }
        fclose(file);
    }
    else{
        perror(filename);
    }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
    FILE * file;
    int ii,jj;
    file = fopen(filename, "w");
    //Numbering from 1 to la
    if (file != NULL){
        for (ii=0;ii<(*la);ii++){
            for (jj=0;jj<(*lab);jj++){
                fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
            }
            fprintf(file,"\n");
        }
        fclose(file);
    }
    else{
        perror(filename);
    }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
    FILE * file;
    int jj;
    file = fopen(filename, "w");
    //Numbering from 1 to la
    if (file != NULL){
        for (jj=1;jj<(*la);jj++){
            fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
        }
        for (jj=0;jj<(*la);jj++){
            fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
        }
        for (jj=0;jj<(*la)-1;jj++){
            fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
        }
        fclose(file);
    }
    else{
        perror(filename);
    }
}

void write_vec(double* vec, int* la, char* filename){
    int jj;
    FILE * file;
    file = fopen(filename, "w");
    // Numbering from 1 to la
    if (file != NULL){
        for (jj=0;jj<(*la);jj++){
            fprintf(file,"%lf\n",vec[jj]);
        }
        fclose(file);
    }
    else{
        perror(filename);
    } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
    int jj;
    FILE * file;
    file = fopen(filename, "w");
    // Numbering from 1 to la
    if (file != NULL){
        for (jj=0;jj<(*la);jj++){
            fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
        }
        fclose(file);
    }
    else{
        perror(filename);
    } 
}  

int indexABCol(int i, int j, int *lab){
    return j*(*lab)+i;
}
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
    ipiv[0]=1;
    for(int i = 1; i < *la; i++){
        if(AB[ *lab*i-2 ]==0){
            *info=1;
            return *info;
        }
        AB[(*lab*i)-1] /= AB[*lab*i-2]; // b(i-1) /= a(i-1)
        AB[*lab*(i+1)-2] -= AB[*lab*i-1] * AB[*lab*(i+1)-3]; // a(i) -= b(i-1) * c(i-1)
        ipiv[i]=i+1;
    }
    return *info;
}

void eig_poisson1D(double* eigval, int *la){
}

double eigmax_poisson1D(int *la){ // lambda_max
    double h = 1.0 / ((*la) +1.0);
    return 4*sin(((*la) * M_PI * h)/2)*sin(((*la) * M_PI * h)/2);
}

double eigmin_poisson1D(int *la){ // lambda_min
    double h = 1.0 / ((*la) +1.0);
    return 4*sin((M_PI * h)/2)*sin(( M_PI * h)/2);
}

double richardson_alpha_opt(int *la){
    // return 0.5 | Démonstration dans le rapport
    return 2/(eigmax_poisson1D(la)+eigmin_poisson1D(la)); 
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    double *B = malloc(*la* sizeof(double));
    double norme_B = cblas_dnrm2(*la,RHS,1);

    for((*nbite) = 0; (*nbite) < *maxit ;(*nbite)++){
        //Copie de RHS dans B
        cblas_dcopy(*la,  RHS, 1, B,1);
        //b = b - Ax
        cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,B,1);
        //Calcul du residu
        resvec[*nbite] = cblas_dnrm2(*la,B,1) / norme_B;
        //x = x + alpha * b | x = x + alpha*(b-Ax)
        cblas_daxpy(*la,*alpha_rich,B,1,X,1);
        if(resvec[*nbite]<=*tol) break;
    }
    free(B);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    //M = D
    for (int i = 0; i < *la; i++)
        MB[i*(*lab)+1] = AB[i*(*lab)+1];
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    //M = D - E 
    for(int i = 0; i < *la; i++){
        MB[*lab*i+1] = AB[*lab*i+1]; // M = D
        MB[*lab*i+2] = AB[*lab*i+2]; // M = M - E
    }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    double *B = malloc(*la * sizeof(double));
    double norme_B = cblas_dnrm2(*la,RHS,1);
    int * ipiv = malloc(*la * sizeof(int));
    int info = 0;
    int NRHS = 1;
    int ku_minus = *ku-1; // nécessaire pour MB

    // LU factorization of MB
    dgbtrf_(la, la, kl, &ku_minus, MB, lab, ipiv, &info);
    for((*nbite) = 0; (*nbite) < *maxit ;(*nbite)++){
        //Copie de RHS pour garder le même RHS à chaque iteration
        cblas_dcopy(*la,  RHS, 1, B,1);
        //b = b - Ax
        cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,B,1);
        //Calcul du residu
        resvec[*nbite] = cblas_dnrm2(*la,B,1) / norme_B;
        //b = b/M | b = (b - Ax)/M
        dgbtrs_("N", la, kl, &ku_minus, &NRHS, MB, lab, ipiv, B, la, &info,1);
        // x = x + b | x = x + (b - Ax)/M
        cblas_daxpy(*la,1,B,1,X,1);
        if(resvec[*nbite]<=*tol) break;
    }
    free(B);
    free(ipiv);
}

//Exercise 10 CSR and CSC (nothing implemented more than that)
void set_CR_operator_poisson1D(double* AA, double* JA, double* IA,  int *nnz, int *la){
    int j=0;
    IA[0]=0;
    for(int i = 0; i < *nnz; i++){
        if(i%3==0)AA[i]=2;
        else AA[i]=-1;
        if(i%3==2)JA[i]=j-1;
        else if(i%3==0){
            JA[i]=j;
            j++;
        }
        else JA[i]=j;
    }
    IA[*la-1]=*nnz-1;
    IA[*la]=*nnz;
}

