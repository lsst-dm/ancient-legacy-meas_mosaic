#include <stdio.h>
#include <math.h>
#include <vector>
#include <ctime>
#include <cstring>
#include <iostream>
#include <strings.h>
#include "fitsio.h"

#include <mkl_lapack.h>
double* solveMatrix_MKL(long size, double *a_data, double *b_data) {
    //char L = 'L';
    MKL_INT n = size;
    MKL_INT nrhs = 1;
    MKL_INT lda = size;
    MKL_INT *ipiv = new MKL_INT[size];
    MKL_INT ldb = size;
    MKL_INT info = 0;

    //double *a = new double[size*size];
    //double *b = new double[size];

    //memcpy(a, a_data, sizeof(double)*size*size);
    //memcpy(b, b_data, sizeof(double)*size);

    //dgesv(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    dgesv(&n, &nrhs, a_data, &lda, ipiv, b_data, &ldb, &info);
    //dposv(&L, &n, &nrhs, a, &lda, b, &ldb, &info);

    double *c_data = new double[size];
    //memcpy(c_data, b, sizeof(double)*size);
    memcpy(c_data, b_data, sizeof(double)*size);
    delete [] ipiv;
    //delete [] a;
    //delete [] b;

    return c_data;
}

int main(int argc, char **argv)
{
    char buf[BUFSIZ];
    char * str;
    int icnt;
    double var;
    FILE *fp = fopen("gslhang.dat", "rt");
    str = fgets(buf, BUFSIZ, fp);
    sscanf(str, "%d", &icnt);
    double *a_data = new double[icnt*icnt];
    double *b_data = new double[icnt];
    for (int i = 0; i < icnt*icnt; i++) {
	str = fgets(buf, BUFSIZ, fp);
	sscanf(str, "%lf", &var);
        a_data[i] = var;
    }
    for (int i = 0; i < icnt; i++) {
	str = fgets(buf, BUFSIZ, fp);
	sscanf(str, "%lf", &var);
        b_data[i] = var;
    }
    std::cout << "calling solveMatrix_MKL" << std::endl;
    solveMatrix_MKL(icnt, a_data, b_data);
    std::cout << "solveMatrix_MKL ran to completion" << std::endl;
    fclose(fp);

    return 0;
}
