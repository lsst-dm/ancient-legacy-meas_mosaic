#include <stdio.h>
#include <math.h>
#include <vector>
#include <ctime>
#include <strings.h>
#include "fitsio.h"

#include <gsl/gsl_linalg.h>
#include "lsst/meas/mosaic/mosaicfit.h"

using namespace lsst::meas::mosaic;

using namespace lsst::meas::mosaic;
using namespace lsst::afw::detection;

double* solveMatrix_GSL(long size, double *a_data, double *b_data) {
    std::cout << "size: " <<size << std::endl;
    std::cout << "a_data" << std::endl;
    for (int i=0; i < size; i++) 
        std::cout <<i << ": " << a_data[i] << std::endl;
    std::cout << "b_data" << std::endl;
    for (int i=0; i < size; i++) 
        std::cout <<i << ": " << b_data[i] << std::endl;

    gsl_matrix_view a = gsl_matrix_view_array(a_data, size, size);
    gsl_vector_view b = gsl_vector_view_array(b_data, size);
    double *c_data = new double[size];
    gsl_vector_view c = gsl_vector_view_array(c_data, size);

    int s;
    gsl_permutation *p = gsl_permutation_alloc(size);

    std::cout << "Calling decomp"   << std::endl;
    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    std::cout << "Calling solve"   << std::endl;
    gsl_linalg_LU_solve(&a.matrix, p, &b.vector, &c.vector);
    std::cout << "Calling free"   << std::endl;
    gsl_permutation_free(p);
    std::cout << "Returning"   << std::endl;

    return c_data;
}
int main(int argc, char **argv)
{
    char buf[BUFSIZ];
    int icnt;
    double var;
    FILE *fp = fopen("gslhang.dat", "rt");
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%d", &icnt);
    std::cout << icnt << std::endl;
    double *a_data = new double[icnt*icnt];
    double *b_data = new double[icnt];
    for (int i = 0; i < icnt*icnt; i++) {
	fgets(buf, BUFSIZ, fp);
	sscanf(buf, "%lf", &var);
        //#std::cout << var << std::endl;
        a_data[i] = var;
    }
    for (int i = 0; i < icnt; i++) {
	fgets(buf, BUFSIZ, fp);
	sscanf(buf, "%lf", &var);
        //#std::cout << var << std::endl;
        b_data[i] = var;
    }
    solveMatrix_GSL(icnt, a_data, b_data);
    fclose(fp);

    return 0;
}
