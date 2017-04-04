#ifndef CUSTOM_H_INCLUDED
#define CUSTOM_H_INCLUDED

#include <gsl/gsl_matrix.h>

#define COLUMNWISE 0
#define ROWWISE    1

void custom_vectorToMatrix(gsl_matrix *m, double y[], int rows, int columns, int shift);
void custom_matrixToVector(double y[], gsl_matrix *m, int rows, int columns, int shift);
void custom_vectorToMatrix_rc(gsl_matrix *m, double y[], int rows, int columns, int shift, int type);
void custom_matrixToVector_rc(double y[], gsl_matrix *m, int rows, int columns, int shift, int type);

//void custom_vectorToMatrix_transp(gsl_matrix *m, double y[], int rows, int columns, int shift);
//void custom_matrixToVector_transp(double y[], gsl_matrix *m, int rows, int columns, int shift);

#endif // CUSTOM_H_INCLUDED
