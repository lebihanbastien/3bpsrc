#include <math.h>
#include <gsl/gsl_matrix.h>
#include "custom.h"
/**
 * \brief custom routine to transform vector to matrix with a 
 * given shift in the initial vector. Storing is COLUMNWISE or ROWWISE,
 * COLUMNWISE being the equivalent of reshape of MATLAB.
 **/
void custom_vectorToMatrix_rc(gsl_matrix *m, double y[], int rows, int columns, int shift, int type)
{
    int i,j,k;
    
    switch(type)
    {
        case ROWWISE:
        {     
            for(i=1; i<=rows ; i++)
                for(j=1; j<=columns; j++)
                {
                    k = rows*(i-1)+j;
                    gsl_matrix_set(m, i-1, j-1, y[shift+k-1]);           
                }
            break;
         }
        
        case COLUMNWISE:
        default:
        {
            for(i=1; i<=rows ; i++)
                for(j=1; j<=columns; j++)
                {
                    k = columns*(j-1)+i;
                    gsl_matrix_set(m, i-1, j-1, y[shift+k-1]);           
                }
            
         break;   
        }  
    }
    
}

/**
 * \brief custom routine to transform matrix to vector with a 
 * given shift in the final vector. Storing is COLUMNWISE or ROWWISE.
 * COLUMNWISE being the equivalent of reshape of MATLAB.
 **/
void custom_matrixToVector_rc(double y[], gsl_matrix *m, int rows, int columns, int shift, int type)
{
    int i,j,k;
    
    switch(type)
    {
        case ROWWISE:
        {     
            for(i=1; i<=rows ; i++)
                for(j=1; j<=columns; j++)
                {
                    k = rows*(i-1)+j;
                    y[shift+k-1] = gsl_matrix_get(m, i-1, j-1);
                }
            break;
         }
        
        case COLUMNWISE:
        default:
        {
            for(i=1; i<=rows ; i++)
                for(j=1; j<=columns; j++)
                {
                    k = columns*(j-1)+i;
                    gsl_matrix_set(m, i-1, j-1, y[shift+k-1]);           
                }
            
         break;   
        }  
    }
}


/**
 * \brief custom routine to transform vector to matrix with a 
 * given shift in the initial vector. Storing is COLUMNWISE or ROWWISE,
 * COLUMNWISE being the equivalent of reshape of MATLAB.
 **/
void custom_vectorToMatrix(gsl_matrix *m, double y[], int rows, int columns, int shift)
{
    custom_vectorToMatrix_rc(m, y, rows, columns, shift, ROWWISE);
}

/**
 * \brief custom routine to transform matrix to vector with a 
 * given shift in the final vector. Storing is COLUMNWISE or ROWWISE.
 * COLUMNWISE being the equivalent of reshape of MATLAB.
 **/
void custom_matrixToVector(double y[], gsl_matrix *m, int rows, int columns, int shift)
{
    custom_matrixToVector_rc(y, m, rows, columns, shift, ROWWISE);
}
