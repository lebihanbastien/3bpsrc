#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

//Custom
#include "cr3bp_derivatives.h"
#include "custom.h"

/**
* \brief compute the derivatives of  the state variables [x, y, z, xp, yp, zp]
* \param t: the current time (unused since the system is autonomous)
* \param y: the current state (dim 6)
* \param f: the derivatives to update
* \param: *params: a pointer towards the integration paramater mu
**/
int cr3bp_derivatives_6 (double t, const double y[], double f[], void *params)
{
    double mu = *(double *)params;

    //---------------------------------------------------------------------
    //Update first & second derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //---------------------------------------------------------------------
    double dU[3];

    dU[0] = (mu*(2*mu + 2*y[0] - 2))/(2*(pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0))) - ((2*mu + 2*y[0])*(mu - 1))/(2*pow(pow(mu + y[0],2.0) + pow(y[1],2.0) + pow(y[2],2.0),3.0/2.0)) - y[0];
    dU[1] = (mu*y[1])/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - y[1] - (y[1]*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),3/2.0);
    dU[2] = (mu*y[2])/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - (y[2]*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),3/2.0);

    //---------------------------------------------------------------------
    //Phase space derivatives
    //---------------------------------------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    f[3] = -dU[0] + 2*y[4];
    f[4] = -dU[1] - 2*y[3];
    f[5] = -dU[2];


    return GSL_SUCCESS;
}

/**
* \brief compute the derivatives of both the state variables [x, y, z, xp, yp, zp] and the STM in a single output f of 42 dim
* \param t: the current time (unused since the system is autonomous)
* \param y: the current state (dim 42)
* \param f: the derivatives to update
* \param: *params: a pointer towards the integration paramater mu
**/
int cr3bp_derivatives_42 (double t, const double y[], double f[], void *params)
{

    double mu = *(double *)params;
    int i,j;

    //---------------------------------------------------------------------
    //Update first & second derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //---------------------------------------------------------------------
    double dU[3];
    double dU2[3][3];

    dU[0] = (mu*(2*mu + 2*y[0] - 2))/(2*(pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0))) - ((2*mu + 2*y[0])*(mu - 1))/(2*pow(pow(mu + y[0],2.0) + pow(y[1],2.0) + pow(y[2],2.0),3.0/2.0)) - y[0];
    dU[1] = (mu*y[1])/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - y[1] - (y[1]*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),3/2.0);
    dU[2] = (mu*y[2])/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - (y[2]*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),3/2.0);

    dU2[0][0] =  mu/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - (mu - 1)/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2), 3/2.0) - (3*mu*pow(2*mu + 2*y[0] - 2,2))/(4*pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0)) + (3*pow(2*mu + 2*y[0],2)*(mu - 1))/(4*pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0)) - 1;
    dU2[0][1] = (3*y[1]*(2*mu + 2*y[0])*(mu - 1))/(2*pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0)) - (3*mu*y[1]*(2*mu + 2*y[0] - 2))/(2*pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0));
    dU2[0][2] = (3*y[2]*(2*mu + 2*y[0])*(mu - 1))/(2*pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0)) - (3*mu*y[2]*(2*mu + 2*y[0] - 2))/(2*pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0));
    dU2[1][0] = dU2[0][1];
    dU2[1][1] = mu/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - (mu - 1)/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2), 3/2.0) + (3*pow(y[1],2)*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0) - (3*mu*pow(y[1],2))/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0) - 1;
    dU2[1][2] = (3*y[1]*y[2]*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0) - (3*mu*y[1]*y[2])/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0);
    dU2[2][0] = dU2[0][2];
    dU2[2][1] = dU2[1][2];
    dU2[2][2] = mu/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - (mu - 1)/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2), 3/2.0) + (3*pow(y[2],2)*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0) - (3*mu*pow(y[2],2))/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0);

    //GSL version
    gsl_matrix *dU2_gsl = gsl_matrix_calloc(3,3);
    for(i=0; i<3 ; i++)
        for(j=0; j<3; j++) gsl_matrix_set(dU2_gsl, i, j, dU2[i][j]);
    gsl_matrix_scale (dU2_gsl, -1); //take the opposite of dU2 for later concatenation

    //---------------------------------------------------------------------
    //Phase space derivatives
    //---------------------------------------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    f[3] = -dU[0] + 2*y[4];
    f[4] = -dU[1] - 2*y[3];
    f[5] = -dU[2];


    //---------------------------------------------------------------------
    //STM derivative
    //---------------------------------------------------------------------
    //STM is updated
    gsl_matrix *STM = gsl_matrix_calloc(6,6);
    custom_vectorToMatrix(STM, (double *)y, 6, 6, 6);

    //Trivial matrices
    gsl_matrix *eye3_gsl = gsl_matrix_calloc(3,3);
    gsl_matrix *z3_gsl = gsl_matrix_calloc(3,3);
    gsl_matrix *Omega2_gsl = gsl_matrix_calloc(3,3);

    gsl_matrix_set_zero (z3_gsl);
    gsl_matrix_set_identity (eye3_gsl);
    gsl_matrix_set_zero(Omega2_gsl);

    gsl_matrix_set(Omega2_gsl, 0, 1, 2);
    gsl_matrix_set(Omega2_gsl, 1, 0, -2);

    // Building Df
    //-----------------------------------------------------------------------
    //
    //Matrix concatenation:
    //
    //      |  O          I3  |
    //Df =  |  -dU2    2Omega |
    gsl_matrix *Df_gsl = gsl_matrix_calloc(6,6);

    gsl_matrix_view Df_ul = gsl_matrix_submatrix (Df_gsl , 0 , 0 , 3 , 3 );
    gsl_matrix_view Df_ur = gsl_matrix_submatrix (Df_gsl , 0 , 3 , 3 , 3 );
    gsl_matrix_view Df_ll = gsl_matrix_submatrix (Df_gsl , 3 , 0 , 3 , 3 );
    gsl_matrix_view Df_lr = gsl_matrix_submatrix (Df_gsl , 3 , 3 , 3 , 3 );

    gsl_matrix_memcpy( &Df_ul.matrix, z3_gsl);
    gsl_matrix_memcpy( &Df_ur.matrix, eye3_gsl);
    gsl_matrix_memcpy( &Df_ll.matrix, dU2_gsl);
    gsl_matrix_memcpy( &Df_lr.matrix, Omega2_gsl);

    //Matrix product
    //-----------------------------------------------------------------------
    gsl_matrix *STMp = gsl_matrix_calloc(6,6);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Df_gsl, STM,0.0, STMp);

    //-----------------------------------------------------------------------
    //Update f (from row 6 to 41)
    //-----------------------------------------------------------------------
    custom_matrixToVector(f, STMp, 6, 6, 6);


    //Memory release
    //---------------------------------------------------------------------
    gsl_matrix_free(dU2_gsl);
    gsl_matrix_free(STM);
    gsl_matrix_free(eye3_gsl);
    gsl_matrix_free(z3_gsl);
    gsl_matrix_free(Omega2_gsl);
    gsl_matrix_free(STMp);
    gsl_matrix_free(Df_gsl);

    return GSL_SUCCESS;
}


/**
 * \brief compute the derivatives of:
 *      - The relative state variable, in y[0:5]
 *      - The absolute state of the target, in y[6:11]
 *      - The STM associated to the relative dynamics, in y[12:47]
 *
 * CAREFUL: transposed versions of custom_vectorToMatrix is used to account 
 * for the use of reshape in the original MATLAB code... 
 * Better to do that everywhere!!
 *
 *
 * \param t: the current time (unused since the system is autonomous)
 * \param y: the current state (dim 48)
 * \param f: the derivatives to update
 * \param: *params: a pointer towards the integration paramater mu
 **/
int cr3bp_derivatives_48 (double t, const double y[], double f[], void *params)
{

    double mu = *(double *)params;
    int i,j;
    
    //---------------------------------------------------------------------
    //Temporary variables
    //---------------------------------------------------------------------
    // Target state
    double x_T = y[6];
    double y_T = y[7];
    double z_T = y[8];

    // Chaser state
    double x_C = y[0] + x_T;
    double y_C = y[1] + y_T;
    double z_C = y[2] + z_T;
    
    // Define distances:
    double dET = sqrt( (x_T + mu)*(x_T + mu) + y_T*y_T + z_T*z_T );          //Target-Earth
    double dMT = sqrt( (x_T - 1 + mu)*(x_T - 1 + mu) + y_T*y_T + z_T*z_T );  //Target-Moon

    double dEC = sqrt( (x_C + mu)*(x_C + mu) + y_C*y_C + z_C*z_C );          //Chaser-Earth
    double dMC = sqrt( (x_C - 1 + mu)*(x_C - 1 + mu) + y_C*y_C + z_C*z_C );  //Chaser-Moon

    // Define constants
    double C1 = (1-mu)/(dET*dET*dET);
    double C2 = mu/(dMT*dMT*dMT);
    double C3 = (1 - mu)/(dEC*dEC*dEC);
    double C4 = mu/(dMC*dMC*dMC);
    
    //---------------------------------------------------------------------
    // Equations of motion
    //---------------------------------------------------------------------
    // Relative velocity in synodical coordinates
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];  

    // Relative acceleration in synodical coordinates
    f[3] =   2*y[4] + y[0] + C1*(x_T + mu) - C3*(x_C + mu) + C2*(x_T + mu - 1) - C4*(x_C + mu - 1);
    f[4] =  -2*y[3] + y[1] + C1*y_T - C3*y_C + C2*y_T - C4*y_C;
    f[5] =                   C1*z_T - C3*z_C + C2*z_T - C4*z_C;

    // Target velocity in synodical coordinates
    f[6] = y[9];
    f[7] = y[10];
    f[8] = y[11];

    // Target acceleration in synodical coordinates
    f[9]  =  2*y[10] + x_T  -  C1*(x_T + mu) - C2*(x_T - 1 + mu);
    f[10] = -2*y[9]  + y_T  -  C1*y_T        - C2* y_T          ;
    f[11] =                  -  C1*z_T       - C2* z_T          ;
            
            
    //---------------------------------------------------------------------
    //Update second derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //---------------------------------------------------------------------
    double dU2[3][3];
    
    dU2[0][0] = 1 - (1 - mu)/pow(dEC, 3.0) - mu/pow(dMC, 3.0) 
                + 3*(1 - mu)*(x_C + mu)*(x_C + mu)/pow(dEC, 5.0) 
                + 3*      mu*(x_C - 1 + mu)*(x_C - 1 + mu)/pow(dMC, 5.0);
            
    dU2[1][1] = 1 - (1 - mu)/pow(dEC, 3.0) - mu/pow(dMC, 3.0) 
                + 3*(1 - mu)*y_C*y_C/pow(dEC, 5.0) 
                + 3*      mu*y_C*y_C/pow(dMC, 5.0);
 
    dU2[2][2] =   - (1 - mu)/pow(dEC, 3.0) - mu/pow(dMC, 3.0) 
                 + 3*(1 - mu)*z_C*z_C/pow(dEC, 5.0) 
                 + 3*      mu*z_C*z_C/pow(dMC, 5.0);
    
    dU2[0][1] =   3*      mu*y_C*(x_C - 1 + mu)/pow(dMC, 5.0) 
                + 3*(1 - mu)*y_C*    (x_C + mu)/pow(dEC, 5.0);

    dU2[0][2] =   3*      mu*z_C*(x_C - 1 + mu)/pow(dMC, 5.0) 
                + 3*(1 - mu)*z_C*    (x_C + mu)/pow(dEC, 5.0);

    dU2[1][2] =   3*      mu*y_C*z_C/pow(dMC, 5.0) 
                + 3*(1 - mu)*y_C*z_C/pow(dEC, 5.0);

    dU2[1][0] = dU2[0][1];
    dU2[2][0] = dU2[0][2];
    dU2[2][1] = dU2[1][2];

    //GSL version
    gsl_matrix *dU2_gsl = gsl_matrix_calloc(3,3);
    for(i=0; i<3 ; i++)
        for(j=0; j<3; j++) gsl_matrix_set(dU2_gsl, i, j, dU2[i][j]);

    //---------------------------------------------------------------------
    //STM derivative
    //---------------------------------------------------------------------
    //STM is updated
    gsl_matrix *STM = gsl_matrix_calloc(6,6);
    custom_vectorToMatrix(STM, (double *)y, 6, 6, 12);

    //Trivial matrices
    gsl_matrix *eye3_gsl = gsl_matrix_calloc(3,3);
    gsl_matrix *z3_gsl = gsl_matrix_calloc(3,3);
    gsl_matrix *Omega2_gsl = gsl_matrix_calloc(3,3);

    gsl_matrix_set_zero (z3_gsl);
    gsl_matrix_set_identity (eye3_gsl);
    gsl_matrix_set_zero(Omega2_gsl);

    gsl_matrix_set(Omega2_gsl, 0, 1, 2);
    gsl_matrix_set(Omega2_gsl, 1, 0, -2);

    // Building Df
    //---------------------------------------------------------------------
    //
    //Matrix concatenation:
    //
    //      |  O          I3  |
    //Df =  |  dU2     2Omega |
    gsl_matrix *Df_gsl = gsl_matrix_calloc(6,6);

    gsl_matrix_view Df_ul = gsl_matrix_submatrix (Df_gsl , 0 , 0 , 3 , 3 );
    gsl_matrix_view Df_ur = gsl_matrix_submatrix (Df_gsl , 0 , 3 , 3 , 3 );
    gsl_matrix_view Df_ll = gsl_matrix_submatrix (Df_gsl , 3 , 0 , 3 , 3 );
    gsl_matrix_view Df_lr = gsl_matrix_submatrix (Df_gsl , 3 , 3 , 3 , 3 );

    gsl_matrix_memcpy( &Df_ul.matrix, z3_gsl);
    gsl_matrix_memcpy( &Df_ur.matrix, eye3_gsl);
    gsl_matrix_memcpy( &Df_ll.matrix, dU2_gsl);
    gsl_matrix_memcpy( &Df_lr.matrix, Omega2_gsl);

    //Matrix product
    //---------------------------------------------------------------------
    gsl_matrix *STMp = gsl_matrix_calloc(6,6);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Df_gsl, STM,0.0, STMp);

    //---------------------------------------------------------------------
    //Update f (from row 13 to 48)
    //---------------------------------------------------------------------
    custom_matrixToVector_rc(f, STMp, 6, 6, 12, COLUMNWISE);


    //Memory release
    //---------------------------------------------------------------------
    gsl_matrix_free(dU2_gsl);
    gsl_matrix_free(STM);
    gsl_matrix_free(eye3_gsl);
    gsl_matrix_free(z3_gsl);
    gsl_matrix_free(Omega2_gsl);
    gsl_matrix_free(STMp);
    gsl_matrix_free(Df_gsl);

    return GSL_SUCCESS;
}


/**
* \brief compute the derivatives of  the state variables [x, y, z, xp, yp, zp], for the bicircular problem
* \param t: the current time (unused since the system is autonomous)
* \param y: the current state (dim 6)
* \param f: the derivatives to update
* \param: *params: a pointer towards the integration paramater mu
**/
int bcp_derivatives_6 (double t, const double y[], double f[], void *params)
{
    //---------------------------------------------------------------------
    //Init
    //---------------------------------------------------------------------
    double *p = (double *)params;
    double mu     = p[0];    //crtbp mass ratio
    double omega0 = p[1];    //initial phase of the fourth body
    double ms     = p[2];    //mass of the Sun
    double as     = p[3];    //semi-major axis of the Sun
    double omegaS = p[4];    //mean motion of the Sun
   

    //Current Sun phase angle
    double theta = omega0 + omegaS*t; 
    
    //---------------------------------------------------------------------
    //Update first & second derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //---------------------------------------------------------------------
    double r1 = sqrt((y[0]+mu)*(y[0]+mu) + y[1]*y[1] + y[2]*y[2]);
    double r2 = sqrt((y[0]-1+mu)*(y[0]-1+mu) + y[1]*y[1] + y[2]*y[2]);
    double rs = sqrt((y[0]-as*cos(theta))*(y[0]-as*cos(theta)) + (y[1]-as*sin(theta))*(y[1]-as*sin(theta)) + y[2]*y[2]);
    
    double dU[3];
    dU[0] = y[0] - (1.0-mu)/pow(r1, 3.0)*(y[0]+mu) - mu/pow(r2, 3.0)*(y[0]-1+mu) - ms/pow(rs, 3.0)*(y[0]-as*cos(theta)) - ms/pow(as, 2.0)*cos(theta);
    dU[1] = y[1] - (1.0-mu)/pow(r1, 3.0)*y[1]      - mu/pow(r2, 3.0)*y[1]        - ms/pow(rs, 3.0)*(y[1]-as*sin(theta)) - ms/pow(as, 2.0)*sin(theta);
    dU[2] =      - (1.0-mu)/pow(r1, 3.0)*y[2]      - mu/pow(r2, 3.0)*y[2]        - ms/pow(rs, 3.0)*y[2];
  
    //---------------------------------------------------------------------
    //Phase space derivatives: careful, change of sign wrt CR3BP case!!
    //---------------------------------------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    f[3] = +dU[0] + 2*y[4];
    f[4] = +dU[1] - 2*y[3];
    f[5] = +dU[2];


    return GSL_SUCCESS;
}



/**
* \brief compute the derivatives of both the state variables [x, y, z, xp, yp, zp] and the STM in a single output f of 42 dim, for the bicircular problem
* \param t: the current time (unused since the system is autonomous)
* \param y: the current state (dim 42)
* \param f: the derivatives to update
* \param: *params: a pointer towards the integration paramater mu
**/
int bcp_derivatives_42 (double t, const double y[], double f[], void *params)
{
    //---------------------------------------------------------------------
    //Init
    //---------------------------------------------------------------------
    double *p = (double *)params;
    double mu     = p[0];    //crtbp mass ratio
    double omega0 = p[1];    //initial phase of the fourth body
    double ms     = p[2];    //mass of the Sun
    double as     = p[3];    //semi-major axis of the Sun
    double omegaS = p[4];    //mean motion of the Sun

    //Current Sun phase angle
    double theta = omega0 + omegaS*t; 
    
    int i,j;

    //---------------------------------------------------------------------
    //Update first derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //---------------------------------------------------------------------
    double r1 = sqrt((y[0]+mu)*(y[0]+mu) + y[1]*y[1] + y[2]*y[2]);
    double r2 = sqrt((y[0]-1+mu)*(y[0]-1+mu) + y[1]*y[1] + y[2]*y[2]);
    double rs = sqrt((y[0]-as*cos(theta))*(y[0]-as*cos(theta)) + (y[1]-as*sin(theta))*(y[1]-as*sin(theta)) + y[2]*y[2]);
      
    double dU[3];
    dU[0] = y[0] - (1.0-mu)/pow(r1, 3.0)*(y[0]+mu) - mu/pow(r2, 3.0)*(y[0]-1+mu) - ms/pow(rs, 3.0)*(y[0]-as*cos(theta)) - ms/pow(as, 2.0)*cos(theta);
    dU[1] = y[1] - (1.0-mu)/pow(r1, 3.0)*y[1]      - mu/pow(r2, 3.0)*y[1]        - ms/pow(rs, 3.0)*(y[1]-as*sin(theta)) - ms/pow(as, 2.0)*sin(theta);
    dU[2] =      - (1.0-mu)/pow(r1, 3.0)*y[2]      - mu/pow(r2, 3.0)*y[2]        - ms/pow(rs, 3.0)*y[2];

    //---------------------------------------------------------------------
    //Update second derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //---------------------------------------------------------------------
    double dU2[3][3];
    //dUx/f,y,z
    dU2[0][0] = 1 + 3.0*(1.0-mu)/pow(r1, 5.0)*(y[0]+mu)*(y[0]+mu)                 - (1.0-mu)/pow(r1, 3.0)
                  + 3.0*mu/pow(r2, 5.0)*(y[0]-1.0+mu)*(y[0]-1.0+mu)               - mu/pow(r2, 3.0)
                  + 3.0*ms/pow(rs, 5.0)*(y[0]-as*cos(theta))*(y[0]-as*cos(theta)) - ms/pow(rs, 3.0);
    
    dU2[0][1] = 0 + 3.0*(1.0-mu)/pow(r1, 5.0)*(y[0]+mu)*y[1]
                  + 3.0*mu/pow(r2, 5.0)*(y[0]-1.0+mu)*y[1]
                  + 3.0*ms/pow(rs, 5.0)*(y[0]-as*cos(theta))*(y[1]-as*sin(theta));
    
    dU2[0][2] = 0 + 3.0*(1.0-mu)/pow(r1, 5.0)*(y[0]+mu)*y[2]
                  + 3.0*mu/pow(r2, 5.0)*(y[0]-1.0+mu)*y[2]
                  + 3.0*ms/pow(rs, 5.0)*(y[0]-as*cos(theta))*y[2];
    
    //dUy/f,y,z
    dU2[1][0] = dU2[0][1];
    
    dU2[1][1] = 1 + 3.0*(1.0-mu)/pow(r1, 5.0)*y[1]*y[1]                           - (1.0-mu)/pow(r1, 3.0)
                  + 3.0*mu/pow(r2, 5.0)*y[1]*y[1]                                 - mu/pow(r2, 3.0)
                  + 3.0*ms/pow(rs, 5.0)*(y[1]-as*sin(theta))*(y[1]-as*sin(theta)) - ms/pow(rs, 3.0);
    
    dU2[1][2] = 0 + 3.0*(1.0-mu)/pow(r1, 5.0)*y[1]*y[2]
                  + 3.0*mu/pow(r2, 5.0)*y[1]*y[2]
                  + 3.0*ms/pow(rs, 5.0)*(y[1]-as*sin(theta))*y[2];
   
   //dUz/f,y,z             
   dU2[2][0]  = dU2[0][2];
   dU2[2][1]  = dU2[1][2];
   dU2[2][2]  = 0 + 3.0*(1.0-mu)/pow(r1, 5.0)*y[2]*y[2] - (1.0-mu)/pow(r1, 3.0)
                  + 3.0*mu/pow(r2, 5.0)*y[2]*y[2]       - mu/pow(r2, 3.0)
                  + 3.0*ms/pow(rs, 5.0)*y[2]*y[2]       - ms/pow(rs, 3.0);    
    
    //GSL version
    gsl_matrix *dU2_gsl = gsl_matrix_calloc(3,3);
    for(i=0; i<3 ; i++)
        for(j=0; j<3; j++) gsl_matrix_set(dU2_gsl, i, j, dU2[i][j]);
    //gsl_matrix_scale (dU2_gsl, -1); //do NOT take the opposite of dU2 here!!
    

    //---------------------------------------------------------------------
    //Phase space derivatives: careful, change of sign wrt CR3BP case!!
    //---------------------------------------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    f[3] = dU[0] + 2*y[4];
    f[4] = dU[1] - 2*y[3];
    f[5] = dU[2];


    //---------------------------------------------------------------------
    //STM derivative
    //---------------------------------------------------------------------
    //STM is updated
    gsl_matrix *STM = gsl_matrix_calloc(6,6);
    custom_vectorToMatrix(STM, (double *)y, 6, 6, 6);

    //Trivial matrices
    gsl_matrix *eye3_gsl = gsl_matrix_calloc(3,3);
    gsl_matrix *z3_gsl = gsl_matrix_calloc(3,3);
    gsl_matrix *Omega2_gsl = gsl_matrix_calloc(3,3);

    gsl_matrix_set_zero (z3_gsl);
    gsl_matrix_set_identity (eye3_gsl);
    gsl_matrix_set_zero(Omega2_gsl);

    gsl_matrix_set(Omega2_gsl, 0, 1, 2);
    gsl_matrix_set(Omega2_gsl, 1, 0, -2);

    // Building Df
    //-----------------------------------------------------------------------
    //
    //Matrix concatenation:
    //
    //      |  O          I3  |
    //Df =  |  -dU2    2Omega |
    gsl_matrix *Df_gsl = gsl_matrix_calloc(6,6);

    gsl_matrix_view Df_ul = gsl_matrix_submatrix (Df_gsl , 0 , 0 , 3 , 3 );
    gsl_matrix_view Df_ur = gsl_matrix_submatrix (Df_gsl , 0 , 3 , 3 , 3 );
    gsl_matrix_view Df_ll = gsl_matrix_submatrix (Df_gsl , 3 , 0 , 3 , 3 );
    gsl_matrix_view Df_lr = gsl_matrix_submatrix (Df_gsl , 3 , 3 , 3 , 3 );

    gsl_matrix_memcpy( &Df_ul.matrix, z3_gsl);
    gsl_matrix_memcpy( &Df_ur.matrix, eye3_gsl);
    gsl_matrix_memcpy( &Df_ll.matrix, dU2_gsl);
    gsl_matrix_memcpy( &Df_lr.matrix, Omega2_gsl);

    //Matrix product
    //-----------------------------------------------------------------------
    gsl_matrix *STMp = gsl_matrix_calloc(6,6);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Df_gsl, STM,0.0, STMp);

    //-----------------------------------------------------------------------
    //Update f (from row 6 to 41)
    //-----------------------------------------------------------------------
    custom_matrixToVector(f, STMp, 6, 6, 6);


    //Memory release
    //-----------------------------------------------------------------------
    gsl_matrix_free(dU2_gsl);
    gsl_matrix_free(STM);
    gsl_matrix_free(eye3_gsl);
    gsl_matrix_free(z3_gsl);
    gsl_matrix_free(Omega2_gsl);
    gsl_matrix_free(STMp);
    gsl_matrix_free(Df_gsl);

    return GSL_SUCCESS;


}

int bcp3dem(double *e1,     double *e2,     double *e3,
            double *e1dot,  double *e2dot,  double *e3dot,
            double *e1ddot, double *e2ddot, double *e3ddot,
            double *xs,     double delta,   double beta,
            double alpha,   Bcp3d *bcp3d)
{
    //=====================================================================
    // Variables from bcp3d structure
    //=====================================================================

    //---------------------------------------------------------------------
    // Angles and rotations
    //---------------------------------------------------------------------
    double tau = alpha + delta;
    double ci  = bcp3d->ci;
    double si  = bcp3d->si;

    //---------------------------------------------------------------------
    // Rotation rates
    //---------------------------------------------------------------------
    double omb   = bcp3d->omb;
    double omtau = bcp3d->ns + bcp3d->omd;

    //---------------------------------------------------------------------
    // Coefficients
    //---------------------------------------------------------------------
    double n1  =  bcp3d->n1;
    double n2  =  bcp3d->n2;
    double n3  =  bcp3d->n3;
    double n4  =  bcp3d->n4;
    double n6  =  bcp3d->n6;
    double n10 =  bcp3d->n10;
    double n11 =  bcp3d->n11;

    double q1  = bcp3d->q1;
    double q2  = bcp3d->q2;
    double q3  = bcp3d->q3;
    double q4  = bcp3d->q4;
    double q5  = bcp3d->q5;
    double q6  = bcp3d->q6;
    double q7  = bcp3d->q7;
    double q8  = bcp3d->q8;
    double q9  = bcp3d->q9;

    double h1  = bcp3d->h1;
    double h2  = bcp3d->h2;
    double h3  = bcp3d->h3;

    double f1  = bcp3d->f1;
    double f2  = bcp3d->f2;
    double f3  = bcp3d->f3;
    double f4  = bcp3d->f4;

    //---------------------------------------------------------------------
    // Temporary variables from current state
    //---------------------------------------------------------------------
    double cd = cos(delta);
    double sd = sin(delta);

    double ct = cos(tau);
    double st = sin(tau);

    double cb = cos(beta);
    double sb = sin(beta);

    double c2b = cos(2*beta);
    double s2b = sin(2*beta);

    // Used in previous versions of the implementation
    // Ac = (ct*cb - st*sb*ci);
    // As = (st*sb - ct*cb*ci);
    // Bc = (ct*sb + st*cb*ci);
    // Bs = (st*cb + ct*sb*ci);

    //=====================================================================
    // Normalization
    //
    // Old version, without temp variables:
    // Nddotc = -omb*omtau*omtau*si*si/(4*N*N*N)*(omb*omtau*omtau*si*si*s2b*s2b + 4*omb*c2b*N*N);
    //=====================================================================
    double N      = sqrt( (omb + omtau*ci)*(omb + omtau*ci) + omtau*omtau*cb*cb*si*si );
    double Ndotc  = -1/(2*N)*omb*omtau*omtau*si*si*s2b;
    double Nddotc = -q9/(4*N*N*N)*(q9*s2b*s2b + 4*omb*N*N*c2b);

    //=====================================================================
    // ei vectors
    //=====================================================================
    e1[0] = ct*cb - ci*st*sb;
    e1[1] = st*cb + ci*ct*sb;
    e1[2] = si*sb;
    
    e2[0] = 1.0/N*(- n1*ct*sb - n2*st*cb);
    e2[1] = 1.0/N*(+ n2*ct*cb - n1*st*sb);
    e2[2] = 1.0/N*(omb * si * cb);
    
    e3[0] = 1.0/N*(+n4*st + n11*st*sb*sb - n6*ct*cb*sb);
    e3[1] = 1.0/N*(-n4*ct - n11*ct*sb*sb - n6*st*cb*sb);
    e3[2] = 1.0/N*(n3 + n10*sb*sb + omtau*cb*cb);

    // Version without the temp variables
    // e3c = 1/N*[+si*(omb*st + omtau*st*sb*sb*ci - omtau*ct*cb*sb);
    //            -si*(omb*ct + omtau*ct*sb*sb*ci + omtau*st*cb*sb);
    //            +omb*ci + omtau*(cb*cb + sb*sb*ci^2)];

    //=====================================================================
    // xs = -as [e1(1); e2(1); e3(1)] but with delta instead of tau.
    //=====================================================================
    xs[0] = - bcp3d->as*(+cd*cb - ci*sd*sb);
    xs[1] = + bcp3d->as/N*(n1*cd*sb + n2*sd*cb);
    xs[2] = - bcp3d->as/N*(n4*sd + n11*sd*sb*sb - n6*cd*cb*sb);
    
    //=====================================================================
    // dot(ei) vectors
    //=====================================================================

    //---------------------------------------------------------------------
    // dot(e1)
    //---------------------------------------------------------------------
    e1dot[0] = - n1*ct*sb - n2*st*cb; 
    e1dot[1] = + n2*ct*cb - n1*st*sb;
    e1dot[2] = omb * si * cb;
    
    //---------------------------------------------------------------------       
    // dot(e2)
    //---------------------------------------------------------------------
    // Version without the temp variables
    // Ne2dotc = [+2*omtau*omb*As - (omtau*omtau + omb^2)*Ac;
    //            -2*omtau*omb*Bc - (omtau*omtau + omb^2)*Bs;
    //            -omb^2 * sb* si];

    // Version with the temp variables
    //     Ne2dotc = [-q1*ct*cb + q2*st*sb;
    //                -q2*ct*sb - q1*st*cb;
    //                -q3 * sb];

    e2dot[0] = (-q1*ct*cb + q2*st*sb - Ndotc*e2[0])/N;
    e2dot[1] = (-q2*ct*sb - q1*st*cb - Ndotc*e2[1])/N;
    e2dot[2] = (-q3 * sb             - Ndotc*e2[2])/N;
    

    //---------------------------------------------------------------------
    // dot(e3)
    //---------------------------------------------------------------------
    // Version without the temp variables
    // Ne3dotc = [  omb*omtau*si*ct + 1/2*omtau*omtau*s2i*ct*sb*sb 
    //            - omtau*omb*si*ct*c2b + 1/2*(omtau*omtau*si + omtau*omb*s2i)*st*s2b;          
    //              omb*omtau*si*st + 1/2*omtau*omtau*s2i*st*sb*sb 
    //            - omtau*omb*si*st*c2b - 1/2*(omtau*omtau*si + omtau*omb*s2i)*ct*s2b;            
    //            - omb*omtau*si*si*s2b];

    // Version with the temp variables
    //     Ne3dotc = [ q4*ct + q5*ct*sb*sb - q6*ct*c2b + q7*st*s2b;         
    //                 q4*st + q5*st*sb*sb - q6*st*c2b - q7*ct*s2b;            
    //               - q8*s2b];

    e3dot[0] = (q4*ct + q5*ct*sb*sb - q6*ct*c2b + q7*st*s2b - Ndotc*e3[0])/N;
    e3dot[1] = (q4*st + q5*st*sb*sb - q6*st*c2b - q7*ct*s2b - Ndotc*e3[1])/N;
    e3dot[2] = (- q8*s2b                                    - Ndotc*e3[2])/N;
    
    //=====================================================================
    // ddot(ei) vectors
    //===================================================================== 

    //---------------------------------------------------------------------
    // ddot(e1)
    //---------------------------------------------------------------------
    //     ddot(e1) = [-q1*ct*cb + q2*st*sb;
    //                 -q2*ct*sb - q1*st*cb;
    //                 -q3 * sb];
    e1ddot[0] = -q1*ct*cb + q2*st*sb;
    e1ddot[1] = -q2*ct*sb - q1*st*cb;
    e1ddot[2] = -q3 * sb;
    
    //---------------------------------------------------------------------
    // ddot(e2)
    //---------------------------------------------------------------------
    //     Ne2ddotc = [+ct*sb*h2 + st*cb*h1;
    //                 -ct*cb*h1 + st*sb*h2;
    //                 -h3*cb];

    e2ddot[0] = (+ct*sb*h2 + st*cb*h1 - Nddotc * e2[0] - 2*Ndotc*e2dot[0])/N;
    e2ddot[1] = (-ct*cb*h1 + st*sb*h2 - Nddotc * e2[1] - 2*Ndotc*e2dot[1])/N;
    e2ddot[2] = (-h3*cb               - Nddotc * e2[2] - 2*Ndotc*e2dot[2])/N;
    
    //---------------------------------------------------------------------
    // ddot(e3)
    //---------------------------------------------------------------------
    e3ddot[0] = (- f1*st + f2*ct*cb*sb + f3*st*cb*cb - Nddotc * e3[0] - 2*Ndotc*e3dot[0])/N;
    e3ddot[0] = (+ f1*ct + f2*st*cb*sb - f3*ct*cb*cb - Nddotc * e3[1] - 2*Ndotc*e3dot[1])/N;
    e3ddot[0] = (- f4*c2b                            - Nddotc * e3[2] - 2*Ndotc*e3dot[2])/N;
    
   return GSL_SUCCESS; 
}

/**
 * \brief Dot product for 3 x 1 vectors
 **/
double dot(double ex[3], double ey[3])
{
    return (ex[0]*ey[0] + ex[1]*ey[1] + ex[2]*ey[2]);
}


/**
 * \brief compute the derivatives of  the state variables [x, y, z, xp, yp, zp], for the three-dimensionnal bicircular problem
 * \param t: the current time (unused since the system is autonomous)
 * \param y: the current state (dim 6)
 * \param f: the derivatives to update
 * \param: *params: a pointer towards the integration parameters
 **/
int tdbcp_derivatives_6 (double t, const double y[], double f[], void *params)
{
    //---------------------------------------------------------------------
    //Init Bcp3d *bcp3d = (Bcp3d *) params;
    //---------------------------------------------------------------------
    Bcp3d *bcp3d = (Bcp3d *) params;

    //---------------------------------------------------------------------
    // Current angles
    //---------------------------------------------------------------------
    double delta = bcp3d->delta0 + bcp3d->omd*t;
    double beta  = bcp3d->beta0  + bcp3d->omb*t;
    double alpha = bcp3d->ns*t;

    //---------------------------------------------------------------------
    //Creating the base
    //---------------------------------------------------------------------
    double e1[3], e2[3], e3[3], xs[3];
    double e1dot[3], e2dot[3], e3dot[3];
    double e1ddot[3], e2ddot[3], e3ddot[3];
    bcp3dem(e1, e2, e3, e1dot, e2dot, e3dot, e1ddot, e2ddot, e3ddot, xs, delta, beta, alpha, bcp3d);
    
    //---------------------------------------------------------------------
    // Getting ddot(B)
    //---------------------------------------------------------------------
    double Bddot[3];
    
    Bddot[0] = - bcp3d->ns*bcp3d->ns*bcp3d->as*(1.0 - bcp3d->muSEM)*cos(alpha);
    Bddot[1] = - bcp3d->ns*bcp3d->ns*bcp3d->as*(1.0 - bcp3d->muSEM)*sin(alpha);
    Bddot[2] = + 0.0;
    
    //---------------------------------------------------------------------
    // Getting the coefficients
    //---------------------------------------------------------------------
    double b1  = - dot(Bddot, e1);
    double b2  = - dot(Bddot, e2);
    double b3  = - dot(Bddot, e3);
    double b5  = 2*dot(e1dot, e2);
    double b6  = 2*dot(e2dot, e3);
    double b7  =   dot(e1dot, e1dot);
    double b8  =   dot(e1ddot, e2);
    double b9  =   dot(e1dot, e3dot);
    double b10 =   dot(e2dot, e2dot);
    double b11 =   dot(e2ddot, e3);
    double b12 =   dot(e3dot, e3dot);
    
    //---------------------------------------------------------------------
    // Getting the potential derivative.
    //---------------------------------------------------------------------
    double re, rs, rm, dU[3];
    
    re = + (y[0] - bcp3d->xe[0])*(y[0] - bcp3d->xe[0]) 
         + (y[1] - bcp3d->xe[1])*(y[1] - bcp3d->xe[1]) 
         + (y[2] - bcp3d->xe[2])*(y[2] - bcp3d->xe[2]);
    re = sqrt(re);
    
    rm = + (y[0] - bcp3d->xm[0])*(y[0] - bcp3d->xm[0]) 
         + (y[1] - bcp3d->xm[1])*(y[1] - bcp3d->xm[1]) 
         + (y[2] - bcp3d->xm[2])*(y[2] - bcp3d->xm[2]);
    rm = sqrt(rm);
    
    rs = + (y[0] - xs[0])*(y[0] - xs[0]) 
         + (y[1] - xs[1])*(y[1] - xs[1]) 
         + (y[2] - xs[2])*(y[2] - xs[2]);
    rs = sqrt(rs);
    
    
    dU[0] = - bcp3d->me/(re*re*re)*(y[0] - bcp3d->xe[0])
            - bcp3d->mm/(rm*rm*rm)*(y[0] - bcp3d->xm[0])
            - bcp3d->ms/(rs*rs*rs)*(y[0] - xs[0]);
    
    dU[1] = - bcp3d->me/(re*re*re)*(y[1] - bcp3d->xe[1])
            - bcp3d->mm/(rm*rm*rm)*(y[1] - bcp3d->xm[1])
            - bcp3d->ms/(rs*rs*rs)*(y[1] - xs[1]);
    
    dU[2] = - bcp3d->me/(re*re*re)*(y[2] - bcp3d->xe[2])
            - bcp3d->mm/(rm*rm*rm)*(y[2] - bcp3d->xm[2])
            - bcp3d->ms/(rs*rs*rs)*(y[2] - xs[2]);


    //---------------------------------------------------------------------
    //Phase space derivatives: careful, change of sign wrt CR3BP case!!
    //---------------------------------------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    
    f[3] = b1 + b5*y[4] + b7*y[0] + b8*y[1]  + b9 *y[2]            + dU[0];
    f[4] = b2 - b5*y[3] + b6*y[5] - b8*y[0]  + b10*y[1] + b11*y[2] + dU[1];
    f[5] = b3 - b6*y[4] + b9*y[0] - b11*y[1] + b12*y[2]            + dU[2];

    return GSL_SUCCESS;
}
