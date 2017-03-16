#ifndef CR3BP_DERIVATIVES_H_INCLUDED
#define CR3BP_DERIVATIVES_H_INCLUDED

#include <gsl/gsl_matrix.h>

typedef struct Bcp3d Bcp3d;
struct Bcp3d
{
    //---------------------------------------------------------------------
    // Primaries
    //---------------------------------------------------------------------
    double mu, muSEM, me, mm, ms, as;
    double xe[3];
    double xm[3];
    double delta0, beta0;
      
    //---------------------------------------------------------------------
    // Angles and coefficients
    //---------------------------------------------------------------------
    double ci, si , c2i, s2i;
    double omb, ns, omd, omtau, inc; 
    double e32, e33, e34, e35;
    double n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12;
    double q1, q2, q3, q4, q5, q6, q7, q8, q9, q10;
    double h1, h2, h3;
    double f1, f2, f3, f4;
};


int cr3bp_derivatives_42 (double t, const double y[], double f[], void *params);
int cr3bp_derivatives_6 (double t, const double y[], double f[], void *params);
int cr3bp_derivatives_48 (double t, const double y[], double f[], void *params);

int bcp_derivatives_6 (double t, const double y[], double f[], void *params);
int bcp_derivatives_42 (double t, const double y[], double f[], void *params);

int tdbcp_derivatives_6 (double t, const double y[], double f[], void *params);
int tdbcp_derivatives_42 (double t, const double y[], double f[], void *params);

#endif // CR3BP_DERIVATIVES_H_INCLUDED
