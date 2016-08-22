/*=========================================================================
 * qbcp_vfn_novar.c mex file to get the QBCP vector field at given time and 
 * state
 *
 * author:  BLB
 * year:    2016
 * version: 1.0
 *=======================================================================*/

//-------------------------------------------------------------------------
// Headers
//-------------------------------------------------------------------------
//Mex
#include "mex.h"

//Custom
#include "Cpp2/FTA.h"
#include "Cpp2/Ofsc.h"
#include "Cpp2/Oftsc.h"
#include "Cpp2/env.h"
#include "Cpp2/pmcoc.h"
#include "Cpp2/ode.h"
#include "Cpp2/vf.h"
#include "Cpp2/single_orbit.h"
#include "Cpp2/parameters.h"
#include "Cpp2/eminsem.h"

//-------------------------------------------------------------------------
//Frameworks
//-------------------------------------------------------------------------
#define F_EM   0
#define F_SEM  1
#define F_VEM  2
#define F_VSEM 3

//-------------------------------------------------------------------------
// Available coordinates system
//-------------------------------------------------------------------------
#define NCSEM  0
#define NCEM   1
#define VNCSEM 2
#define VNCEM  3
#define PSEM   4
#define PEM    5
#define VSEM   6
#define VEM    7

//-------------------------------------------------------------------------
// The gateway function.
// The input must be, in that order:
// 1. t the current time
// 2. double y0(1:6) the current state
// 3. int dsc the current default coordinate system 
//    (either PSEM, VSEM, PEM or VEM)
//-------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //---------------------------------------------------------------------
    // Set the global variables
    //---------------------------------------------------------------------
    OFTS_ORDER = 12;
    OFS_ORDER  = 30;

    //---------------------------------------------------------------------
    // Check for proper number of arguments
    //---------------------------------------------------------------------
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nrhs","3 inputs required.");
    }
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nlhs","1 output required.");
    }

    //---------------------------------------------------------------------
    // Retrieve the variables:
    // 1. t the current time
    // 2. double y0(1:6) the current state
    // 3. int dsc the current default coordinate system 
    //    (either PSEM, VSEM, PEM or VEM)
    //---------------------------------------------------------------------
    double t       = mxGetScalar(prhs[0]);
    double *y0     = mxGetPr(prhs[1]);
    int nvar       = mxGetN(prhs[1]);
    int mvar       = mxGetM(prhs[1]);
    int dcs        = (int) mxGetScalar(prhs[2]);

    //---------------------------------------------------------------------
    // Set the size of the state as the max(nvar, mvar);
    //---------------------------------------------------------------------
    nvar = nvar > mvar? nvar:mvar;

    //---------------------------------------------------------------------
    // Do some checks on the inputs
    //---------------------------------------------------------------------
    //Size of the variable vector
    if(nvar!=6) 
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","The state vector (second input) must be of size 6.");
    
    //Type of default coordinate system (dcs) for integration
    if(dcs!=F_EM && dcs != F_SEM && dcs != F_VSEM && dcs != F_VEM )
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Wrong framework (only F_EM, F_SEM,  F_VSEM, and F_VEM are allowed).");
    
    //---------------------------------------------------------------------
    // Retrieve the configuration parameters
    // For now, by default: EML2 and SEML2 are used
    //---------------------------------------------------------------------
    int model     = M_QBCP;
    int isNorm    = 1;
    int li_EM     = 2;
    int li_SEM    = 2;
    int pms       = PMS_GRAPH;
    int mType_EM  = MAN_CENTER;
    int mType_SEM = MAN_CENTER;
    
    //---------------------------------------------------------------------
    //Define the framework from the default coordinate system
    //---------------------------------------------------------------------
    int fwrk;
    switch(dcs)
    {
        case F_SEM:
        case F_VSEM:
            fwrk = F_SEM;
        break;
        case F_EM:
        case F_VEM:
            fwrk = F_EM;
        break;
    }

    //---------------------------------------------------------------------
    // Computation
    //---------------------------------------------------------------------
    //------------------------------------------
    // Initialization of the environnement
    // Mandatory to perform any computation except qbtbp(int)
    //------------------------------------------
    init_env(li_EM, li_SEM, isNorm, model, fwrk, pms, pms, mType_EM, mType_SEM);


    //---------------------------------------------------------------------
    // Selection of the vector field
    //---------------------------------------------------------------------
    int (*vf)(double, const double*, double*, void*);
    switch(dcs)
    {
        case F_SEM:
        case F_EM:
            vf = qbfbp_vfn_novar;   //vector field with a state (x, px) (default)
            break;
            
        case F_VEM:
        case F_VSEM:
            vf = qbfbp_vfn_novar_xv; //vector field with a state (x, vx)
            break;
    }


   
    //---------------------------------------------------------------------
    // Updating the vector field
    //---------------------------------------------------------------------
    double f[6];
    vf(t, y0, f, &SEML);
    
   
    //---------------------------------------------------------------------
    // Output
    //---------------------------------------------------------------------
    //Create the output matrices
    plhs[0] = mxCreateDoubleMatrix(1, 6, mxREAL);     //tv

    //Get a pointer to the real data in the output
    double *yvout = mxGetPr(plhs[0]);

    //Store the state on the grid [0,..., nGrid]
    for(int i = 0; i < nvar; i++)
    {
       yvout[i] = f[i];
    }

}
