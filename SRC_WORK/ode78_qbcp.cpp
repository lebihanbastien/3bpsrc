/*=========================================================================
 * ode78_qbcp.c mex file to generate
 * a Runge-Kutta 7/8 integrator of the QBCP vector field.
 *
 * @TODO: for now by default: EML2 and SEML2 are used.
 *
 * @TODO: for now, only inputs of the form PSEM, NCSEM, or VNCSEM are accepted
 * (SEM coordinates, with a state of the form (x, px).
 * Need to also accept the other forms: VSEM, NCSEM, PEM, VEM, NCEM.
 *
 * @TODO: for the dcs, only PSEM, VSEM and PEM are accepted (not VEM!!)
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
// 1. [t0, tf] the time span
// 2. double y0[6 or 42], the initial state
// 3. int dsc the current default coordinate system
//    (either PSEM, VSEM, PEM or VEM)
// 4. int nGrid nGrid the number of point to store
// 5. int inputType the type of inputs
// 6. int outputType the type of outputs
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
    if(nrhs!=6)
    {
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nrhs","6 inputs required.");
    }
    if(nlhs != 2)
    {
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nlhs", "2 outputs required.");
    }

    //---------------------------------------------------------------------
    // Retrieve the variables:
    // 1. [t0, tf] the time span
    // 2. double y0[6 or 42], the initial state
    // 3. int dsc the current default coordinate system
    //    (either PSEM, VSEM, PEM or VEM)
    // 4. int nGrid nGrid the number of point to store
    // 5. int inputType the type of inputs
    // 6. int outputType the type of outputs
    //---------------------------------------------------------------------
    double *ts     = mxGetPr(prhs[0]);
    int nts        = mxGetN(prhs[0]);
    double *y0SE   = mxGetPr(prhs[1]);
    int nvar       = mxGetN(prhs[1]);
    int mvar       = mxGetM(prhs[1]);
    int dcs        = (int) mxGetScalar(prhs[2]);
    int nGrid      = (int) mxGetScalar(prhs[3]);
    int inputType  = (int) mxGetScalar(prhs[4]);
    int outputType = (int) mxGetScalar(prhs[5]);

    //---------------------------------------------------------------------
    // Set the size of the state as the max(nvar, mvar);
    //---------------------------------------------------------------------
    nvar = nvar > mvar? nvar:mvar;

    //---------------------------------------------------------------------
    // Do some checks on the inputs
    //---------------------------------------------------------------------
    //Size of the tspan vector
    if(nts!=2)
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","The time vector (first input) must be of size 2: [t0 tf]");
    //Size of the variable vector
    if(nvar!=6 && nvar!=42)
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","The state vector (second input) must be of size 6 or 42.");

    //Type of default coordinate system (dcs) for integration
    if(dcs!=F_EM && dcs != F_SEM && dcs != F_VSEM)
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Wrong framework (only F_EM, F_SEM, and F_VSEM are allowed).");

    //Type of intputs
    if(inputType!=PSEM && inputType!=NCSEM && inputType!=VNCSEM)
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","The input type must be PSEM, NCSEM, or VNCSEM for now. The current inputType is %d", inputType);

    //---------------------------------------------------------------------
    //Type of outputs
    //---------------------------------------------------------------------

    //-------------------------------
    // 1. General check
    //-------------------------------
    if(outputType > 7)
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Unknown outputType. The current outputType is %d", outputType);

    //-------------------------------
    //2. Check the combination dcs/outputType
    //-------------------------------
    switch(dcs)
    {
        //-----------------------------------------------------------------
        // F_SEM
        //-----------------------------------------------------------------
    case F_SEM:
        switch(outputType)
        {
        case NCSEM:
        case VNCSEM:
        case PSEM:
        case VSEM:
        case NCEM:
        case PEM:
        case VEM:
            break;
        default:
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Unknown outputType (%d) for the current dcs (%d)", outputType, dcs);
        }
        break;

        //-------------------------------------------------------------
        // F_VSEM
        //-------------------------------------------------------------
    case F_VSEM:
        switch(outputType)
        {
        case VNCSEM:
        case NCSEM:
        case PSEM:
        case VSEM:
        case NCEM:
        case PEM:
        case VEM:
            break;
        default:
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Unknown outputType (%d) for the current dcs (%d)", outputType, dcs);
        }
        break;

        //-----------------------------------------------------------------
        // F_EM
        //-----------------------------------------------------------------
    case F_EM:
        switch(outputType)
        {
        case NCSEM:
        case PSEM:
        case VSEM:
        case NCEM:
        case PEM:
        case VEM:
            break;
        default:
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Unknown outputType (%d) for the current dcs (%d)", outputType, dcs);
        }
        break;
    }

    //-------------------------------
    //3. Check the combination dcs/inputType
    //-------------------------------
    switch(dcs)
    {
        //-----------------------------------------------------------------
        // From SEM (PSEM) to NCSEM
        //-----------------------------------------------------------------
    case F_SEM:
        switch(inputType)
        {
        case PSEM:
        case NCSEM:
        case VNCSEM:
            break;
        default:
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Unknown inputType (%d) for the current dcs (%d)", inputType, dcs);
        }
        break;

        //-------------------------------------------------------------
        // From SEM (PSEM) to VNCSEM
        //-------------------------------------------------------------
    case F_VSEM:
        switch(inputType)
        {
        case PSEM:
        case NCSEM:
        case VNCSEM:
            break;
        default:
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Unknown inputType (%d) for the current dcs (%d)", inputType, dcs);
        }
        break;

        //-------------------------------------------------------------
        // From SEM (PSEM) to NCEM
        //-------------------------------------------------------------
    case F_EM:
        switch(inputType)
        {
        case PSEM:
            break;
        default:
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Unknown inputType (%d) for the current dcs (%d)", inputType, dcs);
        }
        break;

        //-------------------------------------------------------------
        // From SEM (PSEM) to VNCEM
        //-------------------------------------------------------------
    case F_VEM:
        switch(inputType)
        {
        case PSEM:
            break;
        default:
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Unknown inputType (%d) for the current dcs (%d)", inputType, dcs);
        }
        break;
    }


    //If nvar == 42, the dcs and the outputType must match, otherwise the variational equations make no sense with the outputs
    if(nvar == 42)
    {
        if(dcs == F_EM && outputType != NCEM)
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == F_VEM && outputType != VNCEM)
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == F_SEM && outputType != NCSEM)
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == F_VSEM && outputType != VNCSEM)
            mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
    }

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
    // Build the times
    //---------------------------------------------------------------------
    double t0NC = ts[0];
    double tfNC = ts[1];

    //---------------------------------------------------------------------
    // Selection of the vector field
    //---------------------------------------------------------------------
    int (*vf)(double, const double*, double*, void*);
    switch(nvar)
    {
        //---------------------------------------------------------------------
        // 6 variables: only the state is integrated
        //---------------------------------------------------------------------
    case 6:
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
        break;
        //---------------------------------------------------------------------
        // 42 variables: the state + var. eq. are integrated
        //---------------------------------------------------------------------
    case 42:
        switch(dcs)
        {
        case F_SEM:
        case F_EM:
            vf = qbfbp_vfn_varnonlin; //vector field with a state (x, px) (default)
            break;

        case F_VEM:
        case F_VSEM:
            vf = qbfbp_vfn_varnonlin_xv; //vector field with a state (x, vx)
            break;
        }
        break;
    }


    //---------------------------------------------------------------------
    // ODE system
    //---------------------------------------------------------------------
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver,
                       T,               //stepper: int
                       T_root,          //stepper: root
                       1e-14,           //precision: int_abs
                       1e-14,           //precision: int_rel
                       1e-13,           //precision: root
                       1e-12,           //precision: diffcorr
                       nvar,            //dimension
                       1e-6,            //initial int step
                       vf,              //vector field
                       NULL,            //jacobian
                       &SEML);          //current four-body system


    //---------------------------------------------------------------------
    // to NCFWRK coordinates.
    // Note: this part should be changed when other input types are accepted
    //---------------------------------------------------------------------
    double y0[nvar];
    double t0, tf;
    double ytemp[6], y0NCSE[6];
    switch(dcs)
    {
        //-----------------------------------------------------------------
        // If F_SEM/F_VSEM, the system is focused on the SEM system via
        // SEML
        //-----------------------------------------------------------------

        //-----------------------------------------------------------------
        // From SEM (PSEM) to NCSEM
        //-----------------------------------------------------------------
    case F_SEM:
        switch(inputType)
        {
        case PSEM:
            //From SEM to NC SEM
            SYStoNC(t0NC, y0SE, y0, &SEML);
            //Time is already in SEM units
            t0 = t0NC;
            tf = tfNC;
            break;
b
        case NCSEM:
            //From NCSEM to NCSEM: nothing to do
            for(int i = 0; i < 6; i++) y0[i] = y0SE[i];
            //Time is already in SEM units
            t0 = t0NC;
            tf = tfNC;
            break;

        case VNCSEM:
            //From VNCSEM to NCSEM
            SYSvtoSYSm(t0NC, y0SE, y0, &SEML);
            //Time is already in SEM units
            t0 = t0NC;
            tf = tfNC;
            break;
        }
        break;

        //-------------------------------------------------------------
        // From SEM (PSEM) to VNCSEM
        //-------------------------------------------------------------
    case F_VSEM:
        switch(inputType)
        {
        case PSEM:
            //From SEM to NCSEM
            SYStoNC(t0NC, y0SE, ytemp, &SEML);
            //From NCSEM to VNCSEM
            SYSmtoSYSv(t0NC, ytemp, y0, &SEML);
            //Time is already in SEM units
            t0 = t0NC;
            tf = tfNC;
            break;
        case NCSEM:
            //From NCSEM to VNCSEM
            SYSmtoSYSv(t0NC, y0SE, y0, &SEML);
            //Time is already in SEM units
            t0 = t0NC;
            tf = tfNC;
            break;

        case VNCSEM:
            //From VNCSEM to VNCSEM: nothing to do
            for(int i = 0; i < 6; i++) y0[i] = y0SE[i];
            //Time is already in SEM units
            t0 = t0NC;
            tf = tfNC;
            break;
        }
        break;

        //-------------------------------------------------------------
        // If F_EM/F_VEM, the system is focused on the EM system via
        // SEML
        //-------------------------------------------------------------
        //-------------------------------------------------------------
        // From SEM (PSEM) to NCEM
        //-------------------------------------------------------------
    case F_EM:
        //From SEM to NCSEM
        SEMtoNC(t0NC, y0SE, y0NCSE, &SEML);
        //From  NCSEM to NCEM
        NCSEMmtoNCEMm(t0NC, y0NCSE, y0, &SEML);
        //Time is now in EM units
        t0 = t0NC/SEML.us_em.ns;
        tf = tfNC/SEML.us_em.ns;
        break;

        //-------------------------------------------------------------
        // From SEM (PSEM) to VNCEM
        //-------------------------------------------------------------
    case F_VEM:
        //Time is now in EM units
        t0 = t0NC/SEML.us_em.ns;
        tf = tfNC/SEML.us_em.ns;
        //From  SEM to NCSEM
        SEMtoNC(t0NC, y0SE, y0NCSE, &SEML);
        //From  NCSEM to NCEM
        NCSEMmtoNCEMm(t0NC, y0NCSE, ytemp, &SEML);
        //From NCEM to VNCEM, using EM time
        SYSmtoSYSv(t0, ytemp, y0, &SEML);
        break;
    }


    //---------------------------------------------------------------------
    // Add the identity matrix if necessary
    //---------------------------------------------------------------------
    if(nvar == 42)
    {
        //Identity matrix eye(6)
        gsl_matrix *Id = gsl_matrix_alloc(6,6);
        gsl_matrix_set_identity (Id);
        //Storing eye(6) into the initial vector
        gslc_matrixToVector(y0, Id, 6, 6, 6);
        gsl_matrix_free(Id);
    }

    //---------------------------------------------------------------------
    // Integration, for 2 outputs
    //---------------------------------------------------------------------
    double **yvIN, **yv, *tv, *tvIN, t, y[nvar];
    int nV;
    yvIN  = dmatrix(0, nvar-1, 0, nGrid);
    yv    = dmatrix(0, nvar-1, 0, nGrid);
    tvIN  = dvector(0, nGrid);
    tv    = dvector(0, nGrid);

    //Integration in yvIN/tvIN
    ode78_qbcp_grid(&driver, t0, tf, y0, yvIN, tvIN, nGrid);

    //----------------------------------------------------------------------
    // To the right outputs
    //----------------------------------------------------------------------
    double **yMtemp;
    switch(dcs)
    {
        //-----------------------------------------------------------------
        // If F_SEM/F_VSEM, the system is focused on the SEM system via SEML
        //-----------------------------------------------------------------

        //-----------------------------------------------------------------
        // F_SEM
        //-----------------------------------------------------------------
    case F_SEM:
        switch(outputType)
        {
        case NCSEM:
            //NCSEM -> NCSEM. Nothing to do here, only copy
            for(int ki = 0; ki <= nGrid; ki++)
            {
                for(int kk = 0; kk < nvar; kk++) yv[kk][ki] = yvIN[kk][ki];
                tv[ki] = tvIN[ki];
            }
            break;

        case VNCSEM:
            //NCSEM -> VNCSEM
            NCtoVNC_vec(yvIN, tv, yv, nGrid, &SEML);
            //Copy for the time vector
            for(int ki = 0; ki <= nGrid; ki++) tv[ki] = tvIN[ki];
            break;

        case PSEM:
            //NCSEM -> PSEM
            NCtoSYS_vec(yvIN, tv, yv, nGrid, &SEML);
            //Copy for the time vector
            for(int ki = 0; ki <= nGrid; ki++) tv[ki] = tvIN[ki];
            break;
        case VSEM:
            //NCSEM -> VSEM
            NCtoSYSv_vec(yvIN, tv, yv, nGrid, &SEML);
            //Copy for the time vector
            for(int ki = 0; ki <= nGrid; ki++) tv[ki] = tvIN[ki];
            break;
        case NCEM:
            //NCSEM -> NCEM
            NCSEMmtoNCEMm_vec(yvIN, tvIN, yv, tv, nGrid, &SEML);
            break;
        case PEM:
            //NCSEM -> PEM
            NCSEMmtoEMm_vec(yvIN, tvIN, yv, tv, nGrid, &SEML);
            break;
        case VEM:
            //NCSEM -> VEM
            NCSEMmtoEMv_vec(yvIN, tvIN, yv, tv, nGrid, &SEML);
            break;
        }
        break;

        //-----------------------------------------------------------------
        // F_VSEM
        //-----------------------------------------------------------------
    case F_VSEM:
        switch(outputType)
        {
        case VNCSEM:
            //VNCSEM -> VNCSEM. Nothing to do here, only copy
            for(int ki = 0; ki <= nGrid; ki++)
            {
                for(int kk = 0; kk < nvar; kk++) yv[kk][ki] = yvIN[kk][ki];
                tv[ki] = tvIN[ki];
            }
            break;

        case NCSEM:
            //VNCSEM -> NCSEM
            VNCtoNC_vec(yvIN, tv, yv, nGrid, &SEML);
            //Copy for the time vector
            for(int ki = 0; ki <= nGrid; ki++) tv[ki] = tvIN[ki];
            break;

        case PSEM:
            //Temp variable
            yMtemp = dmatrix(0, nvar-1, 0, nGrid);
            //VNCSEM -> NCSEM
            VNCtoNC_vec(yvIN, tv, yMtemp, nGrid, &SEML);
            //NCSEM -> PSEM
            NCtoSYS_vec(yMtemp, tv, yv, nGrid, &SEML);
            //Copy for the time vector
            for(int ki = 0; ki <= nGrid; ki++) tv[ki] = tvIN[ki];
            //Free temp variable
            free_dmatrix(yMtemp, 0, 5, 0, nGrid);
            break;


        case VSEM:
            //Temp variable
            yMtemp = dmatrix(0, nvar-1, 0, nGrid);
            //VNCSEM -> NCSEM
            VNCtoNC_vec(yvIN, tv, yMtemp, nGrid, &SEML);
            //NCSEM -> VSEM
            NCtoSYSv_vec(yMtemp, tv, yv, nGrid, &SEML);
            //Copy for the time vector
            for(int ki = 0; ki <= nGrid; ki++) tv[ki] = tvIN[ki];
            //Free temp variable
            free_dmatrix(yMtemp, 0, 5, 0, nGrid);
            break;

        case NCEM:
            //Temp variable
            yMtemp = dmatrix(0, nvar-1, 0, nGrid);
            //VNCSEM -> NCSEM
            VNCtoNC_vec(yvIN, tv, yMtemp, nGrid, &SEML);
            //NCSEM -> NCEM
            NCSEMmtoNCEMm_vec(yMtemp, tvIN, yv, tv, nGrid, &SEML);
            //Free temp variable
            free_dmatrix(yMtemp, 0, 5, 0, nGrid);
            break;

        case PEM:
            //Temp variable
            yMtemp = dmatrix(0, nvar-1, 0, nGrid);
            //VNCSEM -> NCSEM
            VNCtoNC_vec(yvIN, tv, yMtemp, nGrid, &SEML);
            //NCSEM -> PEM
            NCSEMmtoEMm_vec(yMtemp, tvIN, yv, tv, nGrid, &SEML);
            //Free temp variable
            free_dmatrix(yMtemp, 0, 5, 0, nGrid);
            break;

        case VEM:
            //Temp variable
            yMtemp = dmatrix(0, nvar-1, 0, nGrid);
            //VNCSEM -> NCSEM
            VNCtoNC_vec(yvIN, tv, yMtemp, nGrid, &SEML);
            //NCSEM -> VEM
            NCSEMmtoEMv_vec(yMtemp, tvIN, yv, tv, nGrid, &SEML);
            //Free temp variable
            free_dmatrix(yMtemp, 0, 5, 0, nGrid);
            break;
        }
        break;

    case F_EM:
        switch(outputType)
        {
        case NCSEM:
            //NCEM -> NCSEM.
            NCEMmtoNCSEMm_vec(yvIN, tvIN, yv, tv, nGrid, &SEML);
            break;
        case PSEM:
            //NCEM -> PSEM
            NCEMmtoSEMm_vec(yvIN, tvIN, yv, tv, nGrid, &SEML);
            break;
        case VSEM:
            //NCEM -> VSEM
            NCEMmtoSEMv_vec(yvIN, tvIN, yv, tv, nGrid, &SEML);
        case NCEM:
            //NCEM -> NCEM. Nothing to do here, only copy
            for(int ki = 0; ki <= nGrid; ki++)
            {
                for(int kk = 0; kk < nvar; kk++) yv[kk][ki] = yvIN[kk][ki];
                tv[ki] = tvIN[ki];
            }
            break;
        case PEM:
            //NCEM -> PEM
            NCtoSYS_vec(yvIN, tv, yv, nGrid, &SEML);
            //Copy for the time vector
            for(int ki = 0; ki <= nGrid; ki++) tv[ki] = tvIN[ki];
            break;
        case VEM:
            //NCSEM -> VSEM
            NCtoSYSv_vec(yvIN, tv, yv, nGrid, &SEML);
            //Copy for the time vector
            for(int ki = 0; ki <= nGrid; ki++) tv[ki] = tvIN[ki];
            break;
        }
        break;
    }

    //---------------------------------------------------------------------
    // Output: the final state
    //---------------------------------------------------------------------
    //Create the output matrices
    plhs[0] = mxCreateDoubleMatrix(nGrid+1, 1, mxREAL);     //tv
    plhs[1] = mxCreateDoubleMatrix(nGrid+1, nvar, mxREAL);  //yv

    //Get a pointer to the real data in the output
    double *tvout = mxGetPr(plhs[0]);
    double *yvout = mxGetPr(plhs[1]);

    //Store the state on the grid [0,..., nGrid]
    int indix = 0;
    for(int i = 0; i < nvar; i++)
    {
        for(int k = 0; k <= nGrid; k++)
        {
            yvout[indix++] = yv[i][k];
        }
    }
    for(int k = 0; k <= nGrid; k++) tvout[k] = tv[k];


    //---------------------------------------------------------------------
    // Free memory
    //---------------------------------------------------------------------
    free_dmatrix(yvIN, 0, nvar-1, 0, nGrid);
    free_dmatrix(yv, 0, nvar-1, 0, nGrid);
    free_dvector(tvIN, 0, nGrid);
    free_dvector(tv, 0, nGrid);
}
