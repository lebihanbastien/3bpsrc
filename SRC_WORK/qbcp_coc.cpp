/*=========================================================================
 * qbcp_coc.c mex file to perform any change of coordinates in the QBCP
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
// Available output/input types
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
// 3. int inputType the type of the inputs
// 4. int outputType the type of the outputs
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
    if(nrhs!=4)
    {
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nrhs","4 inputs required.");
    }
    if(nlhs != 1 && nlhs != 2)
    {
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nlhs","1 (y) or 2 (y, t) outputs required.");
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
    int inputType  = (int) mxGetScalar(prhs[2]);
    int outputType = (int) mxGetScalar(prhs[3]);


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

    //Type of inputs
    if(inputType > 7)
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Unknown inputType: %d.", inputType);


    //Type of outputs
    if(outputType > 7)
        mexErrMsgIdAndTxt("custom:ode78_qbcp:nts","Unknown outputType: %d.", outputType);

    //---------------------------------------------------------------------
    // Define the default framework wrt the inputType
    //---------------------------------------------------------------------
    int fwrk;
    switch(inputType)
    {
    case VNCEM:
    case NCEM:
    case PEM:
    case VEM:
        fwrk = F_EM;
        break;
    case VNCSEM:
    case NCSEM:
    case PSEM:
    case VSEM:
        fwrk = F_SEM;
        break;
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
    // Computation
    //---------------------------------------------------------------------
    //------------------------------------------
    // Initialization of the environnement
    // Mandatory to perform any computation except qbtbp(int)
    //------------------------------------------
    init_env(li_EM, li_SEM, isNorm, model, fwrk, pms, pms, mType_EM, mType_SEM);

    //---------------------------------------------------------------------
    // Updating output
    //---------------------------------------------------------------------
    double yout[6], ytemp[6], ytemp2[6], ytemp3[6];

    switch(inputType)
    {
        //-----------------------------------------------------------------
        // For inputType = VNCEM, NCEM, PEM, VEM, the focus of SEML is on the EM
        // system
        //-----------------------------------------------------------------
        //-----------------------------------------------------------------
        // VNCEM
        //-----------------------------------------------------------------
    case VNCEM:
        switch(outputType)
        {
        case VNCEM:
            //-----------------------------------------------------
            //VNCEM -> VNCEM nothing to do
            //-----------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            break;

        case NCEM:
            //-----------------------------------------------------
            //VNCEM -> NCEM nothing to do
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, yout, &SEML);
            break;

        case PEM:
            //-----------------------------------------------------
            //VNCEM -> PEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCtoSYS(t, ytemp, yout, &SEML);
            break;

        case VEM:
            //-----------------------------------------------------
            //VNCEM -> VEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCtoSYS(t, ytemp, ytemp2, &SEML);
            SYSmtoSYSv(t, ytemp2, yout, &SEML);
            break;

        case NCSEM:
            //-----------------------------------------------------
            //VNCEM -> NCSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCEMmtoNCSEMm(t, ytemp, yout, &SEML);
            break;

        case PSEM:
            //-----------------------------------------------------
            //VNCEM -> PSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCEMmtoSEMm(t, ytemp, yout, &SEML);
            break;

        case VSEM:
            //-----------------------------------------------------
            //VNCEM -> VSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCEMmtoSEMm(t, ytemp, ytemp2, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            SYSmtoSYSv(t*SEML.us_em.ns, ytemp2, yout, &SEML_SEM);
            break;

        case VNCSEM:
            //-----------------------------------------------------
            //VNCEM -> VNCSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCEMmtoNCSEMm(t, ytemp, ytemp2, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            SYSmtoSYSv(t*SEML.us_em.ns, ytemp2, yout, &SEML_SEM);
            break;
        }
        break;


        //-----------------------------------------------------------------
        // NCEM
        //-----------------------------------------------------------------
    case NCEM:
        switch(outputType)
        {
        case NCEM:
            //-----------------------------------------------------
            //NCEM -> NCEM nothing to do
            //-----------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            break;

        case PEM:
            //-----------------------------------------------------
            //NCEM -> PEM
            //-----------------------------------------------------
            NCtoSYS(t, y0, yout, &SEML);
            break;

        case VEM:
            //-----------------------------------------------------
            //NCEM -> VEM
            //-----------------------------------------------------
            NCtoSYS(t, y0, ytemp, &SEML);
            SYSmtoSYSv(t, ytemp, yout, &SEML);
            break;

        case NCSEM:
            //-----------------------------------------------------
            //NCEM -> NCSEM
            //-----------------------------------------------------
            NCEMmtoNCSEMm(t, y0, yout, &SEML);
            break;

        case PSEM:
            //-----------------------------------------------------
            //NCEM -> PSEM
            //-----------------------------------------------------
            NCEMmtoSEMm(t, y0, yout, &SEML);
            break;

        case VSEM:
            //-----------------------------------------------------
            //NCEM -> VSEM
            //-----------------------------------------------------
            NCEMmtoSEMm(t, y0, ytemp, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            SYSmtoSYSv(t*SEML.us_em.ns, ytemp, yout, &SEML_SEM);
            break;

        case VNCSEM:
            //-----------------------------------------------------
            //NCEM -> VNCSEM
            //-----------------------------------------------------
            NCEMmtoNCSEMm(t, y0, ytemp, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            SYSmtoSYSv(t*SEML.us_em.ns, ytemp, yout, &SEML_SEM);
            break;
        }
        break;

        //-----------------------------------------------------------------
        // PEM
        //-----------------------------------------------------------------
    case PEM:
        switch(outputType)
        {
        case NCEM:
            //-----------------------------------------------------
            //PEM -> NCEM
            //-----------------------------------------------------
            SYStoNC(t, y0, yout, &SEML);
            break;

        case PEM:
            //-----------------------------------------------------
            //PEM -> PEM
            //-----------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            break;

        case VEM:
            //-----------------------------------------------------
            // PEM -> VEM
            //-----------------------------------------------------
            SYSmtoSYSv(t, y0, yout, &SEML);
            break;

        case NCSEM:
            //-----------------------------------------------------
            //PEM -> NCSEM
            //-----------------------------------------------------
            EMmtoSEMm(t, y0, ytemp, &SEML);
            //Careful here: the time must be in SEM units.
            SEMtoNC(t*SEML.us_em.ns, ytemp, yout, &SEML);
            break;

        case PSEM:
            //-----------------------------------------------------
            //PEM -> PSEM
            //-----------------------------------------------------
            EMmtoSEMm(t, y0, yout, &SEML);
            break;

        case VSEM:
            //-----------------------------------------------------
            //PEM -> VSEM
            //-----------------------------------------------------
            EMmtoSEMm(t, y0, ytemp, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            SYSmtoSYSv(t*SEML.us_em.ns, ytemp, yout, &SEML_SEM);
            break;

        case VNCSEM:
            //-----------------------------------------------------
            //PEM -> VNCSEM
            //-----------------------------------------------------
            SYStoNC(t, y0, ytemp, &SEML);
            NCEMmtoNCSEMm(t, ytemp, ytemp2, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            SYSmtoSYSv(t*SEML.us_em.ns, ytemp2, yout, &SEML_SEM);
            break;
        }
        break;

        //-----------------------------------------------------------------
        // VEM
        //-----------------------------------------------------------------
    case VEM:
        switch(outputType)
        {
        case NCEM:
            //-----------------------------------------------------
            // VEM -> NCEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            SYStoNC(t, ytemp, yout, &SEML);
            break;

        case PEM:
            //-----------------------------------------------------
            // VEM -> PEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, yout, &SEML);
            break;

        case VEM:
            //-----------------------------------------------------
            // VEM -> VEM
            //-----------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            break;

        case NCSEM:
            //-----------------------------------------------------
            // VEM -> NCSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            EMmtoSEMm(t, ytemp, ytemp2, &SEML);
            //Careful here: the time must be in SEM units.
            SEMtoNC(t*SEML.us_em.ns, ytemp2, yout, &SEML);
            break;

        case PSEM:
            //-----------------------------------------------------
            // VEM -> PSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            EMmtoSEMm(t, ytemp, yout, &SEML);
            break;

        case VSEM:
            //-----------------------------------------------------
            // VEM -> VSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            EMmtoSEMm(t, ytemp, ytemp2, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            SYSmtoSYSv(t*SEML.us_em.ns, ytemp2, yout, &SEML_SEM);
            break;

        case VNCSEM:
            //-----------------------------------------------------
            // VEM -> VNCSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            SYStoNC(t, ytemp, ytemp2, &SEML);
            NCEMmtoNCSEMm(t, ytemp2, ytemp3, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            SYSmtoSYSv(t*SEML.us_em.ns, ytemp3, yout, &SEML_SEM);
            break;
        }
        break;

        //-----------------------------------------------------------------
        // For inputType = VNCSEM, NCSEM, PSEM, VSEM, the focus of SEML is on the SEM
        // system
        //-----------------------------------------------------------------
        //-----------------------------------------------------------------
        // VNCSEM
        //-----------------------------------------------------------------
    case VNCSEM:
        switch(outputType)
        {
        case VNCEM:
            //-----------------------------------------------------
            //VNCSEM -> VNCEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCSEMmtoNCEMm(t, ytemp, ytemp2, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            SYSmtoSYSv(t/SEML.us_em.ns, ytemp2, yout, &SEML_EM);
            break;

        case NCEM:
            //-----------------------------------------------------
            //VNCSEM -> NCEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCSEMmtoNCEMm(t, ytemp, yout, &SEML);
            break;

        case PEM:
            //-----------------------------------------------------
            //VNCSEM -> PEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCSEMmtoEMm(t, ytemp, yout, &SEML);
            break;

        case VEM:
            //-----------------------------------------------------
            //VNCSEM -> VEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCSEMmtoEMm(t, ytemp, ytemp2, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            SYSmtoSYSv(t/SEML.us_em.ns, ytemp2, yout, &SEML_EM);
            break;

        case NCSEM:
            //-----------------------------------------------------
            //VNCSEM -> NCSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, yout, &SEML);
            break;

        case PSEM:
            //-----------------------------------------------------
            //VNCEM -> PSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCtoSYS(t, ytemp, yout, &SEML);
            break;

        case VSEM:
            //-----------------------------------------------------
            //VNCEM -> VSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            NCtoSYS(t, ytemp, ytemp2, &SEML);
            SYSmtoSYSv(t, ytemp2, yout, &SEML);
            break;

        case VNCSEM:
            //-----------------------------------------------------
            //VNCSEM -> VNCSEM nothing to do
            //-----------------------------------------------------
            for(int i = 0; i <6; i++) yout[i] = y0[i];
            break;
        }
        break;

        //-----------------------------------------------------------------
        // NCSEM
        //-----------------------------------------------------------------
    case NCSEM:
        switch(outputType)
        {
        case VNCEM:
            //-----------------------------------------------------
            //NCSEM -> VNCEM
            //-----------------------------------------------------
            NCSEMmtoNCEMm(t, y0, ytemp, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            SYSmtoSYSv(t/SEML.us_em.ns, ytemp, yout, &SEML_EM);
            break;

        case NCEM:
            //-----------------------------------------------------
            //NCSEM -> NCEM
            //-----------------------------------------------------
            NCSEMmtoNCEMm(t, y0, yout, &SEML);
            break;

        case PEM:
            //-----------------------------------------------------
            //NCSEM -> PEM
            //-----------------------------------------------------
            NCSEMmtoEMm(t, y0, yout, &SEML);
            break;

        case VEM:
            //-----------------------------------------------------
            //NCSEM -> VEM
            //-----------------------------------------------------
            NCSEMmtoEMm(t, y0, ytemp, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            SYSmtoSYSv(t/SEML.us_em.ns, ytemp, yout, &SEML_EM);
            break;

        case NCSEM:
            //-----------------------------------------------------
            //NCSEM -> NCSEM nothing to do here
            //-----------------------------------------------------
            for(int i = 0; i <6; i++) yout[i] = y0[i];
            break;

        case PSEM:
            //-----------------------------------------------------
            //NCSEM -> PSEM
            //-----------------------------------------------------
            NCtoSYS(t, y0, yout, &SEML);
            break;

        case VSEM:
            //-----------------------------------------------------
            //NCSEM -> VSEM
            //-----------------------------------------------------
            NCtoSYS(t, y0, ytemp, &SEML);
            SYSmtoSYSv(t, ytemp, yout, &SEML);
            break;

        case VNCSEM:
            //-----------------------------------------------------
            //NCSEM -> VNCSEM
            //-----------------------------------------------------
            SYSmtoSYSv(t, y0, yout, &SEML);
            break;
        }
        break;

        //-----------------------------------------------------------------
        // PSEM
        //-----------------------------------------------------------------
    case PSEM:
        switch(outputType)
        {
        case VNCEM:
            //-----------------------------------------------------
            //PSEM -> VNCEM
            //-----------------------------------------------------
            SYStoNC(t, y0, ytemp, &SEML);
            NCSEMmtoNCEMm(t, ytemp, ytemp2, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            SYSmtoSYSv(t/SEML.us_em.ns, ytemp2, yout, &SEML_EM);
            break;

        case NCEM:
            //-----------------------------------------------------
            //PSEM -> NCEM
            //-----------------------------------------------------
            SYStoNC(t, y0, ytemp, &SEML);
            NCSEMmtoNCEMm(t, ytemp, yout, &SEML);
            break;

        case PEM:
            //-----------------------------------------------------
            //PSEM -> PEM
            //-----------------------------------------------------
            SEMmtoEMm(t, y0, yout, &SEML);
            break;

        case VEM:
            //-----------------------------------------------------
            //PSEM -> VEM
            //-----------------------------------------------------
            SEMmtoEMm(t, y0, ytemp, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            SYSmtoSYSv(t/SEML.us_em.ns, ytemp, yout, &SEML_EM);
            break;

        case NCSEM:
            //-----------------------------------------------------
            //PSEM -> NCSEM
            //-----------------------------------------------------
            SYStoNC(t, y0, yout, &SEML);
            break;

        case PSEM:
            //-----------------------------------------------------
            //PSEM -> PSEM nothing to do here
            //-----------------------------------------------------
            for(int i = 0; i <6; i++) yout[i] = y0[i];
            break;

        case VSEM:
            //-----------------------------------------------------
            //PSEM -> VSEM
            //-----------------------------------------------------
            SYSmtoSYSv(t, y0, yout, &SEML);
            break;

        case VNCSEM:
            //-----------------------------------------------------
            //PSEM -> VNCSEM
            //-----------------------------------------------------
            SYStoNC(t, y0, ytemp, &SEML);
            SYSmtoSYSv(t, ytemp, yout, &SEML);
            break;
        }
        break;

        //-----------------------------------------------------------------
        // VSEM
        //-----------------------------------------------------------------
    case VSEM:
        switch(outputType)
        {
        case VNCEM:
            //-----------------------------------------------------
            //VSEM -> VNCEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            SYStoNC(t, ytemp, ytemp2, &SEML);
            NCSEMmtoNCEMm(t, ytemp2, ytemp3, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            SYSmtoSYSv(t/SEML.us_em.ns, ytemp3, yout, &SEML_EM);
            break;

        case NCEM:
            //-----------------------------------------------------
            //VSEM -> NCEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            SYStoNC(t, ytemp, ytemp2, &SEML);
            NCSEMmtoNCEMm(t, ytemp2, yout, &SEML);
            break;

        case PEM:
            //-----------------------------------------------------
            //VSEM -> PEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            SEMmtoEMm(t, ytemp, yout, &SEML);
            break;

        case VEM:
            //-----------------------------------------------------
            //VSEM -> VEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            SEMmtoEMm(t, ytemp, ytemp2, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            SYSmtoSYSv(t/SEML.us_em.ns, ytemp2, yout, &SEML_EM);
            break;

        case NCSEM:
            //-----------------------------------------------------
            //VSEM -> NCSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            SYStoNC(t, ytemp, yout, &SEML);
            break;

        case PSEM:
            //-----------------------------------------------------
            //VSEM -> PSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, yout, &SEML);
            break;

        case VSEM:
            //-----------------------------------------------------
            //VSEM -> VSEM
            //-----------------------------------------------------
            for(int i = 0; i <6; i++) yout[i] = y0[i];
            break;

        case VNCSEM:
            //-----------------------------------------------------
            //VSEM -> VNCSEM
            //-----------------------------------------------------
            SYSvtoSYSm(t, y0, ytemp, &SEML);
            SYStoNC(t, ytemp, ytemp2, &SEML);
            SYSmtoSYSv(t, ytemp2, yout, &SEML);
            break;
        }
        break;
    }

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
        yvout[i] = yout[i];
    }

}
