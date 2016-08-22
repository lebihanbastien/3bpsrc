/**
 * \file  eminsem.cpp
 * \brief Contains all the routines to perform changes of coordinates between the EM and SEM frameworks. Including
 * \author BLB.
 * \date   2016
 * \version 1.0
 */

#include "eminsem.h"


//-----------------------------------------------------------------------------
// COC: Velocities <--> Momenta
//-----------------------------------------------------------------------------
/**
 *  \brief Change the SEM velocities into SEM momenta
 **/
void SEMvtoSEMm(double t, const double ySEv[], double ySEm[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    QBCP_L* qbp = (QBCP_L*) params_void;
    double n    =  qbp->us_sem.n;

    double delta[3];
    //evaluate the deltas
    evaluateCoef(delta, t, n, qbp->nf, qbp->cs_sem.coeffs, 3);

    //-------------------------------------------------------------------------------
    //Position to position
    //-------------------------------------------------------------------------------
    for(int i =0; i<3; i++) ySEm[i] = ySEv[i];
    //-------------------------------------------------------------------------------
    //Velocity to Momenta
    //-------------------------------------------------------------------------------
    //px
    ySEm[3] = 1.0/delta[0] * (ySEv[3] - delta[1]*ySEv[0] - delta[2]*ySEv[1]);
    //py
    ySEm[4] = 1.0/delta[0] * (ySEv[4] - delta[1]*ySEv[1] + delta[2]*ySEv[0]);
    //pz
    ySEm[5] = 1.0/delta[0] * (ySEv[5] - delta[1]*ySEv[2] );
}

/**
 *  \brief Change the SEM momenta into SEM velocities
 **/
void SEMmtoSEMv(double t, const double ySEm[], double ySEv[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    QBCP_L* qbp = (QBCP_L*) params_void;
    double n    =  qbp->us_sem.n;

    double delta[3];
    //evaluate the deltas
    evaluateCoef(delta, t, n, qbp->nf, qbp->cs_sem.coeffs, 3);

    //-------------------------------------------------------------------------------
    //Position to position
    //-------------------------------------------------------------------------------
    for(int i =0; i<3; i++) ySEv[i] = ySEm[i];

    //-------------------------------------------------------------------------------
    //Momenta to velocities
    //-------------------------------------------------------------------------------
    //vx
    ySEv[3] = delta[0]*ySEm[3] + delta[1]*ySEm[0] + delta[2]*ySEm[1];
    //vy
    ySEv[4] = delta[0]*ySEm[4] + delta[1]*ySEm[1] - delta[2]*ySEm[0];
    //vz
    ySEv[5] = delta[0]*ySEm[5] + delta[1]*ySEm[2];
}

/**
 *  \brief Change the EM velocities into EM momenta
 **/
void EMvtoEMm(double t, const double yEMv[], double yEMm[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    QBCP_L* qbp  = (QBCP_L*) params_void;
    double n     =  qbp->us_em.n;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[3];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 3);

    //-------------------------------------------------------------------------------
    //Position to position
    //-------------------------------------------------------------------------------
    for(int i =0; i<3; i++) yEMm[i] = yEMv[i];


    //-------------------------------------------------------------------------------
    //Velocity to Momenta
    //-------------------------------------------------------------------------------
    //px
    yEMm[3] = 1.0/alpha[0] * (yEMv[3] - alpha[1]*yEMv[0] - alpha[2]*yEMv[1]);
    //py
    yEMm[4] = 1.0/alpha[0] * (yEMv[4] - alpha[1]*yEMv[1] + alpha[2]*yEMv[0]);
    //pz
    yEMm[5] = 1.0/alpha[0] * (yEMv[5] - alpha[1]*yEMv[2] );
}

/**
 *  \brief Change the EM momenta into EM velocities
 **/
void EMmtoEMv(double t, const double yEMm[], double yEMv[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    QBCP_L* qbp  = (QBCP_L*) params_void;
    double n     =  qbp->us_em.n;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[3];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 3);

    //-------------------------------------------------------------------------------
    //Position to position
    //-------------------------------------------------------------------------------
    for(int i =0; i<3; i++) yEMv[i] = yEMm[i];

    //-------------------------------------------------------------------------------
    //Momenta to velocities
    //-------------------------------------------------------------------------------
    //vx
    yEMv[3] = alpha[0]*yEMm[3] + alpha[1]*yEMm[0] + alpha[2]*yEMm[1];
    //vy
    yEMv[4] = alpha[0]*yEMm[4] + alpha[1]*yEMm[1] - alpha[2]*yEMm[0];
    //vz
    yEMv[5] = alpha[0]*yEMm[5] + alpha[1]*yEMm[2];
}

/**
 *  \brief Change the EM velocities into EM momenta
 **/
void SYSvtoSYSm(double t, const double ySYSv[], double ySYSm[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    QBCP_L* qbp  = (QBCP_L*) params_void;
    double n     =  qbp->us.n;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[3];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, 3);

    //-------------------------------------------------------------------------------
    //Position to position
    //-------------------------------------------------------------------------------
    for(int i =0; i<3; i++) ySYSm[i] = ySYSv[i];


    //-------------------------------------------------------------------------------
    //Velocity to Momenta
    //-------------------------------------------------------------------------------
    //px
    ySYSm[3] = 1.0/alpha[0] * (ySYSv[3] - alpha[1]*ySYSv[0] - alpha[2]*ySYSv[1]);
    //py
    ySYSm[4] = 1.0/alpha[0] * (ySYSv[4] - alpha[1]*ySYSv[1] + alpha[2]*ySYSv[0]);
    //pz
    ySYSm[5] = 1.0/alpha[0] * (ySYSv[5] - alpha[1]*ySYSv[2] );
}

/**
 *  \brief Change the EM momenta into EM velocities
 **/
void SYSmtoSYSv(double t, const double ySYSm[], double ySYSv[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    QBCP_L* qbp  = (QBCP_L*) params_void;
    double n     =  qbp->us.n;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[3];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, 3);

    //-------------------------------------------------------------------------------
    //Position to position
    //-------------------------------------------------------------------------------
    for(int i =0; i<3; i++) ySYSv[i] = ySYSm[i];

    //-------------------------------------------------------------------------------
    //Momenta to velocities
    //-------------------------------------------------------------------------------
    //vx
    ySYSv[3] = alpha[0]*ySYSm[3] + alpha[1]*ySYSm[0] + alpha[2]*ySYSm[1];
    //vy
    ySYSv[4] = alpha[0]*ySYSm[4] + alpha[1]*ySYSm[1] - alpha[2]*ySYSm[0];
    //vz
    ySYSv[5] = alpha[0]*ySYSm[5] + alpha[1]*ySYSm[2];
}


//-----------------------------------------------------------------------------
// Change of unit system
//-----------------------------------------------------------------------------
/**
 *   \brief From SEM unit system to EM unit system for a Position/Velocity and time vector in IN coordinates and time
 **/
void ussem2usem(double *tc, double yINv[], QBCP_L *qbcp_l)
{
    //IN[SEM] to IN[EM]
    yINv[0] *= qbcp_l->us_em.as;
    yINv[1] *= qbcp_l->us_em.as;
    yINv[2] *= qbcp_l->us_em.as;
    yINv[3] *= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    yINv[4] *= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    yINv[5] *= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    *tc     /= qbcp_l->us_em.ns;
}

/**
 *   \brief From EM unit system to SEM unit system for a Position/Velocity and time vector in IN coordinates and time
 **/
void usem2ussem(double *tc, double yINv[], QBCP_L *qbcp_l)
{
    //IN[EM] to IN[SEM]
    yINv[0] /= qbcp_l->us_em.as;
    yINv[1] /= qbcp_l->us_em.as;
    yINv[2] /= qbcp_l->us_em.as;
    yINv[3] /= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    yINv[4] /= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    yINv[5] /= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    *tc     *= qbcp_l->us_em.ns;
}

//-----------------------------------------------------------------------------
// COC: IN <--> EM
//-----------------------------------------------------------------------------
/**
 * \brief From EM to IN (in EM units)
 **/
void EMtoIN(double t, const double yEM[], double yIN[], QBCP_L *qbcp_l)
{
    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    //Param
    double n  = qbcp_l->us_em.n;
    double ms = qbcp_l->us_em.ms;
    double ns = qbcp_l->us_em.ns;
    double as = qbcp_l->us_em.as;
    double ni = qbcp_l->us_em.ni;
    double ai = qbcp_l->us_em.ai;
    //r
    double r1 = creal(evz(qbcp_l->cs_em.zt, t, n, ni, ai));
    double r2 = cimag(evz(qbcp_l->cs_em.zt, t, n, ni, ai));
    double r  = sqrt(r1*r1 + r2*r2);

    //R
    double R1 = creal(evz(qbcp_l->cs_em.Zt, t, n, ns, as));
    double R2 = cimag(evz(qbcp_l->cs_em.Zt, t, n, ns, as));
    //rdot
    double r1dot = creal(evzdot(qbcp_l->cs_em.zt, qbcp_l->cs_em.ztdot, t, n, ni, ai));
    double r2dot = cimag(evzdot(qbcp_l->cs_em.zt, qbcp_l->cs_em.ztdot, t, n, ni, ai));
    double rdot  = 1.0/r*(r1*r1dot + r2*r2dot);
    //Rdot
    double R1dot = creal(evzdot(qbcp_l->cs_em.Zt, qbcp_l->cs_em.Ztdot, t, n, ns, as));
    double R2dot = cimag(evzdot(qbcp_l->cs_em.Zt, qbcp_l->cs_em.Ztdot, t, n, ns, as));

    //-------------------------------------------------------------------------------
    // EM to IN: Position
    //-------------------------------------------------------------------------------
    yIN[0] = r1*yEM[0] - r2*yEM[1] - ms/(1.0+ms)*R1;
    yIN[1] = r2*yEM[0] + r1*yEM[1] - ms/(1.0+ms)*R2;
    yIN[2] = r *yEM[2];

    //-------------------------------------------------------------------------------
    // EM to IN: Velocity
    //-------------------------------------------------------------------------------
    yIN[3] = r1dot*yEM[0] - r2dot*yEM[1] + r1*yEM[3] - r2*yEM[4]- ms/(1.0+ms)*R1dot;
    yIN[4] = r2dot*yEM[0] + r1dot*yEM[1] + r2*yEM[3] + r1*yEM[4]- ms/(1.0+ms)*R2dot;
    yIN[5] = rdot *yEM[2] + r *yEM[5];
}

/**
 * \brief From IN to EM (in EM units)
 **/
void INtoEM(double t, const double yIN[], double yEM[],
                                          QBCP_L *qbcp_l)
{
    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    //Param
    double n  = qbcp_l->us_em.n;
    double ms = qbcp_l->us_em.ms;
    double ns = qbcp_l->us_em.ns;
    double as = qbcp_l->us_em.as;
    double ni = qbcp_l->us_em.ni;
    double ai = qbcp_l->us_em.ai;

    //r
    double r1 = creal(evz(qbcp_l->cs_em.zt, t, n, ni, ai));
    double r2 = cimag(evz(qbcp_l->cs_em.zt, t, n, ni, ai));
    double r = sqrt(r1*r1 + r2*r2);
    //R
    double R1 = creal(evz(qbcp_l->cs_em.Zt, t, n, ns, as));
    double R2 = cimag(evz(qbcp_l->cs_em.Zt, t, n, ns, as));
    //rdot
    double r1dot = creal(evzdot(qbcp_l->cs_em.zt, qbcp_l->cs_em.ztdot, t, n, ni, ai));
    double r2dot = cimag(evzdot(qbcp_l->cs_em.zt, qbcp_l->cs_em.ztdot, t, n, ni, ai));
    //Rdot
    double R1dot = creal(evzdot(qbcp_l->cs_em.Zt, qbcp_l->cs_em.Ztdot, t, n, ns, as));
    double R2dot = cimag(evzdot(qbcp_l->cs_em.Zt, qbcp_l->cs_em.Ztdot, t, n, ns, as));

    //Additional parameters
    double a = +pow(r, -4.0)*(r1dot*r*r - 2*r1*(r1*r1dot + r2*r2dot)); //dot(r1/(r*r))
    double b = +pow(r, -4.0)*(r2dot*r*r - 2*r2*(r1*r1dot + r2*r2dot)); //dot(r2/(r*r))
    double c = -pow(r, -3.0)*(r1*r1dot + r2*r2dot);                    //dot(1/r)

    //-------------------------------------------------------------------------------
    // EM to IN: Position
    //-------------------------------------------------------------------------------
    yEM[0] = 1.0/(r*r) * ( +r1*(yIN[0] + ms/(1.0+ms)*R1) +  r2*(yIN[1] + ms/(1.0+ms)*R2) );
    yEM[1] = 1.0/(r*r) * ( -r2*(yIN[0] + ms/(1.0+ms)*R1) +  r1*(yIN[1] + ms/(1.0+ms)*R2) );
    yEM[2] = 1.0/r * yIN[2];

    //-------------------------------------------------------------------------------
    // EM to IN: Velocity
    //-------------------------------------------------------------------------------
    yEM[3] = +a*(yIN[0] + ms/(1.0+ms)*R1) + b*(yIN[1] + ms/(1.0+ms)*R2)
             + 1.0/(r*r) * ( +r1*(yIN[3] + ms/(1.0+ms)*R1dot) +  r2*(yIN[4] + ms/(1.0+ms)*R2dot) );
    yEM[4] = -b*(yIN[0] + ms/(1.0+ms)*R1) + a*(yIN[1] + ms/(1.0+ms)*R2)
             + 1.0/(r*r) * ( -r2*(yIN[3] + ms/(1.0+ms)*R1dot) +  r1*(yIN[4] + ms/(1.0+ms)*R2dot) );
    yEM[5] = c*yIN[2] + 1.0/r * yIN[5];
}

//-----------------------------------------------------------------------------
// COC: IN <--> SEM
//-----------------------------------------------------------------------------
/**
 * \brief From SEM to IN (in SEM units)
 **/
void SEMtoIN(double t, const double ySE[], double yIN[],
                                          QBCP_L *qbcp_l)
{
    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    //Param
    double n  = qbcp_l->us_sem.n;
    double ns = qbcp_l->us_sem.ns;
    double as = qbcp_l->us_sem.as;

    //-------------------------------------------------------------------------------
    //r & R
    //-------------------------------------------------------------------------------
    //R
    double R1 = creal(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    double R2 = cimag(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    double R = sqrt(R1*R1 + R2*R2);

    //-------------------------------------------------------------------------------
    //Derivatives
    //-------------------------------------------------------------------------------
    //Rdot
    double R1dot = creal(evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as));
    double R2dot = cimag(evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as));
    double Rdot  = 1.0/R*(R1*R1dot + R2*R2dot);

    //-------------------------------------------------------------------------------
    //Position & Velocity of the SE barycenter in inertial coordinates (in SE units)
    //-------------------------------------------------------------------------------
    double yesIN[6];
    yesIN[0] = 0.0;
    yesIN[1] = 0.0;
    yesIN[2] = 0.0;
    yesIN[3] = 0.0;
    yesIN[4] = 0.0;
    yesIN[5] = 0.0;


    //-------------------------------------------------------------------------------
    // SE to IN: Position
    //-------------------------------------------------------------------------------
    yIN[0] = R1*ySE[0] - R2*ySE[1] + yesIN[0];
    yIN[1] = R2*ySE[0] + R1*ySE[1] + yesIN[1];
    yIN[2] = R *ySE[2];

    //-------------------------------------------------------------------------------
    // SE to IN: Velocity
    //-------------------------------------------------------------------------------
    yIN[3] = R1dot*ySE[0] - R2dot*ySE[1] + R1*ySE[3] - R2*ySE[4] + yesIN[3];
    yIN[4] = R2dot*ySE[0] + R1dot*ySE[1] + R2*ySE[3] + R1*ySE[4] + yesIN[4];
    yIN[5] = Rdot *ySE[2] + R *ySE[5];
}

/**
 * \brief From IN to SEM (in SEM units)
 **/
void INtoSEM(double t, const double yIN[], double ySE[],
                                          QBCP_L *qbcp_l)
{
    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    //Param
    double n  = qbcp_l->us_sem.n;
    double ns = qbcp_l->us_sem.ns;
    double as = qbcp_l->us_sem.as;

    //-------------------------------------------------------------------------------
    //r & R
    //-------------------------------------------------------------------------------
    //R
    double R1 = creal(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    double R2 = cimag(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    //h
    double h1 = R1;
    double h2 = R2;
    double h  = sqrt(h1*h1 + h2*h2);

    //-------------------------------------------------------------------------------
    //Derivatives
    //-------------------------------------------------------------------------------
    //Rdot
    double R1dot = creal(evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as));
    double R2dot = cimag(evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as));
    //hdot
    double h1dot = R1dot;
    double h2dot = R2dot;

    //-------------------------------------------------------------------------------
    //Position & Velocity of the SE barycenter in inertial coordinates and SE units
    //-------------------------------------------------------------------------------
    double yesIN[6];
    yesIN[0] = 0.0;
    yesIN[1] = 0.0;
    yesIN[2] = 0.0;
    yesIN[3] = 0.0;
    yesIN[4] = 0.0;
    yesIN[5] = 0.0;

    //Additional parameters
    double a = +pow(h, -4.0)*(h1dot*h*h - 2*h1*(h1*h1dot + h2*h2dot)); //dot(h1/(h*h))
    double b = +pow(h, -4.0)*(h2dot*h*h - 2*h2*(h1*h1dot + h2*h2dot)); //dot(h2/(h*h))
    double c = -pow(h, -3.0)*(h1*h1dot  + h2*h2dot);                   //dot(1/h)

    //-------------------------------------------------------------------------------
    // SE to IN: Position
    //-------------------------------------------------------------------------------
    ySE[0] = 1.0/(h*h) * ( +h1*(yIN[0] - yesIN[0]) +  h2*(yIN[1] - yesIN[1]) );
    ySE[1] = 1.0/(h*h) * ( -h2*(yIN[0] - yesIN[0]) +  h1*(yIN[1] - yesIN[1]) );
    ySE[2] = 1.0/h * yIN[2];

    //-------------------------------------------------------------------------------
    // SE to IN: Velocity
    //-------------------------------------------------------------------------------
    ySE[3] = +a*(yIN[0] - yesIN[0]) + b*(yIN[1] - yesIN[1])
             + 1.0/(h*h) * ( +h1*(yIN[3] - yesIN[3]) +  h2*(yIN[4] - yesIN[4]) );
    ySE[4] = -b*(yIN[0] - yesIN[0]) + a*(yIN[1] - yesIN[1])
             + 1.0/(h*h) * ( -h2*(yIN[3] - yesIN[3]) +  h1*(yIN[4] - yesIN[4]) );
    ySE[5] = +c*yIN[2] + 1.0/h * yIN[5];
}

//-----------------------------------------------------------------------------
// COC: SEM <--> IN <--> EM
//-----------------------------------------------------------------------------
/**
 * \brief From SEM to EM (both in position/momenta form)
 **/
void SEMmtoEMm(double t, const double ySEm[], double yEMm[],
              QBCP_L *qbcp_l)
{
    double tc = t;

    //Momenta to velocities
    double ySEv[42];
    SEMmtoSEMv(tc, ySEm, ySEv, qbcp_l);

    //SE --> IN
    double yINv[42];
    SEMtoIN(tc, ySEv, yINv, qbcp_l);

    //IN[SE] to IN[SEM]
    ussem2usem(&tc, yINv, qbcp_l);

    //IN --> EM
    double yEMv[42];
    INtoEM(tc, yINv, yEMv, qbcp_l);

    //Velocities to momenta
    EMvtoEMm(tc, yEMv, yEMm, qbcp_l);
}

/**
 * \brief From EM to SEM (both in position/momenta form)
 **/
void EMmtoSEMm(double t, const double yEMm[], double ySEMm[],
              QBCP_L *qbcp_l)
{
    double tc = t;
    //Momenta to velocities
    double yEMv[42];
    EMmtoEMv(tc, yEMm, yEMv, qbcp_l);

    //EM-->IN
    double yINv[42];
    EMtoIN(tc, yEMv, yINv, qbcp_l);

    //IN[EM] to IN[SEM]
    usem2ussem(&tc, yINv, qbcp_l);

    //IN-->SE
    double ySEMv[42];
    INtoSEM(tc, yINv, ySEMv, qbcp_l);

    //Velocities to momenta
    SEMvtoSEMm(tc, ySEMv, ySEMm, qbcp_l);
}

/**
 * \brief From NC EM to SEM (both in position/momenta form)
 **/
void NCEMmtoSEMm(double t, const double yNCEMm[], double ySEMm[], QBCP_L *qbcp_l)
{
    double yEMm[42];
    //NC to EM
    NCtoEM(t, yNCEMm, yEMm, qbcp_l);
    //EM to SEM
    EMmtoSEMm(t, yEMm, ySEMm, qbcp_l);
}

/**
 * \brief From NC EM to  NC SEM (both in position/momenta form)
 **/
void NCEMmtoNCSEMm(double tEM, const double yNCEMm[], double yNCSEM[], QBCP_L *qbcp_l)
{
    double yEMm[42], ySEMm[42];
    //NC to EM
    NCtoEM(tEM, yNCEMm, yEMm, qbcp_l);
    //EM to SEM
    EMmtoSEMm(tEM, yEMm, ySEMm, qbcp_l);
    //SEM to NC (careful, the time should be set into SEM units!)
    SEMtoNC(tEM*qbcp_l->us_em.ns, ySEMm, yNCSEM, qbcp_l);
}

/**
 * \brief From NC SEM to  NC EM (both in position/momenta form)
 **/
void NCSEMmtoNCEMm(double t, const double yNCSEMm[], double yNCEM[], QBCP_L *qbcp_l)
{
    double ySEMm[42], yEMm[42];
    //NC to EM
    NCtoSEM(t, yNCSEMm, ySEMm, qbcp_l);
    //EM to SEM
    SEMmtoEMm(t, ySEMm, yEMm, qbcp_l);
    //EM to NC (careful, the time should be set into EM units!)
    EMtoNC(t/qbcp_l->us_em.ns, yEMm, yNCEM, qbcp_l);
}

/**
 * \brief From NC SEM to EM (both in position/momenta form)
 **/
void NCSEMmtoEMm(double t, const double yNCSEMm[], double yEMm[], QBCP_L *qbcp_l)
{
    double ySEMm[42];
    //NC to EM
    NCtoSEM(t, yNCSEMm, ySEMm, qbcp_l);
    //EM to SEM
    SEMmtoEMm(t, ySEMm, yEMm, qbcp_l);
}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Change of coordinates on vectors
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief From NC EM to SEM for vectors. Note that the time vector is left unchanged (still in EM units).
 **/
void NCEMmtoSEMm_vec(double **yNCEM, double *tNCEM, double **ySEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to SEM
        NCEMmtoSEMm(tNCEM[p], yNCd, ySEMd, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySEM[k][p] = ySEMd[k];
    }
}

/**
 * \brief From NC EM to SEM for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoSEMm_vec(double **yNCEM, double *tNCEM, double **ySEM, double *tSEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to SEM
        NCEMmtoSEMm(tNCEM[p], yNCd, ySEMd, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySEM[k][p] = ySEMd[k];
        //time
        tSEM[p] = tNCEM[p]*qbcp_l->us_em.ns;
    }
}

/**
 * \brief From NC EM to SEMv for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoSEMv_vec(double **yNCEM, double *tNCEM, double **ySEM, double *tSEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySEMm[6], ySEMv[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Time
        tSEM[p] = tNCEM[p]*qbcp_l->us_em.ns;
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to SEM
        NCEMmtoSEMm(tNCEM[p], yNCd, ySEMm, qbcp_l);
        //SEM to SEMv
        SEMmtoSEMv(tSEM[p], ySEMm, ySEMv, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySEM[k][p] = ySEMv[k];
    }
}

/**
 * \brief From NC EM to NC SEM for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoNCSEMm_vec(double **yNCEM, double *tNCEM, double **yNCSEM, double *tSEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yNCSEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to NC SEM
        NCEMmtoNCSEMm(tNCEM[p], yNCd, yNCSEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yNCSEM[k][p] = yNCSEMd[k];
        //time
        tSEM[p] = tNCEM[p]*qbcp_l->us_em.ns;
    }
}

/**
 * \brief From NC EM to NC SEM for vectors. Note that the time vector is left unchanged (still in EM units).
 **/
void NCEMmtoNCSEMm_vec(double **yNCEM, double *tNCEM, double **yNCSEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yNCSEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to NC SEM
        NCEMmtoNCSEMm(tNCEM[p], yNCd, yNCSEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yNCSEM[k][p] = yNCSEMd[k];
    }
}

/**
 * \brief From NC SEM to EM for vectors. Note that the time vector is left unchanged (still in SEM units).
 **/
void NCSEMmtoEMm_vec(double **yNCSEM, double *tNCSEM, double **yEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCSEM[k][p];
        //NC EM to NC SEM
        NCSEMmtoEMm(tNCSEM[p], yNCd, yEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yEM[k][p] = yEMd[k];
    }
}

/**
 * \brief From NC SEM to EM for vectors. The time vector tNCSEM is transformed into EM units and stored into tEM.
 **/
void NCSEMmtoEMm_vec(double **yNCSEM, double *tNCSEM, double **yEM, double *tEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCSEM[k][p];
        //NC EM to NC SEM
        NCSEMmtoEMm(tNCSEM[p], yNCd, yEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yEM[k][p] = yEMd[k];
        //time
        tEM[p] = tNCSEM[p]/qbcp_l->us_em.ns;
    }
}

/**
 * \brief From NC SEM to EMv for vectors. The time vector tNCSEM is transformed into EM units and stored into tEM.
 **/
void NCSEMmtoEMv_vec(double **yNCSEM, double *tNCSEM, double **yEM, double *tEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yEMm[6], yEMv[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Time
        tEM[p] = tNCSEM[p]/qbcp_l->us_em.ns;
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCSEM[k][p];
        //NC SEM to EM
        NCSEMmtoEMm(tNCSEM[p], yNCd, yEMm, qbcp_l);
        //EM to EMv
        EMmtoEMv(tEM[p], yEMm, yEMv, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yEM[k][p] = yEMv[k];
    }
}



/**
 * \brief From NC SEM to NC EM for vectors. Note that the time vector is left unchanged (still in SEM units).
 **/
void NCSEMmtoNCEMm_vec(double **yNCSEM, double *tNCSEM, double **yNCEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yNCEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCSEM[k][p];
        //NC EM to NC SEM
        NCSEMmtoNCEMm(tNCSEM[p], yNCd, yNCEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yNCEM[k][p] = yNCEMd[k];
    }
}

/**
 * \brief From NC SEM to NC EM for vectors. The time vector tNCSEM is transformed into EM units and stored into tEM.
 **/
void NCSEMmtoNCEMm_vec(double **yNCSEM, double *tNCSEM, double **yNCEM, double *tEM,  int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yNCEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCSEM[k][p];
        //NC EM to NC SEM
        NCSEMmtoNCEMm(tNCSEM[p], yNCd, yNCEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yNCEM[k][p] = yNCEMd[k];
        //time
        tEM[p] = tNCSEM[p]/qbcp_l->us_em.ns;
    }
}

/**
 *  \brief From NC to SYS, in SYS units (either EM or SEM).
 **/
void NCtoSYS_vec(double **yNC, double *tNC, double **ySYS, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySYSd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNC[k][p];
        //NC EM to SEM
        NCtoSYS(tNC[p], yNCd, ySYSd, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySYS[k][p] = ySYSd[k];
    }
}

/**
 *  \brief From NC to SYSv, in SYS units (either EM or SEM).
 **/
void NCtoSYSv_vec(double **yNC, double *tNC, double **ySYS, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySYSm[6], ySYSv[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNC[k][p];
        //NC to SYS
        NCtoSYS(tNC[p], yNCd, ySYSm, qbcp_l);
        //SYSm to SYSv
        SYSmtoSYSv(tNC[p], ySYSm, ySYSv, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySYS[k][p] = ySYSv[k];
    }
}

/**
 *  \brief From SYS to NC, in SYS units (either EM or SEM).
 **/
void SYStoNC_vec(double **ySYS, double *tSYS, double **yNC, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySYSd[6];
    //Loop on all elements in ySYS
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in ySYS
        for(int k = 0; k <6; k++) ySYSd[k] = ySYS[k][p];
        //SYS to NC
        SYStoNC(tSYS[p], ySYSd, yNCd, qbcp_l);
        //Copy result in yNC
        for(int k = 0; k <6; k++) yNC[k][p] = yNCd[k];
    }
}

/**
 *  \brief From NC to SEM, in SEM units.
 **/
void NCtoSEM_vec(double **yNC, double *tNC, double **ySYS, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySYSd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNC[k][p];
        //NC EM to SEM
        NCtoSEM(tNC[p], yNCd, ySYSd, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySYS[k][p] = ySYSd[k];
    }
}

/**
 *  \brief From NC to VNC, in NC units (either EM or SEM).
 **/
void NCtoVNC_vec(double **yNC, double *tNC, double **yVNC, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yNCv[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNC[k][p];
        //NC to VNC
        SYSmtoSYSv(tNC[p], yNCd, yNCv, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) yVNC[k][p] = yNCv[k];
    }
}

/**
 *  \brief From VNC to NC, in NC units (either EM or SEM).
 **/
void VNCtoNC_vec(double **yVNC, double *tNC, double **yNC, int N, QBCP_L *qbcp_l)
{
    double yNCv[6], yNCm[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCv[k] = yVNC[k][p];
        //VNC to NC
        SYSvtoSYSm(tNC[p], yNCv, yNCm, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) yNC[k][p] = yNCm[k];
    }
}

//------------------------------------------------------------------------------------
//   Precisions
//------------------------------------------------------------------------------------
/**
 *  \brief Sets an average precision in cout.
 **/
void coutmp()
{
    cout <<  setw(5) << setprecision(5) << std::showpos  <<  setiosflags(ios::scientific);
}

/**
 *  \brief Sets a big precision in cout.
 **/
void coutlp()
{
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);
}

/**
 *  \brief Sets a small precision in cout.
 **/
void coutsp()
{
    cout <<  setw(5) << setprecision(3) << resetiosflags(ios::scientific);
}

//------------------------------------------------------------------------------------
//   Print
//------------------------------------------------------------------------------------
/**
 *  \brief Prints an array of double using cout.
 **/
void vector_printf(double *y, int n)
{
    coutmp();
    for(int i = 0; i < n; i ++) cout << i << "  " << y[i] << endl;
}

/**
 *  \brief Prints an array of complex double using cout.
 **/
void vector_complex_printf(cdouble *y, int n)
{
    coutmp();
    for(int i = 0; i < n; i ++) cout << i << "  " << creal(y[i]) << "  " << cimag(y[i]) << endl;
}
