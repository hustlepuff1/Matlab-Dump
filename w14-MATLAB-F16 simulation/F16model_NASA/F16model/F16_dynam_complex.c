/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*               6 DOF F-16 FIGHTER AIRCRAFT DYNAMICS                    */
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*                                                                       */
/* Based on the F-16 model created by R. S. Russel in                    */
/* "Nonlinear F-16 simulation using Simulink and Matlab"                 */   
/* 2003 University of Minnesota and on the F-16 model                    */
/* created by Ying Huo in "Model of F-16 Fighter Aircraft"               */
/*                                                                       */
/* Aerodynamic data and the engine model have been obtained from         */
/* "NASA Technical Paper 1538" by Nguyen et al. 1979                     */   
/*                                                                       */
/* File "F16_dynam_complex.c"                                            */
/* Version 0.9d by E.R. van Oort & L. Sonneveldt                         */
/* Created with MATLAB 7.5 R2007B                                        */
/* Feb, 2009                                                             */
/*                                                                       */
/* Difference with F16_dynam.c:                                          */
/* - incorporated cg shifts based on AIAA paper 2007-6306                */
/*   by B.J. Bacon and I.M. Gregory of NASA Langley Research Center      */
/*                                                                       */
/* Notes:                                                                */
/* -All units are SI.                                                    */
/* -Euler rotations are used.                                            */
/* -Flag is used to select between hifi and lofi aerodynamic model.      */
/* -Hifi aerodata is obtained from "aerodata/hifi_f16_aerodata.c"        */
/* -Lofi aerodata is obtained from "aerodata/lofi_f16_aerodata.c"        */
/* -"aerodata/mexndinterp.c" is used for interpolation of the data.      */                                                                    
/* -"aerodata/engine_model.c" contains the engine model.                 */  
/* -"aerodata/ISA_atmos.c" contains the ISA atmosphere model.            */
/*                                                                       */
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*  Used variables:														 */
/*																		 */
/*  Input variables: (port 1)   										 */
/*      dth         throttle setting (between 0 and 1)          [-]      */
/*      de          elevator deflection                         [rad]    */
/*      da          aileron deflection                          [rad]	 */
/*      dr          rudder deflection                           [rad]	 */
/*																		 */
/*  Additonal input variables: (port 2)									 */
/*      dlef        leading edge flap deflection                [rad]    */
/*																		 */
/*  Fidelity flag: (port 3)											     */
/*		fi_flag     hifi/lofi aerodynamic model selection flag  [-]      */
/*																		 */
/*  Center of gravity shift: (port 4)       						     */
/*		delta_x     x direction shift                           [m]      */
/*		delta_y     y direction shift                           [m]      */
/*		delta_z     z direction shift                           [m]      */
/*																		 */
/*  State variables:													 */
/*      Vt          total airspeed                              [m/s]	 */
/*      beta        angle of sideslip                           [rad]	 */
/*      alpha       angle of attack                             [rad]	 */
/*      phi         roll angle                                  [rad]    */
/*      theta       pitch angle                                 [rad]    */
/*      psi         yaw angle                                   [rad]    */
/*      p_body      roll angular rate                           [rad/s]	 */
/*      q_body      pitch angular rate                          [rad/s]	 */
/*      r_body      yaw angular rate                            [rad/s]	 */
/*      x_earth     x position                                  [m]		 */
/*      y_earth     y position                                  [m]		 */
/*      z_earth     z position                                  [m]		 */
/*      power       power level (0-100%)              [-]       [-]      */
/*																		 */
/*  State derivatives:													 */
/*      Vt_dot      rate of change in total airspeed            [m/s^2]	 */
/*      beta_dot    rate of change in angle of sideslip         [rad/s]	 */
/*      alpha_dot   rate of change in angle of attack           [rad/s]	 */
/*      phi_dot     rate of change in roll angle                [rad/s]  */
/*      theta_dot   rate of change in pitch angle               [rad/s]  */
/*      psi_dot     rate of change in yaw angle                 [rad/s]  */
/*      p_body_dot  rate of change in roll angular rate         [rad/s^2]*/
/*      q_body_dot  rate of change in pitch angular rate        [rad/s^2]*/
/*      r_body_dot  rate of change in yaw angular rate          [rad/s^2]*/
/*      x_earth_dot rate of change in x position                [m/s]	 */
/*      y_earth_dot rate of change in y position                [m/s]	 */
/*      z_earth_dot rate of change in z position                [m/s]    */
/*      power_dot   rate of change in power level               [1/s]    */
/*																		 */
/*  Output variables:												     */
/*      Vt          total airspeed                              [m/s]	 */
/*      beta        angle of sideslip                           [rad]	 */
/*      alpha       angle of attack                             [rad]	 */
/*      phi         roll angle                                  [rad]    */
/*      theta       pitch angle                                 [rad]    */
/*      psi         yaw angle                                   [rad]    */
/*      p_body      roll angular rate                           [rad/s]	 */
/*      q_body      pitch angular rate                          [rad/s]	 */
/*      r_body      yaw angular rate                            [rad/s]	 */
/*      x_earth     x position                                  [m]		 */
/*      y_earth     y position                                  [m]		 */
/*      z_earth     z position                                  [m]		 */
/*      power       power level (0-100%)              [-]       [-]      */
/*																		 */
/*  Real work vector: (Persistant memory)                    			 */
/*      C1          coefficients used in moment equations       [-]      */   
/*      C2                                                               */      
/*      C3                                                               */              
/*      C4                                                               */              
/*      C5                                                               */             
/*      C6                                                               */              
/*      C7                                                               */              
/*      C8                                                               */              
/*      C9                                                               */              
/*      Xbar         total force in body fixed x-axis           [N]      */   
/*      Ybar         total force in body fixed y-axis           [N]      */  
/*      Zbar         total force in body fixed z-axis           [N]      */  
/*      Lbar         total moment in body fixed x-axis          [Nm]     */  
/*      Mbar         total moment in body fixed y-axis          [Nm]     */
/*      Nbar         total moment in body fixed z-axis          [Nm]     */
/*      u_body       velocity in body fixed x-axis              [m/s]	 */
/*      v_body       velocity in body fixed y-axis              [m/s]	 */
/*      w_body       velocity in body fixed z-axis              [m/s]	 */
/*      u_body_dot   rate of change in velocity x-axis          [m/s^2]  */
/*      v_body_dot   rate of change in velocity y-axis          [m/s^2]  */
/*      w_body_dot   rate of change in velocity z-axis          [m/s^2]  */
/*      qbar         dynamic pressure                           [N/m]    */
/*      Mach         Mach number                                [-]      */
/*      Thrust       Total engine thrust                        [-]      */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define S_FUNCTION_NAME  F16_dynam_complex
#define S_FUNCTION_LEVEL 2
 
/* include files */
#include <math.h>
#include "simstruc.h"

#include "aerodata/mexndinterp.c"
#include "aerodata/hifi_f16_aerodata.c"     /* hifi lookup tables */
#include "aerodata/lofi_f16_aerodata.c"     /* lofi lookup tables */
#include "aerodata/matrix_functions.c"

#include "aerodata/ISA_atmos.c"             /* ISA atmosphere model */
#include "aerodata/engine_model.c"          /*engine model */

/* input port 1: control inputs */
#define dth             (*u[0])
#define de              (*u[1])
#define da              (*u[2])
#define dr              (*u[3])
/* input port 2: leading edge flap deflection */
#define dlef            (*u2[0])
/* input port 3: fidelity flag, 0 = lofi model, 1 = hifi model */
#define fi_flag         (*u3[0])

/* input port 4: gravity center shifts */
#define delta_x         (*u4[0])
#define delta_y         (*u4[1])
#define delta_z         (*u4[2])

/* 13 States */
#define Vt          x[0]
#define beta        x[1]
#define alpha       x[2]

#define phi          x[3]
#define theta        x[4]
#define psi          x[5]

#define p_body      x[6]
#define q_body      x[7]
#define r_body      x[8]

#define x_earth     x[9]
#define y_earth     x[10]
#define z_earth     x[11]

#define power         x[12]

/* State derivatives */
#define Vt_dot          dx[0]
#define beta_dot        dx[1]
#define alpha_dot       dx[2]

#define phi_dot         dx[3]
#define theta_dot       dx[4]
#define psi_dot         dx[5]

#define p_body_dot      dx[6]
#define q_body_dot      dx[7]
#define r_body_dot      dx[8]

#define x_earth_dot     dx[9]
#define y_earth_dot     dx[10]
#define z_earth_dot     dx[11]
#define power_dot         dx[12]

/* Work Variables */
#define C1              ssGetRWork(S)[0]
#define C2              ssGetRWork(S)[1]
#define C3              ssGetRWork(S)[2]
#define C4              ssGetRWork(S)[3]
#define C5              ssGetRWork(S)[4]
#define C6              ssGetRWork(S)[5]
#define C7              ssGetRWork(S)[6]
#define C8              ssGetRWork(S)[7]
#define C9              ssGetRWork(S)[8]
#define Xbar            ssGetRWork(S)[9]
#define Ybar            ssGetRWork(S)[10]
#define Zbar            ssGetRWork(S)[11]
#define Lbar            ssGetRWork(S)[12]
#define Mbar            ssGetRWork(S)[13]
#define Nbar            ssGetRWork(S)[14]
#define u_body          ssGetRWork(S)[15]
#define v_body          ssGetRWork(S)[16]
#define w_body          ssGetRWork(S)[17]
#define u_body_dot      ssGetRWork(S)[18]
#define v_body_dot      ssGetRWork(S)[19]
#define w_body_dot      ssGetRWork(S)[20]
#define qbar            ssGetRWork(S)[21]
#define Mach            ssGetRWork(S)[22]
#define Thrust          ssGetRWork(S)[23]

/* Aircraft Parameters */
#define mass        9295.44 /* assumed fixed */
#define Ixx         12874.8
#define Iyy         75673.6
#define Izz         85552.1
#define Ixz         1331.4
#define Ixy         0.0
#define Iyz         0.0
#define Sref        27.87
#define bref        9.144
#define cref        3.45
#define xcg         0.3
#define xcgr        0.35
#define heng        216.9 /* engine angular momentum, assumed fixed */

/* Additional parameters */
#define rtd           57.29577951
#define dtr           0.017453293
#define Pi			  3.141592654

/*=============================*/
/* Function: mdlInitalizeSizes */
/*=============================*/
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 1);  /* Number of expected parameters */

    ssSetNumContStates(S, 13);
    ssSetNumDiscStates(S, 0);

    ssSetNumInputPorts(S, 4);
    ssSetInputPortWidth(S, 0, 4);
    ssSetInputPortWidth(S, 1, 1);
    ssSetInputPortWidth(S, 2, 1);
    ssSetInputPortWidth(S, 3, 3);

    /* ssSetInputPortDirectFeedThrough(S, 1, 1); */

    ssSetNumOutputPorts(S, 1);
    ssSetOutputPortWidth(S, 0, 13);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 24);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetOptions(S, 0);
}

/*===================================*/
/* Function: mdlInitalizeSampleTimes */
/*===================================*/
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);

}

/*==================================*/
/* Function: mdlInitalizeConditions */
/*==================================*/

#define MDL_INITIALIZE_CONDITIONS  
#if defined(MDL_INITIALIZE_CONDITIONS)

  static void mdlInitializeConditions(SimStruct *S)
  {
  }
#endif 


/*====================*/
/* Function: mdlStart */
/*====================*/
#define MDL_START 
#if defined(MDL_START) 

  static void mdlStart(SimStruct *S)
  {
        real_T *x = ssGetContStates(S);
        int i;
        
        /* Initialize the state vector */
        for (i = 0; i < ssGetNumContStates(S); i++)
        {
            x[i] = mxGetPr(ssGetSFcnParam(S, 0))[i];    
        }
  }
#endif 



/*======================*/
/* Function: mdlOutputs */
/*======================*/
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T *x = ssGetContStates(S);
    real_T *y = ssGetOutputPortRealSignal(S, 0);
    
	 int i;
	 
	 for (i = 0; i < ssGetNumContStates(S); i++)
	 {
		y[i] = x[i]; /* outputs are the states */
	 }
}


/*=====================*/
/* Function: mdlUpdate */
/*=====================*/
#undef MDL_UPDATE  
#if defined(MDL_UPDATE)
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
  }
#endif


/*==========================*/
/* Function: mdlDerivatives */
/*==========================*/
#define MDL_DERIVATIVES  
#if defined(MDL_DERIVATIVES)
  static void mdlDerivatives(SimStruct *S)
  {
    real_T *x  = ssGetContStates(S);
    real_T *dx = ssGetdX(S);
    
    InputRealPtrsType u = ssGetInputPortRealSignalPtrs(S, 0);
    InputRealPtrsType u2 = ssGetInputPortRealSignalPtrs(S, 1);
    InputRealPtrsType u3 = ssGetInputPortRealSignalPtrs(S, 2);
    InputRealPtrsType u4 = ssGetInputPortRealSignalPtrs(S, 3);
    InputRealPtrsType u5 = ssGetInputPortRealSignalPtrs(S, 4);
    
    /* Declare a lot of variables */
    
    double temp[9] = {0};
    double Gamma;
    double CX_tot, CY_tot, CZ_tot, Cl_tot, Cm_tot, Cn_tot;
	double Cx, Cz, Cm, Cy, Cn, Cl;
    double Cxq, Cyr, Cyp, Czq, Clr, Clp, Cmq, Cnr, Cnp;
    double delta_Cx_lef, delta_Cz_lef,delta_Cm_lef;
    double delta_Cy_lef,delta_Cn_lef,delta_Cl_lef;
    double delta_Cxq_lef, delta_Cyr_lef, delta_Cyp_lef, delta_Czq_lef; 
    double delta_Clr_lef, delta_Clp_lef, delta_Cmq_lef, delta_Cnr_lef;
    double delta_Cnp_lef;
    double delta_Cy_r30, delta_Cn_r30, delta_Cl_r30;
    double delta_Cy_a20, delta_Cy_a20_lef, delta_Cn_a20, delta_Cn_a20_lef;
    double delta_Cl_a20, delta_Cl_a20_lef;
    double delta_Cnbeta, delta_Clbeta, delta_Cm, eta_el, delta_Cm_ds;
	double da_norm, dr_norm, dlef_norm;
	double ALPHA, BETA, DE, DA, DR, DLEF, alpha_lef;
    double cpow, g;
    double Fv[3][3] = {0};
    double Fw[3][3] = {0};
    double Ff[3][1] = {0};
    double Mv[3][3] = {0};
    double Mw[3][3] = {0};
    double Mf[3][1] = {0};    
    double Fv_inv[3][3] = {0};
    double Mw_inv[3][3] = {0};
    double Omega_dot[3][1] = {0};
    double VV_dot[3][1] = {0};
    double MvFv_inv[3][3] = {0};
    double MvFv_invFw[3][3] = {0};
    double inv_w[3][3] = {0};
    double MvFv_invFf[3][1] = {0};   
    double FwMw_inv[3][3] = {0};
    double FwMw_invMv[3][3] = {0};
    double inv_v[3][3] = {0};
    double FwMw_invMf[3][1] = {0};    
    double lol[3][3] = {0};
    double lol1[3][3] = {0};
    double lol2[3][3] = {0};
    double lol3[3][3] = {0};

    double Vx[3][3] = {0};
    double Dx[3][3] = {0};
    double Ox[3][3] = {0};
    double Ix[3][3] = {0};
    double OxDx[3][3] = {0};
    double OxIx[3][3] = {0};
    double VxDx[3][3] = {0};
    
    double AA[6][6] = {0};
    double BB[6][6] = {0};
    double BB_inv[6][6] = {0};
    double AAx[6][1] = {0};
    double AAt[6][1] = {0};
    double statex[6][1] = {0};
    double statex_dot[6][1] = {0};
    
    /* Usefull notation of moments of inertia,
     see e.g. Lewis and Stevens "Aircraft control and simulation" */
    Gamma = Ixx * Izz - (Ixz  * Ixz);
    C1 = ((Iyy - Izz)  * Izz  - (Ixz * Ixz))/ Gamma;
    C2 = ((Ixx - Iyy + Izz ) * Ixz ) / Gamma;
    C3 = Izz / Gamma;
    C4 = Ixz / Gamma;
    C5 = (Izz - Ixx) / Iyy;
    C6 = Ixz / Iyy ;
    C7 = 1 / Iyy;
    C8 = (Ixx * (Ixx - Iyy ) + Ixz * Ixz) / Gamma;
    C9 = Ixx / Gamma;  
    
    /* ISA atmosphere */
    atmos(-z_earth, Vt, temp);
        Mach = temp[0];
        qbar = temp[1];  
        g    = temp[2];
        
    /* Engine model, based on Ying Huo's m-files */
    cpow = tgear (dth);
    power_dot = pdot ( power, cpow );
    Thrust = thrst ( power, -z_earth, Mach ); 
    
    /* Body velocity components */
    u_body = Vt * cos(alpha) * cos(beta);
    v_body = Vt * sin(beta);
    w_body = Vt * sin(alpha) * cos(beta);
	 
    /* Transformation rad to deg for the lookup tables */
    ALPHA = alpha * rtd;
    BETA = beta * rtd;
	DE = de * rtd;
	DA = da * rtd;
	DR = dr * rtd;	 
	DLEF = dlef * rtd;	 
    
	/* LEF tables are only valid up until alpha = 45 degrees*/
	if (ALPHA > 45)
	{
	alpha_lef = 45;
	}
	else
	{
	alpha_lef = ALPHA;
	}
    
    /* Limit of elevator deflection is 25 degrees*/
	if (DE > 25)
	{
	DE = 25;
	}
	if (DE < -25)
	{
	DE = -25;
	}	 
    
    /* Normalizing the control deflections */
	da_norm = DA/21.5;
	dr_norm = DR/30.0;
      
    if (fi_flag == 1)          /* hifi lookup tables */
    {   
        
        dlef_norm = (1 - DLEF/25.0);
        
        Cx = _Cx(ALPHA, BETA, DE);
        Cy = _Cy(ALPHA, BETA);
        Cz = _Cz(ALPHA, BETA, DE);
        Cl = _Cl(ALPHA, BETA, DE);
        Cm = _Cm(ALPHA, BETA, DE);
        Cn = _Cn(ALPHA, BETA, DE);
	 
        Cxq = _CXq(ALPHA);
        Cyp = _CYp(ALPHA);
        Cyr = _CYr(ALPHA);
        Czq = _CZq(ALPHA);
        Clp = _CLp(ALPHA);
        Clr = _CLr(ALPHA);
        Cmq = _CMq(ALPHA);
        Cnp = _CNp(ALPHA);
        Cnr = _CNr(ALPHA);
	 
        delta_Cx_lef = _Cx_lef(alpha_lef, BETA) - _Cx(ALPHA, BETA, 0);
        delta_Cy_lef = _Cy_lef(alpha_lef, BETA) - _Cy(ALPHA, BETA);
        delta_Cz_lef = _Cz_lef(alpha_lef, BETA) - _Cz(ALPHA, BETA, 0);
        delta_Cl_lef = _Cl_lef(alpha_lef, BETA) - _Cl(ALPHA, BETA, 0);
        delta_Cm_lef = _Cm_lef(alpha_lef, BETA) - _Cm(ALPHA, BETA, 0);
        delta_Cn_lef = _Cn_lef(alpha_lef, BETA) - _Cn(ALPHA, BETA, 0);
	 
        delta_Cxq_lef = _delta_CXq_lef(alpha_lef);
        delta_Cyp_lef = _delta_CYp_lef(alpha_lef);
        delta_Cyr_lef = _delta_CYr_lef(alpha_lef);
        delta_Czq_lef = _delta_CZq_lef(alpha_lef);
        delta_Clp_lef = _delta_CLp_lef(alpha_lef);
        delta_Clr_lef = _delta_CLr_lef(alpha_lef);
        delta_Cmq_lef = _delta_CMq_lef(alpha_lef);
        delta_Cnp_lef = _delta_CNp_lef(alpha_lef);
        delta_Cnr_lef = _delta_CNr_lef(alpha_lef);
	 
        delta_Cy_r30 = _Cy_r30(ALPHA, BETA) - _Cy(ALPHA, BETA);
        delta_Cl_r30 = _Cl_r30(ALPHA, BETA) - _Cl(ALPHA, BETA, 0);
        delta_Cn_r30 = _Cn_r30(ALPHA, BETA) - _Cn(ALPHA, BETA, 0);
	 
        delta_Cy_a20     = _Cy_a20(ALPHA, BETA) - _Cy(ALPHA, BETA);
        delta_Cy_a20_lef = _Cy_a20_lef(alpha_lef, BETA) - 
            _Cy_lef(alpha_lef, BETA) - delta_Cy_a20;
        delta_Cn_a20     = _Cn_a20(ALPHA, BETA) - _Cn(ALPHA, BETA, 0);
        delta_Cn_a20_lef = _Cn_a20_lef(alpha_lef, BETA) - 
            _Cn_lef(alpha_lef, BETA) - delta_Cn_a20;
        delta_Cl_a20     = _Cl_a20(ALPHA, BETA) - _Cl(ALPHA, BETA, 0);
        delta_Cl_a20_lef = _Cl_a20_lef(alpha_lef, BETA) - 
            _Cl_lef(alpha_lef, BETA) - delta_Cl_a20;

        delta_Cnbeta = _delta_CNbeta(ALPHA);
        delta_Clbeta = _delta_CLbeta(ALPHA);
        delta_Cm     = _delta_Cm(ALPHA);
        eta_el       = _eta_el(DE);
        delta_Cm_ds  = _delta_Cm_ds(ALPHA,DE);
    }

else if (fi_flag == 0)  /* lofi lookup tables (do not include dlef) */
    {  
        dlef_norm = 0.0;

        damping(ALPHA, temp);
            Cxq = temp[0];
            Cyr = temp[1];
            Cyp = temp[2];
            Czq = temp[3];
            Clr = temp[4];
            Clp = temp[5];
            Cmq = temp[6];
            Cnr = temp[7];
            Cnp = temp[8];
        
        dmomdcon(ALPHA, BETA, temp);
            delta_Cl_a20 = temp[0];    
            delta_Cl_r30 = temp[1];    
            delta_Cn_a20 = temp[2];    
            delta_Cn_r30 = temp[3];    

        clcn(ALPHA, BETA, temp);
            Cl = temp[0];
            Cn = temp[1];

        cxcm(ALPHA, DE, temp);
            Cx = temp[0];
            Cm = temp[1];

            Cy = -.02*BETA + .021*da_norm + .086*dr_norm;

        cz(ALPHA, BETA, DE, temp);
            Cz = temp[0];  
        
        /* All other coeffcients are zero for the lofi model */
        delta_Cx_lef    = 0.0;
        delta_Cz_lef    = 0.0;
        delta_Cm_lef    = 0.0;
        delta_Cy_lef    = 0.0;
        delta_Cn_lef    = 0.0;
        delta_Cl_lef    = 0.0;
        delta_Cxq_lef   = 0.0;
        delta_Cyr_lef   = 0.0;
        delta_Cyp_lef   = 0.0;
        delta_Czq_lef   = 0.0;
        delta_Clr_lef   = 0.0;
        delta_Clp_lef   = 0.0;
        delta_Cmq_lef   = 0.0;
        delta_Cnr_lef   = 0.0;
        delta_Cnp_lef   = 0.0;
        delta_Cy_r30    = 0.0;
        delta_Cy_a20    = 0.0;
        delta_Cy_a20_lef= 0.0;
        delta_Cn_a20_lef= 0.0;
        delta_Cl_a20_lef= 0.0;
        delta_Cnbeta    = 0.0;
        delta_Clbeta    = 0.0;
        delta_Cm        = 0.0;
        eta_el          = 1.0; 
        delta_Cm_ds     = 0.0;
    }

    /* Total force coefficients */
    
    /* Cx_tot */
	CX_tot = Cx + delta_Cx_lef * dlef_norm 
     + (cref/(2*Vt))*(Cxq + delta_Cxq_lef * dlef_norm) * q_body;
    /* Cy_tot */
    CY_tot = Cy + delta_Cy_lef * dlef_norm 
     + (delta_Cy_a20 + delta_Cy_a20_lef * dlef_norm) * da_norm
     + delta_Cy_r30 * dr_norm 
     + (bref / (2*Vt))*(Cyr + delta_Cyr_lef * dlef_norm) * r_body 
     + (bref/(2*Vt))*(Cyp + delta_Cyp_lef * dlef_norm) * p_body;
    /* Cz_tot */
    CZ_tot = Cz + delta_Cz_lef * dlef_norm 
     + (cref/(2*Vt))*(Czq + delta_Czq_lef * dlef_norm) * q_body;
    
    /* Total moment coefficients */
    
    /* Cl_tot */
    Cl_tot = Cl + delta_Cl_lef * dlef_norm 
     + (delta_Cl_a20 + delta_Cl_a20_lef * dlef_norm) * da_norm
     + delta_Cl_r30 * dr_norm 
     + (bref / ((2*Vt))*(Clr + delta_Clr_lef * dlef_norm)) * r_body 
     + ((bref / (2*Vt)) * (Clp + delta_Clp_lef * dlef_norm)) * p_body 
     + delta_Clbeta * beta * rtd;
    /* Cm_tot */
    Cm_tot = Cm * eta_el + CZ_tot * (xcgr - xcg) + delta_Cm_lef * dlef_norm 
     + (cref / (2*Vt))*(Cmq + delta_Cmq_lef * dlef_norm) * q_body 
     + delta_Cm + delta_Cm_ds;
    /* Cn_tot */
    Cn_tot = Cn + delta_Cn_lef * dlef_norm 
     - CY_tot * (xcgr - xcg)*(cref/bref) 
     + (delta_Cn_a20 + delta_Cn_a20_lef * dlef_norm) * da_norm 
     + ((bref / (2*Vt)) * (Cnr + delta_Cnr_lef * dlef_norm))* r_body
     + ((bref / (2*Vt)) * (Cnp + delta_Cnp_lef * dlef_norm)) * p_body 
     + delta_Cn_r30 * dr_norm + delta_Cnbeta * beta * rtd;
	
    /* Total forces */
    Xbar = qbar * Sref * CX_tot;
    Ybar = qbar * Sref * CY_tot;
    Zbar = qbar * Sref * CZ_tot;
       
    /* Total moments */
    Lbar = Cl_tot * qbar * Sref * bref;
    Mbar = Cm_tot * qbar * Sref * cref;
    Nbar = Cn_tot * qbar * Sref * bref; 
       
    /*-------------------------------------------------------------------*/
    /*General Equations of Motion based on AIAA paper 2007-6306          */
    /*by B.J. Bacon and I.M. Gregory of NASA Langley Research Center     */
    /*-------------------------------------------------------------------*/
        
    Vx[0][0] = 0; Vx[0][1] = -w_body; Vx[0][2] = v_body; 
    Vx[1][0] = w_body; Vx[1][1] = 0; Vx[1][2] = -u_body; 
    Vx[2][0] = -v_body; Vx[2][1] = u_body; Vx[2][2] = 0; 
    
    Dx[0][0] = 0; Dx[0][1] = -mass*delta_z; Dx[0][2] = mass*delta_y; 
    Dx[1][0] = mass*delta_z; Dx[1][1] = 0; Dx[1][2] = -mass*delta_x; 
    Dx[2][0] = -mass*delta_y; Dx[2][1] = mass*delta_z; Dx[2][2] = 0; 
    
    Ox[0][0] = 0;Ox[0][1] = -r_body;Ox[0][2] = q_body;
    Ox[1][0] = r_body;Ox[1][1] = 0;Ox[1][2] = -p_body;
    Ox[2][0] = -q_body;Ox[2][1] = p_body;Ox[2][2] = 0;
    
    Ix[0][0] = Ixx;Ix[0][1] = -Ixy;Ix[0][2] = -Ixz;
    Ix[1][0] = -Ixy;Ix[1][1] = Iyy;Ix[1][2] = -Iyz;
    Ix[2][0] = -Ixz;Ix[2][1] = -Iyz;Ix[2][2] = Izz;
    
    mult_matrixmatrix(Ox, Dx, OxDx);
    mult_matrixmatrix(Ox, Ix, OxIx);
    mult_matrixmatrix(Vx, Dx, VxDx);
    matrix_sub(OxIx,VxDx,lol);
        
    AA[0][0] = mass*Ox[0][0]; AA[0][1] = mass*Ox[0][1]; 
    AA[0][2] = mass*Ox[0][2];
    AA[1][0] = mass*Ox[1][0]; AA[1][1] = mass*Ox[1][1]; 
    AA[1][2] = mass*Ox[1][2];
    AA[2][0] = mass*Ox[2][0]; AA[2][1] = mass*Ox[2][1]; 
    AA[2][2] = mass*Ox[2][2];
    AA[3][0] = OxDx[0][0];AA[3][1] = OxDx[0][1];AA[3][2] = OxDx[0][2];
    AA[4][0] = OxDx[1][0];AA[4][1] = OxDx[1][1];AA[4][2] = OxDx[1][2];
    AA[5][0] = OxDx[2][0];AA[5][1] = OxDx[2][1];AA[5][2] = OxDx[2][2];
    AA[0][3] = -OxDx[0][0];AA[0][4] = -OxDx[0][1];AA[0][5] = -OxDx[0][2];
    AA[1][3] = -OxDx[1][0];AA[1][4] = -OxDx[1][1];AA[1][5] = -OxDx[1][2];
    AA[2][3] = -OxDx[2][0];AA[2][4] = -OxDx[2][1];AA[2][5] = -OxDx[2][2];    
    AA[3][3] = lol[0][0];AA[3][4] = lol[0][1];AA[3][5] = lol[0][2];
    AA[4][3] = lol[1][0];AA[4][4] = lol[1][1];AA[4][5] = lol[1][2];
    AA[5][3] = lol[2][0];AA[5][4] = lol[2][1];AA[5][5] = lol[2][2]; 
    
    statex[0][0] = u_body;
    statex[1][0] = v_body;
    statex[2][0] = w_body;
    statex[3][0] = p_body;
    statex[4][0] = q_body;
    statex[5][0] = r_body;
          
    mult_matrixvector6(AA,statex,AAx);   
    
    AAt[0][0] = Xbar-AAx[0][0]-mass*g*sin(theta);
    AAt[1][0] = Ybar-AAx[1][0]+mass*g*cos(theta)*sin(phi);
    AAt[2][0] = Zbar-AAx[2][0]+mass*g*cos(theta)*cos(phi);
    AAt[3][0] = Lbar-AAx[3][0]-mass*g*(-cos(theta)*cos(phi)*delta_y
        +cos(theta)*sin(phi)*delta_z);
    AAt[4][0] = Mbar-AAx[4][0]-mass*g*(cos(theta)*cos(phi)*delta_x
        +sin(theta)*delta_z);
    AAt[5][0] = Nbar-AAx[5][0]-mass*g*(-cos(theta)*sin(phi)*delta_x
        -sin(theta)*delta_y);
    
    BB[0][0] = mass; BB[0][1] = 0; BB[0][2] = 0;
    BB[1][0] = 0; BB[1][1] = mass; BB[1][2] = 0;
    BB[2][0] = 0; BB[2][1] = 0; BB[2][2] = mass;
    BB[3][0] = Dx[0][0];BB[3][1] = Dx[0][1];BB[3][2] = Dx[0][2];
    BB[4][0] = Dx[1][0];BB[4][1] = Dx[1][1];BB[4][2] = Dx[1][2];
    BB[5][0] = Dx[2][0];BB[5][1] = Dx[2][1];BB[5][2] = Dx[2][2];
    BB[0][3] = -Dx[0][0];BB[0][4] = -Dx[0][1];BB[0][5] = -Dx[0][2];
    BB[1][3] = -Dx[1][0];BB[1][4] = -Dx[1][1];BB[1][5] = -Dx[1][2];
    BB[2][3] = -Dx[2][0];BB[2][4] = -Dx[2][1];BB[2][5] = -Dx[2][2];    
    BB[3][3] = Ix[0][0];BB[3][4] = Ix[0][1];BB[3][5] = Ix[0][2];
    BB[4][3] = Ix[1][0];BB[4][4] = Ix[1][1];BB[4][5] = Ix[1][2];
    BB[5][3] = Ix[2][0];BB[5][4] = Ix[2][1];BB[5][5] = Ix[2][2];
    
    matrix_inverse6(BB, BB_inv); 
       
    mult_matrixvector6(BB_inv,AAt,statex_dot);
    
    u_body_dot = statex_dot[0][0];
    v_body_dot = statex_dot[1][0];
    w_body_dot = statex_dot[2][0];
    
    p_body_dot  = statex_dot[3][0];
    q_body_dot  = statex_dot[4][0];
    r_body_dot  = statex_dot[5][0];
    
    /*-------------------------------------------------------------------*/
    
    /* Derivatives */

    
    Vt_dot      = (u_body * u_body_dot + v_body * v_body_dot 
     + w_body * w_body_dot) / Vt;
    beta_dot    = (v_body_dot * Vt - v_body * Vt_dot) 
     / (Vt*sqrt(u_body * u_body + w_body * w_body));
    alpha_dot   = (u_body * w_body_dot - w_body * u_body_dot) 
     / (u_body * u_body + w_body * w_body);
     
    phi_dot     = p_body+tan(theta)*(q_body*sin(phi)+r_body*cos(phi));
    theta_dot   = q_body*cos(phi)-r_body*sin(phi);
    psi_dot     = (q_body*sin(phi)+r_body*cos(phi))/cos(theta);
            
    x_earth_dot = cos(psi)*cos(theta) * u_body 
     + (cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi)) * v_body 
     + (cos(psi)*sin(theta)*cos(phi)+sin(phi)*sin(psi)) * w_body;
    y_earth_dot = sin(psi)*cos(theta) * u_body 
     + (sin(psi)*sin(theta)*sin(phi)+cos(phi)*cos(psi)) * v_body 
     + (sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi)) * w_body;
    z_earth_dot = -sin(theta) * u_body 
     + cos(theta)*sin(phi) * v_body 
     + cos(theta)*cos(phi) * w_body;   

}   /* end mdlDerivatives */

#endif 


/*===================================*/
/* Function: mdlTerminate */
/*===================================*/
static void mdlTerminate(SimStruct *S)
{
}

/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif






