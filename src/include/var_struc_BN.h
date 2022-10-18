/**
 * @file var_struc_BN.h
 * @brief This file is the header file of some globally common variables and structural bodies
 *        for Baer-Nunziato two-phase flow model.
 */

#ifndef VARSTRUC_BN_H
#define VARSTRUC_BN_H

#ifdef DOXYGEN_PREDEFINED
/**
 * @def MULTIPHASE_BASICS
 * @brief Switch whether to compute multi-phase flow. (essential macro: MULTIFLUID_BASICS)
 */
#define MULTIPHASE_BASICS
#endif


#ifdef MULTIPHASE_BASICS
//primitive and conservative variables
struct U_var {
	double     z_s;
	double   rho_s,   u_s,   v_s,   p_s,   rho_g,   u_g,   v_g,   p_g;
	double U_rho_s, U_u_s, U_v_s, U_e_s, U_rho_g, U_u_g, U_v_g, U_e_g;
};

//Riemann invariants
struct RI_var {
	double z_s, rho_s, u_s, Q, P, H, eta_g;
};


//variables at cell centers (including staggered solid cells)
struct center_var {
	//solid cell
	double **Z_sC;
	//gaseous cell
	double **RHO_gC, **P_gC, **U_gC, **V_gC;
	double **RHO_sC, **P_sC, **U_sC, **V_sC;
	double **RHO_U_gC, **RHO_V_gC, **E_gC;
	double **RHO_U_sC, **RHO_V_sC, **E_sC;
	double **ZRHO_gC,  **ZRHO_sC;
	//gaseous cell in x/y direction
	double **Z_sS_xd, **Z_sS_yd; //S-staggered, d-direction
	double **Q_xd, **P_xd, **H_xd, **eta_g_xd;
	double **Q_yd, **P_yd, **H_yd, **eta_g_yd;
};


//slopes at centers of cells (including staggered solid cells)
struct slope_var {
	double **Z_sx, **Z_sy;
        double **RHO_sx, **RHO_sy;
	double **RHO_gx, **P_gx, **U_gx, **V_gx;
	double **RHO_gy, **P_gy, **U_gy, **V_gy;
	double **P_sx, **U_sx, **V_sx;
	double **P_sy, **U_sy, **V_sy;
	double **Z_sS_x, **Z_sS_y;
	double **Q_x, **P_x, **H_x, **eta_g_x;
	double **Q_y, **P_y, **H_y, **eta_g_y;
        double **idx; // idx=1.0 说明这里是孔隙率间断
};


//left and right variables for GRO solver (both Edir and BN-Riemann-invariant solver)
struct GRP_LR_var {
	double rho_g,  p_g,  u_g,  v_g;
	double rho_gx, p_gx, u_gx, v_gx;
	double rho_gy, p_gy, u_gy, v_gy;
	double rho_s,  p_s,  u_s,  v_s;
	double rho_sx, p_sx, u_sx, v_sx;
	double rho_sy, p_sy, u_sy, v_sy;
	double Q,  P,  H,  eta_g;
	double Qx, Px, Hx, eta_gx;
	double Qy, Py, Hy, eta_gy;
};
#endif

#endif
