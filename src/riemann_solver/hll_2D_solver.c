/**
 * @file  hll_2D_solver.c
 * @brief This is a two-dimensional HLL solver for compressible inviscid flow.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/var_struc.h"


/**
 * @brief A HLL approxiamate Riemann solver for unsteady compressible inviscid single-component flow in two space dimension.
 * @param[out] F:          All four fluxes.
 * @param[out] lambda_max: Maximum characteristic velocity.
 * @param[in] ifv_L: Left  States (rho_L, u_L, v_L, p_L, gamma, n_x, n_y).
 * @param[in] ifv_R: Right States (rho_R, u_R, v_R, p_R).
 *                   - gamma: the constant of the perfect gas.
 *                   - (n_x, n_y): unit normal vector coordinates.
 * @sa   Theory is found in Chapter 10 of Reference [1]. \n
 *       [1] E. F. Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics". 
 *           Springer-Verlag, Second Edition, 1999
 */
void HLL_2D_solver(double * F, double * lambda_max, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R)
{
	const double gamma = ifv_L->gamma;
	const double n_x   = ifv_L->n_x, n_y   = ifv_L->n_y;
	const double P_L   = ifv_L->P,   P_R   = ifv_R->P;
	const double RHO_L = ifv_L->RHO, RHO_R = ifv_R->RHO;
	const double U_L   = ifv_L->U,   U_R   = ifv_R->U;
	const double V_L   = ifv_L->V,   V_R   = ifv_R->V;

	double H_L, H_R;
	H_L = gamma/(gamma-1.0)*P_L/RHO_L + 0.5*(U_L*U_L+V_L*V_L);
	H_R = gamma/(gamma-1.0)*P_R/RHO_R + 0.5*(U_R*U_R+V_R*V_R);

	double U_S, V_S, H_S, C_S;
	U_S = (U_L*sqrt(RHO_L)+U_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	V_S = (V_L*sqrt(RHO_L)+V_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	H_S = (H_L*sqrt(RHO_L)+H_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	C_S = sqrt((gamma-1.0)*(H_S-0.5*U_S*U_S-0.5*V_S*V_S));
	double C_L, C_R;
	C_L = sqrt(gamma*P_L/RHO_L);
	C_R = sqrt(gamma*P_R/RHO_R);
	
	double qn_S, qn_L, qn_R;
	qn_S = U_S*n_x + V_S*n_y;
	qn_L = U_L*n_x + V_L*n_y;
	qn_R = U_R*n_x + V_R*n_y;

	double E_L, E_R;
	E_L = 1.0/(gamma-1.0)*P_L/RHO_L + 0.5*(U_L*U_L+V_L*V_L);
	E_R = 1.0/(gamma-1.0)*P_R/RHO_R + 0.5*(U_R*U_R+V_R*V_R);

	double S_L, S_R;
	S_L = fmin(qn_L-C_L, qn_S-C_S);
	S_R = fmax(qn_R+C_R, qn_S+C_S);
	
	S_L = fmin(0,S_L);
	S_R = fmax(0,S_R);
	
	F[0] = (S_R*RHO_L*U_L-S_L*RHO_R*U_R)*n_x + (S_R*RHO_L*V_L-S_L*RHO_R*V_R)*n_y;
	F[0] = F[0]/(S_R-S_L)+S_R*S_L/(S_R-S_L)*(RHO_R - RHO_L);
	F[1] = (S_R*RHO_L*U_L*U_L+S_R*P_L-S_L*RHO_R*U_R*U_R-S_L*P_R)*n_x + (S_R*RHO_L*U_L*V_L-S_L*RHO_R*U_R*V_R)*n_y;
	F[1] = F[1]/(S_R-S_L)+S_R*S_L/(S_R-S_L)*(RHO_R*U_R - RHO_L*U_L);
	F[2] = (S_R*RHO_L*U_L*V_L-S_L*RHO_R*U_R*V_R)*n_x + (S_R*RHO_L*V_L*V_L+S_R*P_L-S_L*RHO_R*V_R*V_R-S_L*P_R)*n_y;
	F[2] = F[2]/(S_R-S_L)+S_R*S_L/(S_R-S_L)*(RHO_R*V_R - RHO_L*V_L);
	F[3] = (S_R*RHO_L*U_L*H_L-S_L*RHO_R*U_R*H_R)*n_x + (S_R*RHO_L*V_L*H_L-S_L*RHO_R*V_R*H_R)*n_y;
	F[3] = F[3]/(S_R-S_L)+S_R*S_L/(S_R-S_L)*(RHO_R*E_R - RHO_L*E_L);

	* lambda_max = fabs(qn_S)+C_S;	  
}


