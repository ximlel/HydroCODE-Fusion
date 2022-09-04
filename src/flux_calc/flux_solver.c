/**
 * @file flux_solver.c
 * @brief This file is a set of functions to calculate interfacial fluxes and demanded variables
 *        according to the left and right state of the cell interface by certain solver.
 */
#include <stdio.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/riemann_solver.h"


void Roe_flux(struct i_f_var * ifv, struct i_f_var * ifv_R)
{
	const int dim = (int)config[0];
	const double delta = 0.2;

	double F[4];
	double lambda_max;
	if (dim == 1)
		{
			Roe_solver(F, &lambda_max, ifv, ifv_R, delta);
			ifv->F_rho = F[0];
			ifv->F_u   = F[1];
			ifv->F_e   = F[2];
		}
	else if (dim == 2)
		{
			Roe_2D_solver(F, &lambda_max, ifv, ifv_R, delta);
			ifv->F_rho = F[0];
			ifv->F_u   = F[1];
			ifv->F_v   = F[2];
			ifv->F_e   = F[3];
		}
}


void HLL_flux(struct i_f_var * ifv, struct i_f_var * ifv_R)
{
	double F[4];
	double lambda_max;
	HLL_2D_solver(F, &lambda_max, ifv, ifv_R);
	ifv->F_rho = F[0];
	ifv->F_u   = F[1];
	ifv->F_v   = F[2];
	ifv->F_e   = F[3];
}


/**
 * @brief This function calculate Eulerian fluxes of 2-D Euler equations by Riemann solver.
 * @param[in,out] ifv: Structure pointer of interfacial evaluated variables and fluxes and left state.
 * @param[in] ifv_R:   Structure pointer of interfacial right state.
 * @param[in] tau:     The length of the time step.
 * @return    miscalculation indicator.
 *   @retval  0: Successful calculation.
 *   @retval  1: < 0.0 error.
 *   @retval  2: NAN or INFinite error of mid[].
 */
int Riemann_exact_flux(struct i_f_var * ifv, struct i_f_var * ifv_R)
{
	const int dim = (int)config[0];
	const double eps = config[4];
	const double n_x = ifv->n_x, n_y = ifv->n_y;
	double gamma_mid = ifv->gamma;
	ifv->lambda_u = 0.0;  ifv->lambda_v = 0.0;

	if (dim == 2)
		{
			double u, u_R;
			u        =  ifv->U  *n_x + ifv->V  *n_y;
			u_R      =  ifv_R->U*n_x + ifv_R->V*n_y;
			ifv->V   = -ifv->U  *n_y + ifv->V  *n_x;
			ifv_R->V = -ifv_R->U*n_y + ifv_R->V*n_x;
			ifv->U   =  u;
			ifv_R->U =  u_R;
		}

	double wave_speed[2], dire[6], mid[6], star[6];

	linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, ifv, ifv_R, eps, INFINITY);

	if(mid[3] < eps || mid[0] < eps)
	    return 1;
	if(!isfinite(mid[1])|| !isfinite(mid[2])|| !isfinite(mid[0])|| !isfinite(mid[3]))
	    return 2;

	double rho_mid = mid[0], p_mid = mid[3], u_mid = mid[1], v_mid;
#ifdef MULTIFLUID_BASICS
	double phi_mid = mid[4], z_a_mid = mid[5];
	gamma_mid = mid[1] > 0.0 ? ifv->gamma : ifv_R->gamma;
#endif
	if (dim == 1)
		{
			u_mid = mid[1];
			ifv->F_rho = rho_mid*u_mid;
			ifv->F_u   = ifv->F_rho*u_mid + p_mid;
		}
	if (dim == 2)
		{
			u_mid  = mid[1]*n_x - mid[2]*n_y;
			v_mid  = mid[1]*n_y + mid[2]*n_x;
			ifv->F_rho = rho_mid*(u_mid*n_x + v_mid*n_y);
			ifv->F_u   = ifv->F_rho*u_mid + p_mid*n_x;
			ifv->F_v   = ifv->F_rho*v_mid + p_mid*n_y;
		}
	ifv->F_e   = (gamma_mid/(gamma_mid-1.0))*p_mid/rho_mid + 0.5*u_mid*u_mid;
	if (dim == 2)
	    ifv->F_e += 0.5*v_mid*v_mid;
	ifv->F_e = ifv->F_rho*ifv->F_e;

#ifdef MULTIFLUID_BASICS
	ifv->F_phi = ifv->F_rho * phi_mid;
	if (!isinf(config[60]))
		ifv->F_gamma = ifv->F_rho*gamma_mid;
	ifv->F_e_a  = z_a_mid/(config[6]-1.0)*p_mid/rho_mid + 0.5*phi_mid*u_mid*u_mid;
	if (dim == 2)
	    ifv->F_e_a += 0.5*phi_mid*v_mid*v_mid;
	ifv->F_e_a  = ifv->F_rho*ifv->F_e_a;
#endif
	
#ifdef MULTIFLUID_BASICS
	ifv->U_qt_add_c = ifv->F_rho*u_mid*phi_mid;
	if (dim == 2)
	    ifv->V_qt_add_c = ifv->F_rho*v_mid*phi_mid;
	ifv->U_qt_star  = p_mid*n_x;
	ifv->V_qt_star  = p_mid*n_y;
	ifv->P_star     = p_mid/rho_mid*ifv->F_rho;
#endif
	return 0;
}


/**
 * @brief This function calculate Eulerian fluxes of 2-D Euler equations by 2-D GRP solver.
 * @param[in,out] ifv: Structure pointer of interfacial evaluated variables and fluxes and left state.
 * @param[in] ifv_R:   Structure pointer of interfacial right state.
 * @param[in] tau:     The length of the time step.
 * @return    miscalculation indicator.
 *   @retval  0: Successful calculation.
 *   @retval  1: < 0.0 error.
 *   @retval  2: NAN or INFinite error of mid[].
 *   @retval  3: NAN or INFinite error of dire[].
 */
int GRP_2D_flux(struct i_f_var * ifv, struct i_f_var * ifv_R, const double tau)
{
	const double eps = config[4];
	const double n_x = ifv->n_x, n_y = ifv->n_y;
	double gamma_mid = ifv->gamma;
	ifv->lambda_u = 0.0;  ifv->lambda_v = 0.0;

	double u, u_R, d_u, d_u_R, t_u, t_u_R;
	u          =  ifv->U    *n_x + ifv->V    *n_y;
	u_R        =  ifv_R->U  *n_x + ifv_R->V  *n_y;
	d_u        =  ifv->d_u  *n_x + ifv->d_v  *n_y;
	d_u_R      =  ifv_R->d_u*n_x + ifv_R->d_v*n_y;
	t_u        =  ifv->t_u  *n_x + ifv->t_v  *n_y;
	t_u_R      =  ifv_R->t_u*n_x + ifv_R->t_v*n_y;
	ifv->V     = -ifv->U    *n_y + ifv->V    *n_x;
	ifv_R->V   = -ifv_R->U  *n_y + ifv_R->V  *n_x;
	ifv->d_v   = -ifv->d_u  *n_y + ifv->d_v  *n_x;
	ifv_R->d_v = -ifv_R->d_u*n_y + ifv_R->d_v*n_x;
	ifv->t_v   = -ifv->t_u  *n_y + ifv->t_v  *n_x;
	ifv_R->t_v = -ifv_R->t_u*n_y + ifv_R->t_v*n_x;
	ifv->U     =  u;
	ifv_R->U   =  u_R;
	ifv->d_u   =  d_u;
	ifv_R->d_u =  d_u_R;
	ifv->t_u   =  t_u;
	ifv_R->t_u =  t_u_R;
	
	double wave_speed[2], dire[6], mid[6], star[6];

	// linear_GRP_solver_Edir_G2D(wave_speed, dire, mid, star, ifv, ifv_R, eps, eps);
	// linear_GRP_solver_Edir_G2D(wave_speed, dire, mid, star, ifv, ifv_R, eps, INFINITY);
	linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, ifv, ifv_R, eps, eps);
	// linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, ifv, ifv_R, eps, INFINITY);

	if(mid[3] < eps || mid[0] < eps)
	    return 1;
	if(!isfinite(mid[1])|| !isfinite(mid[2])|| !isfinite(mid[0])|| !isfinite(mid[3]))
	    return 2;
	if(!isfinite(dire[1])|| !isfinite(dire[2])|| !isfinite(dire[0])|| !isfinite(dire[3]))
	    return 3;

	double rho_mid, p_mid, u_mid, v_mid;
	rho_mid =  mid[0] + 0.5*tau*dire[0];
	u_mid   = (mid[1] + 0.5*tau*dire[1])*n_x - (mid[2] + 0.5*tau*dire[2])*n_y;
	v_mid   = (mid[1] + 0.5*tau*dire[1])*n_y + (mid[2] + 0.5*tau*dire[2])*n_x;
	p_mid   =  mid[3] + 0.5*tau*dire[3];
#ifdef MULTIFLUID_BASICS
	double phi_mid, z_a_mid;
	phi_mid =  mid[5] + 0.5*tau*dire[5];
	z_a_mid =  mid[4] + 0.5*tau*dire[4];
	gamma_mid = 1.0/(z_a_mid/(config[6]-1.0)+(1.0-z_a_mid)/(config[106]-1.0))+1.0;
#endif

	ifv->F_rho = rho_mid*(u_mid*n_x + v_mid*n_y);
	ifv->F_u   = ifv->F_rho*u_mid + p_mid*n_x;
	ifv->F_v   = ifv->F_rho*v_mid + p_mid*n_y;
	ifv->F_e   = (gamma_mid/(gamma_mid-1.0))*p_mid/rho_mid + 0.5*(u_mid*u_mid + v_mid*v_mid);
	ifv->F_e   = ifv->F_rho*ifv->F_e;

	ifv->U_int   = (mid[1] + tau*dire[1])*n_x - (mid[2] + tau*dire[2])*n_y;
	ifv->V_int   = (mid[1] + tau*dire[1])*n_y + (mid[2] + tau*dire[2])*n_x;
	ifv->RHO_int =  mid[0] + tau*dire[0];
	ifv->P_int   =  mid[3] + tau*dire[3];

#ifdef MULTIFLUID_BASICS
	ifv->F_phi = ifv->F_rho*phi_mid;
	if (!isinf(config[60]))
		ifv->F_gamma = ifv->F_rho*gamma_mid;
	ifv->F_e_a = z_a_mid/(config[6]-1.0)*p_mid/rho_mid + 0.5*phi_mid*(u_mid*u_mid + v_mid*v_mid);
	ifv->F_e_a = ifv->F_rho*ifv->F_e_a;	
	ifv->PHI = mid[5] + tau*dire[5];
	ifv->Z_a = mid[4] + tau*dire[4];
#endif

#ifdef MULTIFLUID_BASICS
	ifv->U_qt_add_c = ifv->F_rho*u_mid*phi_mid;
	ifv->V_qt_add_c = ifv->F_rho*v_mid*phi_mid;
	ifv->U_qt_star  = p_mid*n_x;
	ifv->V_qt_star  = p_mid*n_y;
	ifv->P_star     = p_mid/rho_mid*ifv->F_rho;
#endif
	return 0;
}
