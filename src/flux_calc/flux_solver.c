#include <stdio.h>
#include <math.h>

#include "../include/Riemann_solver.h"
#include "../include/var_struc.h"


void GRP_2D_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, const double tau)
{
	const double eps = config[4];
	double gamma_mid = config[6];
	const double n_x = ifv->n_x, n_y = ifv->n_y;

	double u, u_R, d_u, d_u_R, t_u, t_u_R, v, v_R, d_v, d_v_R, t_v, t_v_R;
	u     =  ifv->U    *n_x + ifv->V    *n_y;
	u_R   =  ifv_R->U  *n_x + ifv_R->V  *n_y;
	d_u   =  ifv->d_u  *n_x + ifv->d_v  *n_y;
	d_u_R =  ifv_R->d_u*n_x + ifv_R->d_v*n_y;
	t_u   =  ifv->t_u  *n_x + ifv->t_v  *n_y;
	t_u_R =  ifv_R->t_u*n_x + ifv_R->t_v*n_y;
	v     = -ifv->U    *n_y + ifv->V    *n_x;
	v_R   = -ifv_R->U  *n_y + ifv_R->V  *n_x;
	d_v   = -ifv->d_u  *n_y + ifv->d_v  *n_x;
	d_v_R = -ifv_R->d_u*n_y + ifv_R->d_v*n_x;
	t_v   = -ifv->t_u  *n_y + ifv->t_v  *n_x;
	t_v_R = -ifv_R->t_u*n_y + ifv_R->t_v*n_x;
	
	double wave_speed[2], dire[6], mid[6], star[6];
	double rho_mid, p_mid, u_mid, v_mid;

#ifdef MULTIFLUID_BASICS
	double phi_mid, z_a_mid;

	linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, ifv->t_rho, ifv_R->t_rho, u, u_R, d_u, d_u_R, t_u, t_u_R, v, v_R, d_v, d_v_R, t_v, t_v_R, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, ifv->t_p, ifv_R->t_p, ifv->Z_a, ifv_R->Z_a, ifv->d_z_a, ifv_R->d_z_a, ifv->t_z_a, ifv_R->t_z_a, ifv->PHI, ifv_R->PHI, ifv->d_phi, ifv_R->d_phi, ifv->t_phi, ifv_R->t_phi, ifv->gamma, ifv_R->gamma, eps, -0.0);
	// linear_GRP_solver_Edir_G2D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, ifv->t_rho, ifv_R->t_rho, u, u_R, d_u, d_u_R, t_u, t_u_R, v, v_R, d_v, d_v_R, t_v, t_v_R, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, ifv->t_p, ifv_R->t_p, ifv->Z_a, ifv_R->Z_a, ifv->d_z_a, ifv_R->d_z_a, ifv->t_z_a, ifv_R->t_z_a, ifv->PHI, ifv_R->PHI, ifv->d_phi, ifv_R->d_phi, ifv->t_phi, ifv_R->t_phi, ifv->gamma, ifv_R->gamma, eps, eps);
	// linear_GRP_solver_Edir_G2D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, ifv->t_rho, ifv_R->t_rho, u, u_R, d_u, d_u_R, t_u, t_u_R, v, v_R, d_v, d_v_R, t_v, t_v_R, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, ifv->t_p, ifv_R->t_p, ifv->Z_a, ifv_R->Z_a, ifv->d_z_a, ifv_R->d_z_a, ifv->t_z_a, ifv_R->t_z_a, ifv->PHI, ifv_R->PHI, ifv->d_phi, ifv_R->d_phi, ifv->t_phi, ifv_R->t_phi, ifv->gamma, ifv_R->gamma, eps, -0.0);
// Acoustic approximation
	// linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, ifv->t_rho, ifv_R->t_rho, u, u_R, d_u, d_u_R, t_u, t_u_R, v, v_R, d_v, d_v_R, t_v, t_v_R, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, ifv->t_p, ifv_R->t_p, ifv->Z_a, ifv_R->Z_a, ifv->d_z_a, ifv_R->d_z_a, ifv->t_z_a, ifv_R->t_z_a, ifv->PHI, ifv_R->PHI, ifv->d_phi, ifv_R->d_phi, ifv->t_phi, ifv_R->t_phi, ifv->gamma, ifv_R->gamma, eps, 1.0/0.0);
#else
	linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, ifv->t_rho, ifv_R->t_rho, u, u_R, d_u, d_u_R, t_u, t_u_R, v, v_R, d_v, d_v_R, t_v, t_v_R, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, ifv->t_p, ifv_R->t_p, 1.0, 1.0, -0.0, -0.0, -0.0, -0.0, 1.0, 1.0, -0.0, -0.0, -0.0, -0.0, gamma_mid, gamma_mid, eps, -0.0);
#endif

	rho_mid =  mid[0] + 0.5*tau*dire[0];
	u_mid   = (mid[1] + 0.5*tau*dire[1])*n_x - (mid[2] + 0.5*tau*dire[2])*n_y;
	v_mid   = (mid[1] + 0.5*tau*dire[1])*n_y + (mid[2] + 0.5*tau*dire[2])*n_x;
	p_mid   =  mid[3] + 0.5*tau*dire[3];

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
	phi_mid = mid[5] + 0.5*tau*dire[5];
	z_a_mid = mid[4] + 0.5*tau*dire[4];
	gamma_mid = 1.0/(z_a_mid/(config[6]-1.0)+(1.0-z_a_mid)/(config[106]-1.0))+1.0;
	ifv->F_phi = ifv->F_rho*phi_mid;
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
}
