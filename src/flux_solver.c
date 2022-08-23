#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/Riemann_solver.h"
#include "../include/var_struc.h"

void Roe_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R)
{
	const int dim = (int)config[0];
	const double delta = 0.2;

	double F[4];
	double lambda_max;
	if (dim == 1)
		{
			Roe_solver(F, ifv->gamma, ifv->P, ifv->RHO, ifv->U, ifv_R->P, ifv_R->RHO, ifv_R->U, &lambda_max, delta);
			ifv->F_rho = F[0];
			ifv->F_u   = F[1];
			ifv->F_e   = F[2];
		}
	else if (dim == 2)
		{
			Roe_2D_solver(F, ifv->gamma, ifv->P, ifv->RHO, ifv->U, ifv->V, ifv->n_x, ifv->n_y, ifv_R->P, ifv_R->RHO, ifv_R->U, ifv_R->V, &lambda_max, delta);
			ifv->F_rho = F[0];
			ifv->F_u   = F[1];
			ifv->F_v   = F[2];
			ifv->F_e   = F[3];
		}
}

void HLL_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R)
{
	double F[4];
	double lambda_max;
	HLL_solver(F, ifv->gamma, ifv->P, ifv->RHO, ifv->U, ifv->V, ifv->n_x, ifv->n_y, ifv_R->P, ifv_R->RHO, ifv_R->U, ifv_R->V, &lambda_max);
	ifv->F_rho = F[0];
	ifv->F_u   = F[1];
	ifv->F_v   = F[2];
	ifv->F_e   = F[3];
}

void Riemann_exact_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R)
{
	const int dim = (int)config[0];
	const double eps = config[4];

	const double n_x = ifv->n_x, n_y = ifv->n_y;
	double u, u_R;
	if (dim == 1)
		{
			u   = ifv->U;
			u_R = ifv_R->U;
		}
	if (dim == 2)
		{
			u   = ifv->U  *n_x + ifv->V  *n_y;
			u_R = ifv_R->U*n_x + ifv_R->V*n_y;
		}

	double wave_speed[2], dire[6], mid[6], star[6];
/*	
	  linear_GRP_solver_Edir(wave_speed, dire, mid, star, 0.0, 0.000504, 1.4112, 0.0, 0.0, 1.24, 0.8787, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 1.0/1.4, 1.6272/1.4, 0.0, 0.0, 1.093, 1.4, eps);
	  printf("%.16lf,%.16lf,%.16lf,%.16lf,%.16lf,%.16lf\n",star[0],star[1],star[2],star[3],wave_speed[0],wave_speed[1]);
	  exit(1);
*/	
	linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, -0.0, -0.0, -0.0, -0.0, u, u_R, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, ifv->P, ifv_R->P, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, ifv->PHI, ifv_R->PHI, -0.0, -0.0, -0.0, -0.0, ifv->gamma, ifv_R->gamma, eps, -0.0);

	double rho_mid, p_mid, u_mid, v_mid, mid_qt;
	rho_mid = mid[0];
	p_mid = mid[3];
	double gamma_mid, z_a_mid;
	if(mid[1] > 0)
		{					
			gamma_mid = ifv->gamma;
			z_a_mid   = ifv->Z_a;
		}
	else
		{					
			gamma_mid = ifv_R->gamma;
			z_a_mid   = ifv_R->Z_a;
		}
	if (dim == 1)
		{
			u_mid = mid[1];

			ifv->F_rho = rho_mid*u_mid;
			ifv->F_u   = ifv->F_rho*u_mid + p_mid;
			ifv->F_e   = (gamma_mid/(gamma_mid-1.0))*p_mid/rho_mid + 0.5*u_mid*u_mid;
			ifv->F_e   = ifv->F_rho*ifv->F_e;
		}
	if (dim == 2)
		{
			if(mid[1] > 0)
				mid_qt = -ifv->U  *n_y + ifv->V  *n_x;
			else
				mid_qt = -ifv_R->U*n_y + ifv_R->V*n_x;
			u_mid = mid[1]*n_x - mid_qt*n_y;
			v_mid = mid[1]*n_y + mid_qt*n_x;

			ifv->F_rho = rho_mid*(u_mid*n_x + v_mid*n_y);
			ifv->F_u   = ifv->F_rho*u_mid + p_mid*n_x;
			ifv->F_v   = ifv->F_rho*v_mid + p_mid*n_y;
			ifv->F_e   = (gamma_mid/(gamma_mid-1.0))*p_mid/rho_mid + 0.5*(u_mid*u_mid + v_mid*v_mid);
			ifv->F_e   = ifv->F_rho*ifv->F_e;
		}
	double phi_mid;
	if ((int)config[2] == 2)
		{
			if(mid[1] > 0)
				phi_mid = ifv->PHI;
			else
				phi_mid = ifv_R->PHI;
			ifv->F_phi = ifv->F_rho * phi_mid;
		}
	if (!isinf(config[60]))
		ifv->F_gamma = ifv->F_rho*gamma_mid;
	
	ifv->F_e_a  = z_a_mid/(config[6]-1.0)*p_mid/rho_mid + 0.5*phi_mid*(u_mid*u_mid + v_mid*v_mid);
	ifv->F_e_a  = ifv->F_rho*ifv->F_e_a;
	ifv->P_star = p_mid/rho_mid*ifv->F_rho;
	ifv->U_qt_add_c = ifv->F_rho*u_mid*phi_mid;
	ifv->V_qt_add_c = ifv->F_rho*v_mid*phi_mid;
	ifv->U_qt_star  = p_mid*n_x;
	ifv->V_qt_star  = p_mid*n_y;
}


void GRP_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, const double tau)
{
	const int dim = (int)config[0];
	const double eps = config[4];
	const double n_x = ifv->n_x, n_y = ifv->n_y;

	double u, u_R, d_u, d_u_R, v, v_R, d_v, d_v_R;
	if (dim == 1)
		{
			u     = ifv->U;
			u_R   = ifv_R->U;
			d_u   = ifv->d_u;
			d_u_R = ifv_R->d_u;
		}
	else if (dim == 2)
		{
			u     =  ifv->U    *n_x + ifv->V    *n_y;
			u_R   =  ifv_R->U  *n_x + ifv_R->V  *n_y;
			d_u   =  ifv->d_u  *n_x + ifv->d_v  *n_y;
			d_u_R =  ifv_R->d_u*n_x + ifv_R->d_v*n_y;
			v     = -ifv->U    *n_y + ifv->V    *n_x;
			v_R   = -ifv_R->U  *n_y + ifv_R->V  *n_x;
			d_v   = -ifv->d_u  *n_y + ifv->d_v  *n_x;
			d_v_R = -ifv_R->d_u*n_y + ifv_R->d_v*n_x;
		}

	double wave_speed[2], dire[6], mid[6], star[6];
	double gamma = ifv->gamma;
	double gamma_mid,z_a_mid;

	double rho_mid, p_mid, u_mid, v_mid, phi_mid, mid_qt;
	if (dim == 1)
		{
			linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, -0.0, -0.0, u, u_R, d_u, d_u_R, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, ifv->PHI, ifv_R->PHI, ifv->d_phi, ifv_R->d_phi, -0.0, -0.0, ifv->gamma, ifv_R->gamma, eps, eps);

			rho_mid = mid[0] + 0.5*tau*dire[0];
			u_mid   = mid[1] + 0.5*tau*mid[0]*dire[1]/rho_mid;
			p_mid   = mid[3] + 0.5*tau*dire[3] + (gamma-1.0)*0.5*(tau*(dire[0]*0.5*mid[1]*mid[1]+mid[0]*mid[1]*dire[1]) + mid[0]*mid[1]*mid[1]-rho_mid*u_mid*u_mid);

			ifv->F_rho = rho_mid*u_mid;
			ifv->F_u   = ifv->F_rho*u_mid + p_mid;
			ifv->F_e   = (gamma/(gamma-1.0))*p_mid/rho_mid + 0.5*u_mid*u_mid;
			ifv->F_e   = ifv->F_rho*ifv->F_e;
			if ((int)config[2] == 2)
				{
					phi_mid = mid[5] + 0.5*tau*mid[0]*dire[5]/rho_mid;
					ifv->F_phi = ifv->F_rho*phi_mid;
				}
 -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0,
			ifv->RHO = mid[0] + tau*dire[0];
			ifv->U   = mid[1] + tau*mid[0]*dire[1]/ifv->RHO;
			ifv->P   = mid[3] + tau*dire[3] + (gamma-1.0)*0.5*(tau*(dire[0]*mid[1]*mid[1]+2*mid[0]*mid[1]*dire[1]) + mid[0]*mid[1]*mid[1]-ifv->RHO*ifv->U*ifv->U);
			if ((int)config[2] == 2)
				ifv->PHI = mid[5] + tau*mid[0]*dire[5]/ifv->RHO;
			ifv->gamma = gamma;
		}
	else if (dim == 2)
		{
			linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, -0.0, -0.0, u, u_R, d_u, d_u_R, -0.0, -0.0, v, v_R, d_v, d_v_R, -0.0, -0.0, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, -0.0, -0.0, ifv->Z_a, ifv_R->Z_a, ifv->d_z_a, ifv_R->d_z_a, -0.0, -0.0, ifv->PHI, ifv_R->PHI, ifv->d_phi, ifv_R->d_phi, -0.0, -0.0, ifv->gamma, ifv_R->gamma, eps, eps);

			z_a_mid = mid[4] + 0.5*tau*dire[4];
			gamma_mid = 1.0/(z_a_mid/(config[6]-1.0)+(1.0-z_a_mid)/(config[106]-1.0))+1.0;
			/*
			  rho_mid = mid[0] + 0.5*tau*dire[0];
			  u_mid   = (mid[1] + 0.5*tau*mid[0]*dire[1]/rho_mid)*n_x - (mid[2] + 0.5*tau*mid[0]*dire[2]/rho_mid)*n_y;
			  v_mid   = (mid[1] + 0.5*tau*mid[0]*dire[1]/rho_mid)*n_y + (mid[2] + 0.5*tau*mid[0]*dire[2]/rho_mid)*n_x;
			  p_mid   = mid[3] + 0.5*tau*dire[3] + (gamma_mid-1.0)*0.5*(tau*(dire[0]*0.5*(mid[1]*mid[1]+mid[2]*mid[2])+mid[0]*(mid[1]*dire[1]+mid[2]*dire[2])) + mid[0]*mid[1]*mid[1]+mid[0]*mid[2]*mid[2]-rho_mid*u_mid*u_mid-rho_mid*v_mid*v_mid);
			*/
			rho_mid = mid[0] + 0.5*tau*dire[0];
			u_mid   = (mid[1] + 0.5*tau*dire[1])*n_x - (mid[2] + 0.5*tau*dire[2])*n_y;
			v_mid   = (mid[1] + 0.5*tau*dire[1])*n_y + (mid[2] + 0.5*tau*dire[2])*n_x;
			p_mid   = mid[3] + 0.5*tau*dire[3];

			ifv->F_rho = rho_mid*(u_mid*n_x + v_mid*n_y);
			ifv->F_u   = ifv->F_rho*u_mid + p_mid*n_x;
			ifv->F_v   = ifv->F_rho*v_mid + p_mid*n_y;
			ifv->F_e   = (gamma_mid/(gamma_mid-1.0))*p_mid/rho_mid + 0.5*(u_mid*u_mid + v_mid*v_mid);
			ifv->F_e   = ifv->F_rho*ifv->F_e;

			ifv->U   = (mid[1] + tau*dire[1])*n_x - (mid[2] + tau*dire[2])*n_y;
			ifv->V   = (mid[1] + tau*dire[1])*n_y + (mid[2] + tau*dire[2])*n_x;
			if ((int)config[2] == 2)
				{
					//phi_mid = mid[5] + 0.5*tau*mid[0]*dire[5]/rho_mid;
					phi_mid = mid[5] + 0.5*tau*dire[5];
					ifv->F_phi = ifv->F_rho*phi_mid;
				}
			ifv->F_e_a = z_a_mid/(config[6]-1.0)*p_mid/rho_mid + 0.5*phi_mid*(u_mid*u_mid + v_mid*v_mid);
			ifv->F_e_a = ifv->F_rho*ifv->F_e_a;
			ifv->P_star = p_mid/rho_mid*ifv->F_rho;
			ifv->U_qt_add_c = ifv->F_rho*u_mid*phi_mid;
			ifv->V_qt_add_c = ifv->F_rho*v_mid*phi_mid;
			ifv->U_qt_star  = p_mid*n_x;
			ifv->V_qt_star  = p_mid*n_y;

			ifv->RHO = mid[0] + tau*dire[0];
			ifv->P   = mid[3] + tau*dire[3];
			ifv->Z_a = mid[4] + tau*dire[4];
			ifv->PHI = mid[5] + tau*dire[5];
		}
}


void GRP_2D_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, const double tau)
{
	const int dim = (int)config[0];
	const double eps = config[4];
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
	double gamma_mid,z_a_mid;
	double rho_mid, p_mid, u_mid, v_mid, phi_mid, mid_qt;

	linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, ifv->t_rho, ifv_R->t_rho, u, u_R, d_u, d_u_R, t_u, t_u_R, v, v_R, d_v, d_v_R, t_v, t_v_R, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, ifv->t_p, ifv_R->t_p, ifv->Z_a, ifv_R->Z_a, ifv->d_z_a, ifv_R->d_z_a, ifv->t_z_a, ifv_R->t_z_a, ifv->PHI, ifv_R->PHI, ifv->d_phi, ifv_R->d_phi, ifv->t_phi, ifv_R->t_phi, ifv->gamma, ifv_R->gamma, eps, -0.0);
//	linear_GRP_solver_Edir_G2D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, ifv->t_rho, ifv_R->t_rho, u, u_R, d_u, d_u_R, t_u, t_u_R, v, v_R, d_v, d_v_R, t_v, t_v_R, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, ifv->t_p, ifv_R->t_p, ifv->Z_a, ifv_R->Z_a, ifv->d_z_a, ifv_R->d_z_a, ifv->t_z_a, ifv_R->t_z_a, ifv->PHI, ifv_R->PHI, ifv->d_phi, ifv_R->d_phi, ifv->t_phi, ifv_R->t_phi, ifv->gamma, ifv_R->gamma, eps, eps);
//	linear_GRP_solver_Edir_G2D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, ifv->t_rho, ifv_R->t_rho, u, u_R, d_u, d_u_R, t_u, t_u_R, v, v_R, d_v, d_v_R, t_v, t_v_R, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, ifv->t_p, ifv_R->t_p, ifv->Z_a, ifv_R->Z_a, ifv->d_z_a, ifv_R->d_z_a, ifv->t_z_a, ifv_R->t_z_a, ifv->PHI, ifv_R->PHI, ifv->d_phi, ifv_R->d_phi, ifv->t_phi, ifv_R->t_phi, ifv->gamma, ifv_R->gamma, eps, -0.0);

// Acoustic approximation
//	linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, ifv->RHO, ifv_R->RHO, ifv->d_rho, ifv_R->d_rho, ifv->t_rho, ifv_R->t_rho, u, u_R, d_u, d_u_R, t_u, t_u_R, v, v_R, d_v, d_v_R, t_v, t_v_R, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, ifv->t_p, ifv_R->t_p, ifv->Z_a, ifv_R->Z_a, ifv->d_z_a, ifv_R->d_z_a, ifv->t_z_a, ifv_R->t_z_a, ifv->PHI, ifv_R->PHI, ifv->d_phi, ifv_R->d_phi, ifv->t_phi, ifv_R->t_phi, ifv->gamma, ifv_R->gamma, eps, 1.0/0.0);

	rho_mid = mid[0] + 0.5*tau*dire[0];
	u_mid   = (mid[1] + 0.5*tau*dire[1])*n_x - (mid[2] + 0.5*tau*dire[2])*n_y;
	v_mid   = (mid[1] + 0.5*tau*dire[1])*n_y + (mid[2] + 0.5*tau*dire[2])*n_x;
	p_mid   = mid[3] + 0.5*tau*dire[3];
	z_a_mid = mid[4] + 0.5*tau*dire[4];
	gamma_mid = 1.0/(z_a_mid/(config[6]-1.0)+(1.0-z_a_mid)/(config[106]-1.0))+1.0;
	phi_mid = mid[5] + 0.5*tau*dire[5];

	ifv->F_rho = rho_mid*(u_mid*n_x + v_mid*n_y);
	ifv->F_u   = ifv->F_rho*u_mid + p_mid*n_x;
	ifv->F_v   = ifv->F_rho*v_mid + p_mid*n_y;
	ifv->F_e   = (gamma_mid/(gamma_mid-1.0))*p_mid/rho_mid + 0.5*(u_mid*u_mid + v_mid*v_mid);
	ifv->F_e   = ifv->F_rho*ifv->F_e;
	ifv->F_phi = ifv->F_rho*phi_mid;
	
	ifv->U_int   = (mid[1] + tau*dire[1])*n_x - (mid[2] + tau*dire[2])*n_y;
	ifv->V_int   = (mid[1] + tau*dire[1])*n_y + (mid[2] + tau*dire[2])*n_x;

	ifv->F_e_a = z_a_mid/(config[6]-1.0)*p_mid/rho_mid + 0.5*phi_mid*(u_mid*u_mid + v_mid*v_mid);
	ifv->F_e_a = ifv->F_rho*ifv->F_e_a;
	ifv->P_star = p_mid/rho_mid*ifv->F_rho;
	ifv->U_qt_add_c = ifv->F_rho*u_mid*phi_mid;
	ifv->V_qt_add_c = ifv->F_rho*v_mid*phi_mid;
	ifv->U_qt_star  = p_mid*n_x;
	ifv->V_qt_star  = p_mid*n_y;

	ifv->RHO_int = mid[0] + tau*dire[0];
	ifv->P_int   = mid[3] + tau*dire[3];
	ifv->Z_a = mid[4] + tau*dire[4];
	ifv->PHI = mid[5] + tau*dire[5];
}
