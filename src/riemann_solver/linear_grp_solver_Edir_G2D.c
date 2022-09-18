/**
 * @file  linear_GRP_solver_Edir_G2D.c
 * @brief This is a Genuinely-2D direct Eulerian GRP solver for compressible inviscid flow in Li's paper.
 */

#include <math.h>
#include <stdio.h>

#include "../include/var_struc.h"
#include "../include/riemann_solver.h"

/**
 * @def EXACT_TANGENT_DERIVATIVE
 * @brief Switch whether the tangential derivatives are accurately computed.
 */
#ifdef DOXYGEN_PREDEFINED
#define EXACT_TANGENT_DERIVATIVE
#endif


/**
 * @brief A Genuinely-2D direct Eulerian GRP solver for unsteady compressible inviscid two-component flow in two space dimension.
 * @param[out] wave_speed: the velocity of left and right waves.
 * @param[out] D: the temporal derivative of fluid variables. \n
 *                   [rho, u, v, p, phi, z_a]_t
 * @param[out] U:  the intermediate Riemann solutions at t-axis. \n
 *                   [rho_mid, u_mid, v_mid, p_mid, phi_mid, z_a_mid]
 * @param[out] U_star: the Riemann solutions in star region. \n
 *                   [rho_star_L, u_star, rho_star_R, p_star, c_star_L, c_star_R]
 * @param[in] ifv_L: Left  States (rho/u/v/p/phi/z, d_, t_, gammaL).
 * @param[in] ifv_R: Right States (rho/u/v/p/phi/z, d_, t_, gammaR).
 *                   - s_: normal derivatives.
 *                   - t_: tangential derivatives.
 *                   - gamma: the constant of the perfect gas.
 * @param[in] eps: the largest value could be seen as zero.
 * @param[in] atc: Parameter that determines the solver type.
 *              - INFINITY: acoustic approximation
 *                - ifv_.s_, ifv_.t_ = -0.0: exact Riemann solver 
 *              - eps:      Genuinely-2D GRP solver(nonlinear + acoustic case)
 *                - ifv_.t_ = -0.0: Planar-1D GRP solver
 *              - -0.0:     Genuinely-2D GRP solver(only nonlinear case)
 *                - ifv_.t_ = -0.0: Planar-1D GRP solver
 * @remark macro definition \b EXACT_TANGENT_DERIVATIVE: \n
 *         Switch whether the tangential derivatives are accurately computed.
 * @sa   Theory is found in Reference [1]. \n
 *       [1] 齐进, 二维欧拉方程广义黎曼问题数值建模及其应用, Ph.D Thesis, Beijing Normal University, 2017.
 */
void linear_GRP_solver_Edir_G2D
(double *wave_speed, double *D, double *U, double *U_star, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R, const double eps, const double atc)
{
	const double lambda_u = ifv_L->lambda_u, lambda_v = ifv_R->lambda_v;
	const double  gammaL = ifv_L->gamma,  gammaR = ifv_R->gamma;
	const double   rho_L = ifv_L->RHO,     rho_R = ifv_R->RHO;
	const double d_rho_L = ifv_L->d_rho, d_rho_R = ifv_R->d_rho;
	const double t_rho_L = ifv_L->t_rho, t_rho_R = ifv_R->t_rho;
	const double     u_L = ifv_L->U,         u_R = ifv_R->U;
	const double   d_u_L = ifv_L->d_u,     d_u_R = ifv_R->d_u;
	const double   t_u_L = ifv_L->t_u,     t_u_R = ifv_R->t_u;
	const double     v_L = ifv_L->V,         v_R = ifv_R->V;
	const double   d_v_L = ifv_L->d_v,     d_v_R = ifv_R->d_v;
	const double   t_v_L = ifv_L->t_v,     t_v_R = ifv_R->t_v;
	const double     p_L = ifv_L->P,         p_R = ifv_R->P;
	const double   d_p_L = ifv_L->d_p,     d_p_R = ifv_R->d_p;
	const double   t_p_L = ifv_L->t_p,     t_p_R = ifv_R->t_p;
#ifdef MULTIFLUID_BASICS
	const double     z_L = ifv_L->Z_a,       z_R = ifv_R->Z_a;
	const double   d_z_L = ifv_L->d_z_a,   d_z_R = ifv_R->d_z_a;
	const double   t_z_L = ifv_L->t_z_a,   t_z_R = ifv_R->t_z_a;
	const double   phi_L = ifv_L->PHI,     phi_R = ifv_R->PHI;
	const double d_phi_L = ifv_L->d_phi, d_phi_R = ifv_R->d_phi;
	const double t_phi_L = ifv_L->t_phi, t_phi_R = ifv_R->t_phi;
#else
	const double     z_L =  0.0,     z_R =  0.0;
	const double   d_z_L = -0.0,   d_z_R = -0.0;
	const double   t_z_L = -0.0,   t_z_R = -0.0;
	const double   phi_L =  0.0,   phi_R =  0.0;
	const double d_phi_L = -0.0, d_phi_R = -0.0;
	const double t_phi_L = -0.0, t_phi_R = -0.0;
#endif

	_Bool CRW[2];
	double dist;
	double c_L, c_R, C, c_frac = 1.0;

	double d_Phi, d_Psi, TdS, VAR;
	double D_rho, D_u, D_v, D_p, D_z, D_phi, T_rho, T_u, T_v, T_p, T_z, T_phi; 
	double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R;
	double Q;

	double H1, H2, H3;
	double a_L, b_L, d_L, a_R, b_R, d_R, detA;
	double L_u, L_p, L_rho, L_v;

	double u_t_mat, p_t_mat, D0_p_tau, D0_u_tau;
	double SmUs, SmUL, SmUR;
  
	const double zetaL = (gammaL-1.0)/(gammaL+1.0);
	const double zetaR = (gammaR-1.0)/(gammaR+1.0);
 
	double rho_x, f;
	double speed_L, speed_R;
#ifdef EXACT_TANGENT_DERIVATIVE
	double da_y = 0.05*config[11];
	double gammaL_up, gammaR_up, gammaL_dn, gammaR_dn;
	double mid_up[6], star_up[6], mid_dn[6], star_dn[6];
	double wave_speed_tmp[2], dire_tmp[6];
#endif

	c_L = sqrt(gammaL * p_L / rho_L);
	c_R = sqrt(gammaR * p_R / rho_R);

	dist = sqrt((rho_L-rho_R)*(rho_L-rho_R) + (u_L-u_R)*(u_L-u_R) + (v_L-v_R)*(v_L-v_R) + (p_L-p_R)*(p_L-p_R));
	if (dist < atc && atc < 2*eps)
	    {
		u_star = 0.5*(u_R+u_L);
		p_star = 0.5*(p_R+p_L);
		rho_star_L = rho_L;
		c_star_L = c_L;
		speed_L = u_star - c_star_L;
		rho_star_R = rho_R;
		c_star_R = c_R;
		speed_R = u_star + c_star_R;
	    }
	else //=========Riemann solver==========
	    {
		Riemann_solver_exact(&u_star, &p_star, gammaL, gammaR, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, eps, 500);
		if(CRW[0])
		    {
			rho_star_L = rho_L*pow(p_star/p_L, 1.0/gammaL);
			c_star_L = c_L*pow(p_star/p_L, 0.5*(gammaL-1.0)/gammaL);
			speed_L = u_L - c_L;
		    }
		else
		    {
			rho_star_L = rho_L*(p_star+zetaL*p_L)/(p_L+zetaL*p_star);
			c_star_L = sqrt(gammaL * p_star / rho_star_L);
			speed_L = u_L - c_L*sqrt(0.5*((gammaL+1.0)*(p_star/p_L) + (gammaL-1.0))/gammaL);
		    }
		if(CRW[1])
		    {
			rho_star_R = rho_R*pow(p_star/p_R,1.0/gammaR);
			c_star_R = c_R*pow(p_star/p_R, 0.5*(gammaR-1.0)/gammaR);
			speed_R = u_R + c_R;
		    }
		else
		    {
			rho_star_R = rho_R*(p_star+zetaR*p_R)/(p_R+zetaR*p_star);
			c_star_R = sqrt(gammaR * p_star / rho_star_R);
			speed_R = u_R + c_R*sqrt(0.5*((gammaR+1.0)*(p_star/p_R) + (gammaR-1.0))/gammaR);
		    }
	    }
	wave_speed[0] = speed_L;
	wave_speed[1] = speed_R;

	//=========acoustic case==========
	if(dist < atc)
		{
			if(speed_L > lambda_u) //the direction is on the left side of all the three waves
				{
					U[0] = rho_L;
					U[1] =   u_L;
					U[2] =   v_L;
					U[3] =   p_L;
					U[4] =   z_L;
					U[5] = phi_L;
					D[0] = -(u_L-lambda_u)*d_rho_L - (v_L-lambda_v)*t_rho_L - rho_L*(d_u_L+t_v_L);
					D[1] = -(u_L-lambda_u)*d_u_L   - (v_L-lambda_v)*t_u_L   - d_p_L/rho_L;
					D[2] = -(u_L-lambda_u)*d_v_L   - (v_L-lambda_v)*t_v_L   - t_p_L/rho_L;
					D[3] = -(u_L-lambda_u)*d_p_L   - (v_L-lambda_v)*t_p_L   - rho_L*c_L*c_L*(d_u_L+t_v_L) ;
					D[4] = -(u_L-lambda_u)*d_z_L   - (v_L-lambda_v)*t_z_L;
					D[5] = -(u_L-lambda_u)*d_phi_L - (v_L-lambda_v)*t_phi_L;
				}
			else if(speed_R < lambda_u) //the direction is on the right side of all the three waves
				{
					U[0] = rho_R;
					U[1] =   u_R;
					U[2] =   v_R;
					U[3] =   p_R;
					U[4] =   z_R;
					U[5] = phi_R;
					D[0] = -(u_R-lambda_u)*d_rho_R - (v_R-lambda_v)*t_rho_R - rho_R*(d_u_R+t_v_R);
					D[1] = -(u_R-lambda_u)*d_u_R   - (v_R-lambda_v)*t_u_R   - d_p_R/rho_R;
					D[2] = -(u_R-lambda_u)*d_v_R   - (v_R-lambda_v)*t_v_R   - t_p_R/rho_R;
					D[3] = -(u_R-lambda_u)*d_p_R   - (v_R-lambda_v)*t_p_R   - rho_R*c_R*c_R*(d_u_R+t_v_R);
					D[4] = -(u_R-lambda_u)*d_z_R   - (v_R-lambda_v)*t_z_R;
					D[5] = -(u_R-lambda_u)*d_phi_R - (v_R-lambda_v)*t_phi_R;
				}
			else
				{
					if(CRW[0] && ((u_star-c_star_L) > lambda_u)) // the direction is in a 1-CRW
						{
							U[1] = zetaL*(u_L+2.0*(c_L+lambda_u)/(gammaL-1.0));
							C = U[1] - lambda_u;
							U[3] = pow(C/c_L, 2.0*gammaL/(gammaL-1.0)) * p_L;
							U[0] = gammaL*U[3]/C/C;
							U[2] = v_L;
							U[4] = z_L;
							U[5] = phi_L;
						}
					else if(CRW[1] && ((u_star+c_star_R) < lambda_u)) // the direction is in a 3-CRW
						{
							U[1] = zetaR*(u_R-2.0*(c_R-lambda_u)/(gammaR-1.0));
							C = lambda_u-U[1];
							U[3] = pow(C/c_R, 2.0*gammaR/(gammaR-1.0)) * p_R;
							U[0] = gammaR*U[3]/C/C;
							U[2] = v_R;
							U[4] = z_R;
							U[5] = phi_R;
						}	
					else if(u_star > lambda_u) //the direction is between the 1-wave and the contact discontinuety
						{
							U[0] = rho_star_L;
							U[1] =   u_star;
							U[2] =        v_L;
							U[3] =   p_star;
							U[4] =        z_L;
							U[5] =      phi_L;
							C    =   c_star_L;
						}
					else //the direction is between the contact discontinuety and the 3-wave
						{
							U[0] = rho_star_R;
							U[1] =   u_star;
							U[2] =        v_R;
							U[3] =   p_star;
							U[4] =        z_R;
							U[5] =      phi_R;
							C    =   c_star_R;
						}			

					D_p = 0.5*((d_u_L*(U[0]*C) + d_p_L) - (d_u_R*(U[0]*C) - d_p_R));
					T_p = 0.5*((t_u_L*(U[0]*C) + t_p_L) - (t_u_R*(U[0]*C) - t_p_R));
					D_u = 0.5*(d_u_L + d_p_L/(U[0]*C) + d_u_R - d_p_R/(U[0]*C));
					T_u = 0.5*(t_u_L + t_p_L/(U[0]*C) + t_u_R - t_p_R/(U[0]*C));
					if(u_star > lambda_u)
						{
							D_v = d_v_L;
							T_v = t_v_L;
							D_z = d_z_L;
							T_z = t_z_L;
							D_phi = d_phi_L;
							T_phi = t_phi_L;
							D_rho = d_rho_L - d_p_L/(C*C) + D_p/(C*C);
							T_rho = t_rho_L - t_p_L/(C*C) + T_p/(C*C);
						}
					else
						{
							D_v = d_v_R;
							T_v = t_v_R;
							D_z = d_z_R;
							T_z = t_z_R;
							D_phi = d_phi_R;
							T_phi = t_phi_R;
							D_rho = d_rho_R - d_p_R/(C*C) + D_p/(C*C);
							T_rho = t_rho_R - t_p_R/(C*C) + T_p/(C*C);
						}
					D[0] = -(U[1]-lambda_u)*D_rho - (U[2]-lambda_v)*T_rho - U[0]*(D_u+T_v);
					D[1] = -(U[1]-lambda_u)*D_u   - (U[2]-lambda_v)*T_u   - D_p/U[0];
					D[2] = -(U[1]-lambda_u)*D_v   - (U[2]-lambda_v)*T_v   - T_p/U[0];
					D[3] = -(U[1]-lambda_u)*D_p   - (U[2]-lambda_v)*T_p   - U[0]*C*C*(D_u+T_v);
					D[4] = -(U[1]-lambda_u)*D_z   - (U[2]-lambda_v)*T_z;
					D[5] = -(U[1]-lambda_u)*D_phi - (U[2]-lambda_v)*T_phi;
				}				
			U_star[0] = rho_star_L;
			U_star[1] = u_star;
			U_star[2] = rho_star_R;
			U_star[3] = p_star;
			U_star[4] = c_star_L;
			U_star[5] = c_star_R;
			return;
		}

	//=========non-acoustic case==========
	//------trivial case------
	if(speed_L > lambda_u) //the direction is on the left side of all the three waves
		{
			U[0] = rho_L;
			U[1] =   u_L;
			U[2] =   v_L;
			U[3] =   p_L;
			U[4] =   z_L;
			U[5] = phi_L;
			D[0] = -(u_L-lambda_u)*d_rho_L - (v_L-lambda_v)*t_rho_L - rho_L*(d_u_L+t_v_L);
			D[1] = -(u_L-lambda_u)*d_u_L   - (v_L-lambda_v)*t_u_L   - d_p_L/rho_L;
			D[2] = -(u_L-lambda_u)*d_v_L   - (v_L-lambda_v)*t_v_L   - t_p_L/rho_L;
			D[3] = -(u_L-lambda_u)*d_p_L   - (v_L-lambda_v)*t_p_L   - rho_L*c_L*c_L*(d_u_L+t_v_L) ;
			D[4] = -(u_L-lambda_u)*d_z_L   - (v_L-lambda_v)*t_z_L;
			D[5] = -(u_L-lambda_u)*d_phi_L - (v_L-lambda_v)*t_phi_L;
		}
	else if(speed_R < lambda_u) //the direction is on the right side of all the three waves
		{
			U[0] = rho_R;
			U[1] =   u_R;
			U[2] =   v_R;
			U[3] =   p_R;
			U[4] =   z_R;
			U[5] = phi_R;
			D[0] = -(u_R-lambda_u)*d_rho_R - (v_R-lambda_v)*t_rho_R - rho_R*(d_u_R+t_v_R);
			D[1] = -(u_R-lambda_u)*d_u_R   - (v_R-lambda_v)*t_u_R   - d_p_R/rho_R;
			D[2] = -(u_R-lambda_u)*d_v_R   - (v_R-lambda_v)*t_v_R   - t_p_R/rho_R;
			D[3] = -(u_R-lambda_u)*d_p_R   - (v_R-lambda_v)*t_p_R   - rho_R*c_R*c_R*(d_u_R+t_v_R);
			D[4] = -(u_R-lambda_u)*d_z_R   - (v_R-lambda_v)*t_z_R;
			D[5] = -(u_R-lambda_u)*d_phi_R - (v_R-lambda_v)*t_phi_R;
		}
	else//----non-trivial case----
		{
			// calculate T_rho, T_u, T_v, T_p, T_z, T_phi
#ifdef EXACT_TANGENT_DERIVATIVE
			gammaL_up = 1.0/((z_L+da_y*t_z_L)/(config[6]-1.0)+(1.0-(z_L+da_y*t_z_L))/(config[106]-1.0))+1.0;
			gammaR_up = 1.0/((z_R+da_y*t_z_R)/(config[6]-1.0)+(1.0-(z_R+da_y*t_z_R))/(config[106]-1.0))+1.0;
			linear_GRP_solver_Edir_Q1D(wave_speed_tmp, dire_tmp, mid_up, star_up, 0.0, 0.0, rho_L+da_y*t_rho_L, rho_R+da_y*t_rho_R, -0.0, -0.0, -0.0, -0.0, u_L+da_y*t_u_L, u_R+da_y*t_u_R, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, p_L+da_y*t_p_L, p_R+da_y*t_p_R, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gammaL_up, gammaR_up, eps*da_y, -0.0);
			gammaL_dn = 1.0/((z_L-da_y*t_z_L)/(config[6]-1.0)+(1.0-(z_L-da_y*t_z_L))/(config[106]-1.0))+1.0;
			gammaR_dn = 1.0/((z_R-da_y*t_z_R)/(config[6]-1.0)+(1.0-(z_R-da_y*t_z_R))/(config[106]-1.0))+1.0;
			linear_GRP_solver_Edir_Q1D(wave_speed_tmp, dire_tmp, mid_dn, star_dn, 0.0, 0.0, rho_L-da_y*t_rho_L, rho_R-da_y*t_rho_R, -0.0, -0.0, -0.0, -0.0, u_L-da_y*t_u_L, u_R-da_y*t_u_R, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, p_L-da_y*t_p_L, p_R-da_y*t_p_R, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gammaL_dn, gammaR_dn, eps*da_y, -0.0);

			if (CRW[0] && ((u_star-c_star_L) > lambda_u||(star_up[1]-star_up[4]) > lambda_u||(star_dn[1]-star_dn[4]) > lambda_u)) //the direction is in a 1-CRW
				{
					T_u = (mid_up[1]-mid_dn[1])/da_y*0.5;
					T_p = (mid_up[3]-mid_dn[3])/da_y*0.5;					
					T_rho = (mid_up[0]-mid_dn[0])/da_y*0.5;					
				}
			else if (CRW[1] && ((u_star+c_star_R) < lambda_u||(star_up[1]+star_up[5]) < lambda_u||(star_dn[1]+star_dn[5]) < lambda_u)) //the direction is in a 3-CRW
				{
					T_u = (mid_up[1]-mid_dn[1])/da_y*0.5;
					T_p = (mid_up[3]-mid_dn[3])/da_y*0.5;					
					T_rho = (mid_up[0]-mid_dn[0])/da_y*0.5;					
				}
			else
				{
					T_u = (star_up[1]-star_dn[1])/da_y*0.5;
					T_p = (star_up[3]-star_dn[3])/da_y*0.5;					
					if(u_star < lambda_u)
						T_rho = (star_up[2]-star_dn[2])/da_y*0.5;	
					else
						T_rho = (star_up[0]-star_dn[0])/da_y*0.5;	
				}
#else
			if(u_star < lambda_u)
				{		   
					T_p = 0.5*((t_u_L-t_u_R)*rho_star_R*c_star_R+t_p_L+t_p_R);
					T_u = 0.5*(t_u_L+t_u_R+(t_p_L-t_p_R)/rho_star_R/c_star_R);
					T_rho = t_rho_R - t_p_R/(c_star_R*c_star_R) + T_p/(c_star_R*c_star_R);
				}
			else
				{	
					T_p = 0.5*((t_u_L-t_u_R)*rho_star_L*c_star_L+t_p_L+t_p_R);
					T_u = 0.5*(t_u_L+t_u_R+(t_p_L-t_p_R)/rho_star_L/c_star_L);
					T_rho = t_rho_L - t_p_L/(c_star_L*c_star_L) + T_p/(c_star_L*c_star_L);
				}
#endif
			if(CRW[0] && ((u_star-c_star_L) > lambda_u)) // the direction is in a 1-CRW
				{
					U[1] = zetaL*(u_L+2.0*(c_L+lambda_u)/(gammaL-1.0));
					C = U[1] - lambda_u;
					U[3] = pow(C/c_L, 2.0*gammaL/(gammaL-1.0)) * p_L;
					U[0] = gammaL*U[3]/C/C;
					U[2] = v_L;
					U[4] = z_L;
					U[5] = phi_L;

					c_frac = C/c_L;
					TdS = (d_p_L - d_rho_L*c_L*c_L)/(gammaL-1.0)/rho_L;
					d_Psi = d_u_L + (gammaL*d_p_L/c_L - c_L*d_rho_L)/(gammaL-1.0)/rho_L;
					D[1] = ((1.0+zetaL)*pow(c_frac, 0.5/zetaL) + zetaL*pow(c_frac, (1.0+zetaL)/zetaL));
					D[1] = D[1]/(1.0+2.0*zetaL) * TdS;
					D[1] = D[1] - c_L*pow(c_frac, 0.5/zetaL) * d_Psi;
					if (gammaL<3.0-eps || gammaL>3.0+eps)
						Q = (c_frac*(zetaL-1.0)+pow(c_frac, 0.5/zetaL)*zetaL)/(2.0*zetaL-1.0);
					else
						Q = 0.5*c_frac+pow(c_frac, 0.5/zetaL)*(0.5-0.25/zetaL*log(c_frac));
					D[1] = D[1] - c_L*t_v_L*Q;
					D[3] = U[0]*(U[1] - lambda_u)*D[1];

					D[0] = U[0]*(U[1] - lambda_u)*pow(c_frac, (1.0+zetaL)/zetaL)*TdS*(gammaL-1.0);
					D[0] = (D[0] + D[3]) / C/C - (U[2]-lambda_v)*T_rho;

					D[2] = -(U[1] - lambda_u)*d_v_L*U[0]/rho_L - (U[2]-lambda_v)*t_v_L - T_p/U[0];
					D[2] = D[2] - (zetaL-1.0)*(pow(c_frac, 2.0/zetaL-1.0)-1.0)/(zetaL-2.0)/U[0] * t_p_L;
					D[4] = -(U[1] - lambda_u)*d_z_L*U[0]/rho_L - (U[2]-lambda_v)*t_z_L;
					D[5] = -(U[1] - lambda_u)*d_phi_L*U[0]/rho_L - (U[2]-lambda_v)*t_phi_L;

					D[3] = D[3] - (U[2]-lambda_v)*T_p;
					D[1] = D[1] - (U[2]-lambda_v)*T_u + U[1]*t_v_L;					
				}
			else if(CRW[1] && ((u_star+c_star_R) < lambda_u)) // the direction is in a 3-CRW
				{
					U[1] = zetaR*(u_R-2.0*(c_R-lambda_u)/(gammaR-1.0));
					C = lambda_u-U[1];
					U[3] = pow(C/c_R, 2.0*gammaR/(gammaR-1.0)) * p_R;
					U[0] = gammaR*U[3]/C/C;
					U[2] = v_R;
					U[4] = z_R;
					U[5] = phi_R;

					c_frac = C/c_R;
					TdS = (d_p_R - d_rho_R*c_R*c_R)/(gammaR-1.0)/rho_R;
					d_Phi = d_u_R - (gammaR*d_p_R/c_R - c_R*d_rho_R)/(gammaR-1.0)/rho_R;
					D[1] = ((1.0+zetaR)*pow(c_frac, 0.5/zetaR) + zetaR*pow(c_frac, (1.0+zetaR)/zetaR));
					D[1] = D[1]/(1.0+2.0*zetaR) * TdS;
					D[1] = D[1] + c_R*pow(c_frac, 0.5/zetaR)*d_Phi;
					if (gammaR<3.0-eps || gammaR>3.0+eps)
						Q = (c_frac*(zetaR-1.0)+pow(c_frac, 0.5/zetaR)*zetaR)/(2.0*zetaR-1.0);
					else
						Q = 0.5*c_frac+pow(c_frac, 0.5/zetaR)*(0.5-0.25/zetaR*log(c_frac));
					D[1] = D[1] + c_R*t_v_R*Q;
					D[3] = U[0]*(U[1]-lambda_u)*D[1];

					D[0] = U[0]*(U[1]-lambda_u)*pow(c_frac, (1.0+zetaR)/zetaR)*TdS*(gammaR-1.0);
					D[0] = (D[0] + D[3]) / C/C - (U[2]-lambda_v)*T_rho;

					D[2] = -(U[1]-lambda_u)*d_v_R*U[0]/rho_R - (U[2]-lambda_v)*t_v_R - T_p/U[0];
					D[2] = D[2] - (zetaR-1.0)*(pow(c_frac, 2.0/zetaR-1.0)-1.0)/(zetaR-2.0)/U[0] * t_p_R;
					D[4] = -(U[1]-lambda_u)*d_z_R*U[0]/rho_R - (U[2]-lambda_v)*t_z_R;
					D[5] = -(U[1]-lambda_u)*d_phi_R*U[0]/rho_R - (U[2]-lambda_v)*t_phi_R;

					D[3] = D[3] - (U[2]-lambda_v)*T_p;
					D[1] = D[1] - (U[2]-lambda_v)*T_u + U[1]*t_v_R;
				}
			else//--non-sonic case--
				{
					if(u_star < lambda_u) //the direction is between the contact discontinuety and the 3-wave
						{
							U[0] = rho_star_R;
							U[1] =   u_star;
							U[2] =   v_R;
							U[3] =   p_star;
							U[4] =   z_R;
							U[5] = phi_R;
							C = c_star_R;
							T_v = t_v_R;
						}
					else //the direction is between the 1-wave and the contact discontinuety
						{
							U[0] = rho_star_L;
							U[1] =   u_star;
							U[2] =   v_L;
							U[3] =   p_star;
							U[4] =   z_L;
							U[5] = phi_L;
							C = c_star_L;
							T_v = t_v_L;
						}

					//determine a_L, b_L and d_L
					if(CRW[0]) //the 1-wave is a CRW
						{
							a_L = 1.0;
							b_L = 1.0 / rho_star_L / c_star_L;
							c_frac = c_star_L/c_L;
							TdS = (d_p_L - d_rho_L*c_L*c_L)/(gammaL-1.0)/rho_L;
							d_Psi = d_u_L + (gammaL*d_p_L/c_L - c_L*d_rho_L)/(gammaL-1.0)/rho_L;
							d_L = ((1.0+zetaL)*pow(c_frac, 0.5/zetaL) + zetaL*pow(c_frac, (1.0+zetaL)/zetaL));
							d_L = d_L/(1.0+2.0*zetaL) * TdS;
							d_L = d_L - c_L*pow(c_frac, 0.5/zetaL) * d_Psi;
							if (gammaL<3.0-eps || gammaL>3.0+eps)
								Q = (c_frac*(zetaL-1.0)+pow(c_frac, 0.5/zetaL)*zetaL)/(2.0*zetaL-1.0);
							else
								Q = 0.5*c_frac+pow(c_frac, 0.5/zetaL)*(0.5-0.25/zetaL*log(c_frac));
							d_L = d_L - c_L*t_v_L*Q;
						}
					else //the 1-wave is a shock
						{
							SmUs = -sqrt(0.5*((gammaL+1.0)*p_L   +(gammaL-1.0)*p_star)/rho_star_L);
							SmUL = -sqrt(0.5*((gammaL+1.0)*p_star+(gammaL-1.0)*p_L   )/rho_L);

							VAR = sqrt((1-zetaL)/(rho_L*(p_star+zetaL*p_L)));

							H1 =  0.5*VAR * (p_star+(1.0+2.0*zetaL)*p_L)/(p_star+zetaL*p_L);
							H2 = -0.5*VAR * ((2.0+zetaL)*p_star + zetaL*p_L)/(p_star+zetaL*p_L);
							H3 = -0.5*VAR * (p_star-p_L) / rho_L;

							L_p = -1.0/rho_L - SmUL*H2;
							L_u = SmUL + rho_L*(c_L*c_L*H2 + H3);
							L_rho = -SmUL * H3;
							L_v = SmUs + rho_L*(c_L*c_L*H2 + H3);

							a_L = 1.0 - rho_star_L* SmUs * H1;
							b_L = -SmUs/(rho_star_L*c_star_L*c_star_L)+ H1;
							d_L = L_rho*d_rho_L + L_u*d_u_L + L_p*d_p_L + L_v*t_v_L;
						}
					d_L = d_L - a_L*v_L*T_u - b_L*v_L*T_p;
					//determine a_R, b_R and d_R
					if(CRW[1]) //the 3-wave is a CRW
						{
							a_R = 1.0;
							b_R = -1.0 / rho_star_R / c_star_R;
							c_frac = c_star_R/c_R;
							TdS = (d_p_R - d_rho_R*c_R*c_R)/(gammaR-1.0)/rho_R;
							d_Phi = d_u_R - (gammaR*d_p_R/c_R - c_R*d_rho_R)/(gammaR-1.0)/rho_R;
							d_R = ((1.0+zetaR)*pow(c_frac, 0.5/zetaR) + zetaR*pow(c_frac, (1.0+zetaR)/zetaR));
							d_R = d_R/(1.0+2.0*zetaR) * TdS;
							d_R = d_R + c_R*pow(c_frac, 0.5/zetaR) * d_Phi;
							if (gammaR<3.0-eps || gammaR>3.0+eps)
								Q = (c_frac*(zetaR-1.0)+pow(c_frac, 0.5/zetaR)*zetaR)/(2.0*zetaR-1.0);
							else
								Q = 0.5*c_frac+pow(c_frac, 0.5/zetaR)*(0.5-0.25/zetaR*log(c_frac));
							d_R = d_R + c_R*t_v_R*Q;
						}
					else //the 3-wave is a shock
						{
							SmUs = sqrt(0.5*((gammaR+1.0)*p_R   + (gammaR-1.0)*p_star)/rho_star_R);
							SmUR = sqrt(0.5*((gammaR+1.0)*p_star+ (gammaR-1.0)*p_R   )/rho_R);

							VAR  = sqrt((1.0-zetaR)/(rho_R*(p_star+zetaR*p_R)));

							H1 =  0.5*VAR * (p_star+(1+2.0*zetaR)*p_R)/(p_star+zetaR*p_R);
							H2 = -0.5*VAR * ((2.0+zetaR)*p_star+zetaR*p_R)/(p_star+zetaR*p_R);
							H3 = -0.5*VAR * (p_star-p_R) /rho_R;

							L_p = -1.0/rho_R + SmUR*H2;
							L_u = SmUR - rho_R*(c_R*c_R*H2 + H3);
							L_rho = SmUR * H3;
							L_v = SmUs - rho_R*(c_R*c_R*H2 + H3);

							a_R = 1.0 +rho_star_R* SmUs * H1;
							b_R = -(SmUs/(rho_star_R*c_star_R*c_star_R) + H1);
							d_R = L_rho*d_rho_R + L_u*d_u_R + L_p*d_p_R + L_v*t_v_R;
						}
					d_R = d_R - a_R*v_R*T_u - b_R*v_R*T_p;

					detA = a_L*b_R - b_L*a_R;
					u_t_mat = (b_R*d_L - b_L*d_R)/detA;
					p_t_mat = (a_L*d_R - a_R*d_L)/detA;
					D0_p_tau = p_t_mat + U[2]*T_p;
					D0_u_tau = u_t_mat + U[2]*T_u;

					//already total D!
					D[1] = u_t_mat + (u_star-lambda_u)/U[0]/C/C * D0_p_tau + (u_star-lambda_u)*T_v;
					D[3] = p_t_mat + (u_star-lambda_u)*U[0] * D0_u_tau;
	
					if(u_star < lambda_u) //the direction is between the contact discontinuety and the 3-wave
						{
							if(CRW[1]) //the 3-wave is a CRW
								{
									//already total D!
									D[0] = rho_star_R*(u_star-lambda_u)*pow(c_star_R/c_R, (1.0+zetaR)/zetaR)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
									D[0] = (D[0] + D[3] + U[2]*T_p) / c_star_R/c_star_R - (U[2]-lambda_v)*T_rho;

									D[2] = -U[1]*d_v_R*U[0]/rho_R - (U[2]-lambda_v)*t_v_R - T_p/U[0];
									D[2] = D[2] + lambda_u*d_v_R;
									D[2] = D[2] + u_star/c_star_R*(zetaR-1.0)*(pow(c_frac, 2.0/zetaR-1.0)-1.0)/(zetaR-2.0)/U[0] * t_p_R;
									D[4] = -U[1]*d_z_R*U[0]/rho_R - (U[2]-lambda_v)*t_z_R;
									D[4] = D[4] + lambda_u*d_z_R;
									D[5] = -U[1]*d_phi_R*U[0]/rho_R - (U[2]-lambda_v)*t_phi_R;
									D[5] = D[5] + lambda_u*d_phi_R;
								}
							else //the 3-wave is a shock
								{
									SmUs = sqrt(0.5*((gammaR+1.0)*p_R   + (gammaR-1.0)*p_star)/rho_star_R);
									SmUR = sqrt(0.5*((gammaR+1.0)*p_star+ (gammaR-1.0)*p_R   )/rho_R);

									VAR = p_R + zetaR*p_star;
									H1 = rho_R * p_R    * (1.0 - zetaR*zetaR) / VAR/VAR;
									H2 = rho_R * p_star * (zetaR*zetaR - 1.0) / VAR/VAR;
									H3 = (p_star + zetaR*p_R)/VAR;

									L_rho = SmUR * H3 * d_rho_R;
									L_u = -rho_R * (H2*c_R*c_R + H3) * d_u_R;
									L_p = H2 * SmUR * d_p_R;
									L_v = -rho_R * (H2*c_R*c_R + H3) * t_v_R;

									D[0] = ((u_star+SmUs)/c_star_R/c_star_R - u_star*H1)*D0_p_tau + rho_star_R*u_star*SmUs*H1*D0_u_tau;
									D[0] = (D[0] - u_star*(L_p+L_rho+L_u+L_v)) / SmUs;
									D[0] = D[0] - (U[2]-lambda_v)*T_rho;

									f = SmUR*(H2*d_p_R + H3*d_rho_R) - rho_R*(H2*c_R*c_R+H3)*d_u_R;
									rho_x = (f + H1*(p_t_mat - rho_star_R*SmUs*u_t_mat) - D[0]) / (SmUR+u_R);//shk_spd;
									D[0] = D[0] + lambda_u*rho_x;

									D[2] = -(U[1]*(SmUR * d_v_R - t_p_R/rho_R)+(u_star+SmUs)*T_p/U[0]) / SmUs;
									D[2] = D[2] + lambda_u*d_v_R - (U[2]-lambda_v)*t_v_R;
									D[4] = -U[1] * SmUR * d_z_R / SmUs;
									D[4] = D[4] + lambda_u*d_z_R - (U[2]-lambda_v)*t_z_R;
									D[5] = -U[1] * SmUR * d_phi_R / SmUs;
									D[5] = D[5] + lambda_u*d_phi_R - (U[2]-lambda_v)*t_phi_R;
								}
						}
					else //the direction is between the 1-wave and the contact discontinuety
						{
							if(CRW[0]) //the 1-wave is a CRW
								{
									//already total D!
									D[0] = rho_star_L*(u_star-lambda_u)*pow(c_star_L/c_L, (1.0+zetaL)/zetaL)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
									D[0] = (D[0] + D[3] + U[2]*T_p) / c_star_L/c_star_L - (U[2]-lambda_v)*T_rho;

									D[2] = -U[1]*d_v_L*U[0]/rho_L - (U[2]-lambda_v)*t_v_L - T_p/U[0];
									D[2] = D[2] + lambda_u*d_v_L;
									D[2] = D[2] - u_star/c_star_L*(zetaL-1.0)*(pow(c_frac, 2.0/zetaL-1.0)-1.0)/(zetaL-2.0)/U[0] * t_p_L;
									D[4] = -U[1]*d_z_L*U[0]/rho_L - (U[2]-lambda_v)*t_z_L;
									D[4] = D[4] + lambda_u*d_z_L;
									D[5] = -U[1]*d_phi_L*U[0]/rho_L - (U[2]-lambda_v)*t_phi_L;
									D[5] = D[5] + lambda_u*d_phi_L;
								}
							else //the 1-wave is a shock
								{
									SmUs = -sqrt(0.5*((gammaL+1.0)*p_L   +(gammaL-1.0)*p_star)/rho_star_L);
									SmUL = -sqrt(0.5*((gammaL+1.0)*p_star+(gammaL-1.0)*p_L   )/rho_L);

									VAR = p_L + zetaL*p_star;
									H1 = rho_L * p_L    * (1.0 - zetaL*zetaL) / VAR/VAR;
									H2 = rho_L * p_star * (zetaL*zetaL - 1.0) / VAR/VAR;
									H3 = (p_star + zetaL*p_L)/VAR;

									L_rho = SmUL * H3 * d_rho_L;
									L_u = -rho_L*(H2*c_L*c_L + H3) * d_u_L;
									L_p = H2 * SmUL * d_p_L;
									L_v = -rho_L*(H2*c_L*c_L + H3) * t_v_L;

									D[0] = ((u_star+SmUs)/c_star_L/c_star_L - u_star*H1)*D0_p_tau + rho_star_L*u_star*SmUs*H1*D0_u_tau;
									D[0] = (D[0] - u_star*(L_p+L_rho+L_u+L_v))/ SmUs;
									D[0] = D[0] - (U[2]-lambda_v)*T_rho;

									f = SmUL*(H2*d_p_L + H3*d_rho_L) - rho_L*(H2*c_L*c_L+H3)*d_u_L;
									rho_x = (f + H1*(p_t_mat - rho_star_L*SmUs*u_t_mat) - D[0]) / (SmUL+u_L);
									D[0] = D[0] + lambda_u*rho_x;

									D[2] = -(U[1]*(SmUL * d_v_L - t_p_L/rho_L)+(u_star+SmUs)*T_p/U[0]) / SmUs;
									D[2] = D[2] + lambda_u*d_v_L - (U[2]-lambda_v)*t_v_L;
									D[4] = -U[1] * SmUL * d_z_L / SmUs;
									D[4] = D[4] + lambda_u*d_z_L - (U[2]-lambda_v)*t_z_L;
									D[5] = -U[1] * SmUL * d_phi_L / SmUs;
									D[5] = D[5] + lambda_u*d_phi_L - (U[2]-lambda_v)*t_phi_L;
								}
						}
					D[1] = D[1] + lambda_v*T_u;
					D[3] = D[3] + lambda_v*T_p;
					//--end of non-sonic case--
				}
			//----end of non-trivial case----
		}
	U_star[0] = rho_star_L;
	U_star[1] = u_star;
	U_star[2] = rho_star_R;
	U_star[3] = p_star;
	U_star[4] = c_star_L;
	U_star[5] = c_star_R;
}
