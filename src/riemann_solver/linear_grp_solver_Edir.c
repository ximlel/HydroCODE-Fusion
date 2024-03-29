/**
 * @file  linear_GRP_solver_Edir.c
 * @brief This is a direct Eulerian GRP solver for compressible inviscid flow in Li's paper.
 */

#include <math.h>
#include <stdio.h>

#include "../include/var_struc.h"
#include "../include/riemann_solver.h"


/**
 * @brief A direct Eulerian GRP solver for unsteady compressible inviscid flow in one space dimension.
 * @param[out] D: the temporal derivative of fluid variables. \n
 *                   [rho, u, p]_t
 * @param[out] U:  the intermediate Riemann solutions at t-axis. \n
 *                   [rho_mid, u_mid, p_mid]
 * @param[in] ifv_L: Left  States (rho_L, u_L, p_L, s_rho_L, s_u_L, s_p_L, gamma).
 * @param[in] ifv_R: Right States (rho_R, u_R, p_R, s_rho_R, s_u_R, s_p_R).
 *                   - s_rho, s_u, s_p: x-spatial derivatives.
 *                   - gamma: the constant of the perfect gas.
 * @param[in] eps: the largest value could be seen as zero.
 * @param[in] atc: Parameter that determines the solver type.
 *              - INFINITY: acoustic approximation
 *                - ifv_.s_ = -0.0: exact Riemann solver 
 *              - eps:      1D GRP solver(nonlinear + acoustic case)
 *              - -0.0:     1D GRP solver(only nonlinear case)
 * @sa   Theory is found in Reference [1]. \n
 *       [1] M. Ben-Artzi, J. Li & G. Warnecke, A direct Eulerian GRP scheme for compressible fluid flows.
 *           Journal of Computational Physics, 218.1: 19-43, 2006.
 */
void linear_GRP_solver_Edir(double * D, double * U, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R, const double eps, const double atc)
{
  const double   rho_L = ifv_L->RHO,     rho_R = ifv_R->RHO;
  const double s_rho_L = ifv_L->d_rho, s_rho_R = ifv_R->d_rho;
  const double     u_L = ifv_L->U,         u_R = ifv_R->U;
  const double   s_u_L = ifv_L->d_u,     s_u_R = ifv_R->d_u;
  const double     p_L = ifv_L->P,         p_R = ifv_R->P;
  const double   s_p_L = ifv_L->d_p,     s_p_R = ifv_R->d_p;
  const double   gamma = ifv_L->gamma;

  double dist;
  double c_L, c_R;
  _Bool CRW[2];
  double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R;

  double PI, H1, H2, H3;
  double a_L, b_L, d_L, a_R, b_R, d_R;
  double L_u, L_p, L_rho;
  double u_t_mat, p_t_mat;
  double shk_spd, zeta = (gamma-1.0)/(gamma+1.0), zts = zeta*zeta;
  double g_rho, g_u, g_p, f;
  double speed_L, speed_R;

  c_L = sqrt(gamma * p_L / rho_L);
  c_R = sqrt(gamma * p_R / rho_R);

  dist = sqrt((u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R));
  if (dist < atc && atc < 2*eps)
      {
	  rho_star_L = rho_L;
	  rho_star_R = rho_R;
	  c_star_L = c_L;
	  c_star_R = c_R;
	  u_star = 0.5*(u_R+u_L);
	  p_star = 0.5*(p_R+p_L);
      }
  else
      {
	  Riemann_solver_exact_single(&u_star, &p_star, gamma, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, eps, 50);

	  if(p_star > p_L)
	      rho_star_L = rho_L*(p_star+zeta*p_L)/(p_L+zeta*p_star);
	  else
	      rho_star_L = rho_L*pow(p_star/p_L,1.0/gamma);
	  if(p_star > p_R)
	      rho_star_R = rho_R*(p_star+zeta*p_R)/(p_R+zeta*p_star);
	  else
	      rho_star_R = rho_R*pow(p_star/p_R,1.0/gamma);
	  c_star_L = sqrt(gamma * p_star / rho_star_L);
	  c_star_R = sqrt(gamma * p_star / rho_star_R);
      }

//=========acoustic case==========
  if(dist < atc)
  {
    //------trivial case------
    if(u_L-c_L > 0.0) //the t-axe is on the left side of all the three waves
    {
      D[0] = -s_rho_L*u_L - rho_L*s_u_L;
      D[1] = (D[0]*u_L + s_rho_L*u_L*u_L + 2.0*rho_L*u_L*s_u_L + s_p_L) / -rho_L;
      D[2] = -(gamma-1.0) * (0.5*D[0]*u_L*u_L + rho_L*u_L*D[1]);
      D[2] = D[2] - s_u_L * (gamma*p_L + 0.5*(gamma-1.0)*rho_L*u_L*u_L);
      D[2] = D[2] - u_L * (gamma * s_p_L + (gamma-1.0)*(0.5*s_rho_L*u_L*u_L + rho_L*u_L*s_u_L));

      U[0] = rho_L;
      U[1] =   u_L;
      U[2] =   p_L;
    }
    else if(u_R+c_R < 0.0) //the t-axe is on the right side of all the three waves
    {
      D[0] = -s_rho_R*u_R - rho_R*s_u_R;
      D[1] = (D[0]*u_R + s_rho_R*u_R*u_R + 2.0*rho_R*u_R*s_u_R + s_p_R) / -rho_R;
      D[2] = -(gamma-1.0) * (0.5*D[0]*u_R*u_R + rho_R*u_R*D[1]);
      D[2] = D[2] - s_u_R * (gamma*p_R + 0.5*(gamma-1.0)*rho_R*u_R*u_R);
      D[2] = D[2] - u_R * (gamma * s_p_R + (gamma-1.0)*(0.5*s_rho_R*u_R*u_R + rho_R*u_R*s_u_R));

      U[0] = rho_R;
      U[1] =   u_R;
      U[2] =   p_R;
    }
    //------non-trivial case------
    else
    {
      if(u_star > 0.0)
      {
	U[0] = rho_star_L;
	U[1] =   u_star;
	U[2] =   p_star;

	PI = (u_star+c_star_R)*rho_star_L*c_star_L*c_star_L / (u_star-c_star_L)/rho_star_R/c_star_R/c_star_R;
	D[1] = (s_p_L/rho_L+c_L*s_u_L)*PI/(1.0-PI) + (s_p_R/rho_R-c_R*s_u_R)/(PI-1.0);
	D[2] = ((u_star+c_star_R)/rho_star_R/c_star_R/c_star_R) - ((u_star-c_star_L)/rho_star_L/c_star_L/c_star_L);
	D[2] = (s_p_R/rho_R-c_R*s_u_R-s_p_L/rho_L-c_L*s_u_L) / D[2];
	D[2] = D[2] * (1.0 - (u_star*u_star/c_star_L/c_star_L)) + rho_star_L*u_star*D[1];
	D[0] = (u_star*(s_p_L - s_rho_L*c_star_L*c_star_L) + D[2])/c_star_L/c_star_L;
      }
      else
      {
	U[0] = rho_star_R;
	U[1] =   u_star;
	U[2] =   p_star;

	PI = (u_star+c_star_R)*rho_star_L*c_star_L*c_star_L / (u_star-c_star_L)/rho_star_R/c_star_R/c_star_R;
	D[1] = (s_p_L/rho_L+c_L*s_u_L)*PI/(1.0-PI) + (s_p_R/rho_R-c_R*s_u_R)/(PI-1.0);
	D[2] = ((u_star+c_star_R)/rho_star_R/c_star_R/c_star_R) - ((u_star-c_star_L)/rho_star_L/c_star_L/c_star_L);
	D[2] = (s_p_R/rho_R-c_R*s_u_R-s_p_L/rho_L-c_L*s_u_L) / D[2];
	D[2] = D[2] * (1.0 - (u_star*u_star/c_star_R/c_star_R)) + rho_star_R*u_star*D[1];
	D[0] = (u_star*(s_p_R - s_rho_R*c_star_R*c_star_R) + D[2])/c_star_R/c_star_R;
      }
    }
    return;
  }

//=========non-acoustic case==========
//----------solving the LINEAR GRP----------
  if(CRW[0])
    speed_L = u_L - c_L;
  else
    speed_L = (rho_star_L*u_star - rho_L*u_L) / (rho_star_L - rho_L);
  if(CRW[1])
    speed_R = u_R + c_R;
  else
    speed_R = (rho_star_R*u_star - rho_R*u_R) / (rho_star_R - rho_R);

  //------trivial case------
  if(speed_L > 0.0) //the t-axe is on the left side of all the three waves
  {
    D[0] = -s_rho_L*u_L - rho_L*s_u_L;
    D[1] = (D[0]*u_L + s_rho_L*u_L*u_L + 2.0*rho_L*u_L*s_u_L + s_p_L) / -rho_L;
    D[2] = (s_u_L*p_L + u_L*s_p_L)*gamma/(1.0-gamma) - 0.5*s_rho_L*u_L*u_L*u_L - 1.5*rho_L*u_L*u_L*s_u_L;
    D[2] = D[2] - 0.5*D[0]*u_L*u_L - rho_L*u_L*D[1];
    D[2] = D[2] * (gamma-1.0);

    U[0] = rho_L;
    U[1] =   u_L;
    U[2] =   p_L;
  }
  else if(speed_R < 0.0) //the t-axe is on the right side of all the three waves
  {
    D[0] = -s_rho_R*u_R - rho_R*s_u_R;
    D[1] = (D[0]*u_R + s_rho_R*u_R*u_R + 2.0*rho_R*u_R*s_u_R + s_p_R) / -rho_R;
    D[2] = -(gamma-1.0) * (0.5*D[0]*u_R*u_R + rho_R*u_R*D[1]);
    D[2] = D[2] - s_u_R * (gamma*p_R + 0.5*(gamma-1.0)*rho_R*u_R*u_R);
    D[2] = D[2] - u_R * (gamma * s_p_R + (gamma-1.0)*(0.5*s_rho_R*u_R*u_R + rho_R*u_R*s_u_R));

    U[0] = rho_R;
    U[1] =   u_R;
    U[2] =   p_R;
  }
  //----non-trivial case----
  else
  {
    if((CRW[0]) && ((u_star-c_star_L) > 0.0)) // the t-axe is in a 1-CRW
    {
      // shk_spd = (rho_star_L*u_star - rho_L*u_L)/(rho_star_L - rho_L);

      U[1] = zeta*(u_L+2.0*c_L/(gamma-1.0));
      U[2] = U[1]*U[1]*rho_L/gamma/pow(p_L, 1.0/gamma);
      U[2] = pow(U[2], gamma/(gamma-1.0));
      U[0] = gamma*U[2]/U[1]/U[1];

      D[1] = 0.5*(pow(U[1]/c_L, 0.5/zeta)*(1.0+zeta) + pow(U[1]/c_L, (1.0+zeta)/zeta)*zeta)/(0.5+zeta);
      D[1] = D[1] * (s_p_L - s_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
      D[1] = D[1] - c_L*pow(U[1]/c_L, 0.5/zeta)*(s_u_L + (gamma*s_p_L/c_L - c_L*s_rho_L)/(gamma-1.0)/rho_L);

      D[2] = U[0]*U[1]*D[1];

      D[0] = U[0]*U[1]*pow(U[1]/c_L, (1.0+zeta)/zeta)*(s_p_L - s_rho_L*c_L*c_L)/rho_L;
      D[0] = (D[0] + D[2]) / U[1]/U[1];
    }
    else if((CRW[1]) && ((u_star+c_star_R) < 0.0)) // the t-axe is in a 3-CRW
    {
      // shk_spd = (rho_star_R*u_star - rho_R*u_R)/(rho_star_R - rho_R);

      U[1] = zeta*(u_R-2.0*c_R/(gamma-1.0));
      U[2] = U[1]*U[1]*rho_R/gamma/pow(p_R, 1.0/gamma);
      U[2] = pow(U[2], gamma/(gamma-1.0));
      U[0] = gamma*U[2]/U[1]/U[1];

      D[1] = 0.5*(pow(-U[1]/c_R, 0.5/zeta)*(1.0+zeta) + pow(-U[1]/c_R, (1.0+zeta)/zeta)*zeta)/(0.5+zeta);
      D[1] = D[1] * (s_p_R - s_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
      D[1] = D[1] + c_R*pow(-U[1]/c_R, 0.5/zeta)*(s_u_R - (gamma*s_p_R/c_R - c_R*s_rho_R)/(gamma-1.0)/rho_R);

      D[2] = U[0]*U[1]*D[1];

      D[0] = U[0]*U[1]*pow(-U[1]/c_R, (1.0+zeta)/zeta)*(s_p_R - s_rho_R*c_R*c_R)/rho_R;
      D[0] = (D[0] + D[2]) / U[1]/U[1];
    }
    //--non-sonic case--
    else
    {
    //determine a_L, b_L and d_L
      if(CRW[0]) //the 1-wave is a CRW
      {
	a_L = 1.0;
        b_L = 1.0 / rho_star_L / c_star_L;
	d_L = 0.5*(pow(c_star_L/c_L, 0.5/zeta)*(1.0+zeta) + pow(c_star_L/c_L, (1.0+zeta)/zeta)*zeta)/(0.5+zeta);
	d_L = d_L * (s_p_L - s_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
	d_L = d_L - c_L*pow(c_star_L/c_L, 0.5/zeta)*(s_u_L + (gamma*s_p_L/c_L - c_L*s_rho_L)/(gamma-1.0)/rho_L);
      }
      else //the 1-wave is a shock
      {
	H1 = 0.5*sqrt((1.0-zeta)/(rho_L*(p_star+zeta*p_L))) * (p_star + (1.0+2.0*zeta)*p_L)/(p_star+zeta*p_L);
	H2 = -0.5*sqrt((1.0-zeta)/(rho_L*(p_star+zeta*p_L))) * ((2.0+zeta)*p_star + zeta*p_L)/(p_star+zeta*p_L);
	H3 = -0.5*sqrt((1.0-zeta)/(rho_L*(p_star+zeta*p_L))) * (p_star-p_L) / rho_L;
	shk_spd = (rho_star_L*u_star - rho_L*u_L)/(rho_star_L - rho_L);

	a_L = 1.0 - rho_star_L*(shk_spd-u_star)*H1;
	b_L = (u_star - shk_spd)/rho_star_L/c_star_L/c_star_L + H1;

	L_rho = (u_L-shk_spd) * H3;
	L_u = shk_spd - u_L + rho_L*c_L*c_L*H2 + rho_L*H3;
	L_p = (u_L-shk_spd)*H2 - 1.0/rho_L;

	d_L = L_rho*s_rho_L + L_u*s_u_L + L_p*s_p_L;
      }
    //determine a_R, b_R and d_R
      if(CRW[1]) //the 3-wave is a CRW
      {
	a_R = 1.0;
        b_R = -1.0 / rho_star_R / c_star_R;
	d_R = 0.5*(pow(c_star_R/c_R, 0.5/zeta)*(1.0+zeta) + pow(c_star_R/c_R, (1.0+zeta)/zeta)*zeta)/(0.5+zeta);
	d_R = d_R * (s_p_R - s_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
	d_R = d_R + c_R*pow(c_star_R/c_R, 0.5/zeta)*(s_u_R - (gamma*s_p_R/c_R - c_R*s_rho_R)/(gamma-1.0)/rho_R);
      }
      else //the 3-wave is a shock
      {
	H1 = 0.5*sqrt((1.0-zeta)/(rho_R*(p_star+zeta*p_R))) * (p_star + (1.0+2.0*zeta)*p_R)/(p_star+zeta*p_R);
	H2 = -0.5*sqrt((1.0-zeta)/(rho_R*(p_star+zeta*p_R))) * ((2.0+zeta)*p_star + zeta*p_R)/(p_star+zeta*p_R);
	H3 = -0.5*sqrt((1.0-zeta)/(rho_R*(p_star+zeta*p_R))) * (p_star-p_R) / rho_R;
	shk_spd = (rho_star_R*u_star - rho_R*u_R)/(rho_star_R - rho_R);

	a_R = 1.0 + rho_star_R*(shk_spd-u_star)*H1;
	b_R = (u_star - shk_spd)/rho_star_R/c_star_R/c_star_R - H1;

	L_rho = (shk_spd-u_R) * H3;
	L_u = shk_spd - u_R - rho_R*c_R*c_R*H2 - rho_R*H3;
	L_p = (shk_spd-u_R)*H2 - 1.0/rho_R;

	d_R = L_rho*s_rho_R + L_u*s_u_R + L_p*s_p_R;
      }

      p_t_mat = (d_L*a_R/a_L-d_R)/(b_L*a_R/a_L-b_R);
      u_t_mat = (d_L - b_L*p_t_mat)/a_L;

      if(u_star < 0.0) //the t-axi is between the contact discontinuety and the 3-wave
      {
	U[0] = rho_star_R;
	U[1] =   u_star;
	U[2] =   p_star;
        D[1] = u_t_mat + u_star*p_t_mat/rho_star_R/c_star_R/c_star_R;
        D[2] = p_t_mat + rho_star_R*u_star * u_t_mat;

	if(CRW[1]) //the 3-wave is a CRW
	{
	  D[0] = rho_star_R*u_star*pow(c_star_R/c_R, (1.0+zeta)/zeta)*(s_p_R - s_rho_R*c_R*c_R)/rho_R;
	  D[0] = (D[0] + D[2]) / c_star_R/c_star_R;
	}
	else //the 3-wave is a shock
	{
	  shk_spd = (rho_star_R*u_star - rho_R*u_R)/(rho_star_R - rho_R);
	  H1 = rho_R * p_R    * (1.0 - zts) / (p_R + zeta*p_star) / (p_R + zeta*p_star);
	  H2 = rho_R * p_star * (zts - 1.0) / (p_R + zeta*p_star) / (p_R + zeta*p_star);
	  H3 = (p_star + zeta*p_R) / (p_R + zeta*p_star);

	  g_rho = u_star-shk_spd;
	  g_u   = u_star*rho_star_R*(shk_spd-u_star)*H1;
	  g_p   = shk_spd/c_star_R/c_star_R - u_star*H1;
	  f = (shk_spd-u_R)*(H2*s_p_R + H3*s_rho_R) - rho_R*(H2*c_R*c_R+H3)*s_u_R;

	  D[0] = (f*u_star - g_p*p_t_mat - g_u*u_t_mat) / g_rho;
	}
      }
      else //the t-axi is between the 1-wave and the contact discontinuety
      {
	U[0] = rho_star_L;
	U[1] =   u_star;
	U[2] =   p_star;
        D[1] = u_t_mat + u_star*p_t_mat/rho_star_L/c_star_L/c_star_L;
        D[2] = p_t_mat + rho_star_L*u_star * u_t_mat;
	if(CRW[0]) //the 1-wave is a CRW
	{
	  D[0] = rho_star_L*u_star*pow(c_star_L/c_L, (1.0+zeta)/zeta)*(s_p_L - s_rho_L*c_L*c_L)/rho_L;
	  D[0] = (D[0] + D[2]) / c_star_L/c_star_L;
	}
	else //the 1-wave is a shock
	{
	  shk_spd = (rho_star_L*u_star - rho_L*u_L)/(rho_star_L - rho_L);
	  H1 = rho_L * p_L    * (1.0 - zts) / (p_L + zeta*p_star) / (p_L + zeta*p_star);
	  H2 = rho_L * p_star * (zts - 1.0) / (p_L + zeta*p_star) / (p_L + zeta*p_star);
	  H3 = (p_star + zeta*p_L) / (p_L + zeta*p_star);

	  g_rho = u_star-shk_spd;
	  g_u   = u_star*rho_star_L*(shk_spd-u_star)*H1;
	  g_p   = shk_spd/c_star_L/c_star_L - u_star*H1;
	  f = (shk_spd-u_L)*(H2*s_p_L + H3*s_rho_L) - rho_L*(H2*c_L*c_L+H3)*s_u_L;

	  D[0] = (f*u_star - g_p*p_t_mat - g_u*u_t_mat) / g_rho;
	}
      }
    //--end of non-sonic case--
    }
  //----end of non-trivial case----
  }
}
