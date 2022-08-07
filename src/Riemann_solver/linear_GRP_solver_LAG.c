#include <math.h>
#include <stdio.h>

#include "../include/Riemann_solver.h"

/**
 * @param[in] atc: Parameter that determines the solver type.
 *              - INFINITY: acoustic approximation
 *              - eps:      GRP solver(nonlinear + acoustic case)
 *              - -0.0:     GRP solver(only nonlinear case)
*/

void linear_GRP_solver_LAG
(double * dire, double * mid,
 const double rho_L, const double rho_R, const double s_rho_L, const double s_rho_R,
 const double   u_L, const double   u_R, const double   s_u_L, const double   s_u_R,
 const double   p_L, const double   p_R, const double   s_p_L, const double   s_p_R,
 const double gamma, const double eps, const double  atc)
{
  const double zeta = (gamma-1.0)/(gamma+1.0);

  double dist;
  int CRW[2];

  double c_L, c_R, g_L, g_R;
  c_L = sqrt(gamma * p_L / rho_L);
  c_R = sqrt(gamma * p_R / rho_R);
  g_L = rho_L*c_L;
  g_R = rho_R*c_R;
  double speed_L, speed_R;
  double c_star_L, c_star_R, g_star_L, g_star_R;
  double u_star, p_star, rho_star_L, rho_star_R;
  double beta_star;

  double a_L, b_L, d_L, a_R, b_R, d_R, L_rho, L_u, L_p, A, B, sigma;  

  Riemann_solver_exact(&u_star, &p_star, gamma, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, eps, 500);

  if(CRW[0])
      {
	  rho_star_L = rho_L*pow(p_star/p_L, 1.0/gamma);
	  c_star_L = c_L*pow(p_star/p_L, 0.5*(gamma-1.0)/gamma);
	  speed_L = u_L - c_L;
      }
  else
      {
	  rho_star_L = rho_L*(p_star+zeta*p_L)/(p_L+zeta*p_star);
	  c_star_L = sqrt(gamma * p_star / rho_star_L);
	  speed_L = u_L - c_L*sqrt(0.5*((gamma+1.0)*(p_star/p_L) + (gamma-1.0))/gamma);
      }
  if(CRW[1])
      {
	  rho_star_R = rho_R*pow(p_star/p_R,1.0/gamma);
	  c_star_R = c_R*pow(p_star/p_R, 0.5*(gamma-1.0)/gamma);
	  speed_R = u_R + c_R;
      }
  else
      {
	  rho_star_R = rho_R*(p_star+zeta*p_R)/(p_R+zeta*p_star);
	  c_star_R = sqrt(gamma * p_star / rho_star_R);
	  speed_R = u_R + c_R*sqrt(0.5*((gamma+1.0)*(p_star/p_R) + (gamma-1.0))/gamma);
      }
  g_star_R = rho_star_R*c_star_R;
  g_star_L = rho_star_L*c_star_L;
 
  dist = sqrt((u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R));
  if(dist < atc) // acoustic Case
      {
	  a_L =  1.0;
	  b_L =  1.0 / g_star_L;
	  d_L = - g_L*s_u_L - s_p_L;
	
	  a_R = -1.0;
	  b_R =  1.0 / g_star_R;
	  d_R = - g_R*s_u_R + s_p_R;		   
      }
  else // nonlinear case
      {
	  //determine a_L, b_L and d_L
	  if(CRW[0]) //the 1-wave is a CRW
	      {
		  beta_star = g_star_L/g_L;
		  a_L = 1.0;
		  b_L = 1.0 / g_star_L;
		  d_L = (s_u_L+s_p_L/g_L) + 1.0/g_L/(3.0*gamma-1.0)*(c_L*c_L*s_rho_L-s_p_L)*(pow(beta_star,(3.0*gamma-1.0)/2.0/(gamma+1.0))-1.0); 
		  d_L = - 1.0 * sqrt(g_L*g_star_L)*d_L;
	      }
	  else //the 1-wave is a shock
	      {
		  
	      }
	  //determine a_R, b_R and d_R
	  if(CRW[1]) //the 3-wave is a CRW
	      {
		  beta_star = g_star_R/g_R;
		  a_R = -1.0;
		  b_R = 1.0 / g_star_R;
		  d_R = (s_u_R-s_p_R/g_R) + 1.0/g_R/(3.0*gamma-1.0)*(-c_L*c_L*s_rho_L+s_p_L)*(pow(beta_star,(3.0*gamma-1.0)/2.0/(gamma+1.0))-1.0);
		  d_R = - 1.0 * sqrt(g_R*g_star_R)*d_R;
	      }
	  else //the 3-wave is a shock
	      {
		  sigma = (p_star-p_R)/(u_star-u_R);
		  a_R = 2.0 - 0.5*(p_star-p_R)/(p_star+zeta*p_R);
		  b_R = -1.0*sigma/g_star_R/g_star_R - (a_R-1)/sigma;
		  L_rho = -1.0*(p_star-p_R)/2.0/rho_R;
		  L_u = sigma + rho_R*(u_star-u_R)/2.0 + g_R*g_R/sigma*(1.0+zeta*(p_star-p_R)/2.0/(p_star+zeta*p_R));
		  L_p = -1.0*(2+zeta/2.0*(p_star-p_R)/(p_star+zeta*p_R));
		  d_R = L_u*s_u_R + L_p*s_p_R + L_rho*s_rho_R;
	      }
      }
  
  mid[1] =   u_star;
  mid[2] =   p_star;
  mid[0] = rho_star_L;
  mid[3] = rho_star_R;
  dire[1] = (d_L*b_R-d_R*b_L)/(a_L*b_R-a_R*b_L);
  dire[2] = (d_L*a_R-d_R*a_L)/(b_L*a_R-b_R*a_L);
  dire[0] = 1.0/c_star_L/c_star_L*dire[2];
  dire[3] = 1.0/c_star_R/c_star_R*dire[2];
}
