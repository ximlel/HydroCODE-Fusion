#include <math.h>
#include <stdio.h>

#include "Riemann_solver.h"



void linear_GRP_solver_LAG
(double * direvative, double * mid,
 double rho_L, double rho_R, double s_rho_L, double s_rho_R,
 double   u_L, double   u_R, double   s_u_L, double   s_u_R,
 double   p_L, double   p_R, double   s_p_L, double   s_p_R,
 double gamma, double eps)
{
  double dist;
  int CRW[2];

  double c_L, c_R, g_L, g_R;
  double rho_star, u_star, p_star,  rho_star_L, rho_star_R, c_star_L, c_star_R, g_star_L, g_star_R, beta_star;

  double a_L, b_L, d_L, a_R, b_R, d_R, L_rho, L_u, L_p, sigma;  

  double zeta = (gamma-1.0)/(gamma+1.0);

  c_L = sqrt(gamma * p_L / rho_L);
  c_R = sqrt(gamma * p_R / rho_R);

  g_L = rho_L*c_L;
  g_R = rho_R*c_R;

  Riemann_solver_exact(&u_star, &p_star, gamma, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, 500);
  
  if(p_star > p_L)
    rho_star_L = rho_L*(p_star+zeta*p_L)/(p_L+zeta*p_star);
  else
    rho_star_L = rho_L*pow(p_star/p_L,1.0/gamma);
  if(p_star > p_R)
    rho_star_R = rho_R*(p_star+zeta*p_R)/(p_R+zeta*p_star);
  else
    rho_star_R = rho_R*pow(p_star/p_R,1.0/gamma);

    mid[1] =   u_star;
    mid[2] =   p_star;
    mid[0] = rho_star_L;
    mid[3] = rho_star_R;
  
  c_star_L = sqrt(gamma * p_star / rho_star_L);
  c_star_R = sqrt(gamma * p_star / rho_star_R);

  g_star_R = rho_star_R*c_star_R;
  g_star_L = rho_star_L*c_star_L;

  beta_star = g_star_L/g_L;
 
  dist = sqrt( (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R));

  if((dist < eps*1000))//Acoustic Case
	  {
		  a_L = 1.0;
		  b_L = 1.0 / g_star_L;
		  d_L = -1.0/g_star_L*g_L*(g_L*s_u_L + s_p_L);

	
		  a_R = -1.0;
		  b_R = 1.0 / g_star_R;
		  d_R = -1.0/g_star_R*g_R*(g_R*s_u_R - s_p_R);		   
	  }
  else
	  {
		  a_L = 1.0;
		  b_L = 1.0/ g_star_L;
		  d_L = (s_u_L+s_p_L/g_L) + 1.0/g_L/(3.0*gamma-1.0)*(c_L*c_L*s_rho_L-s_p_L)*(pow(beta_star,(3.0*gamma-1.0)/2.0/(gamma+1.0))-1.0); 
		  d_L = (-1.0)*sqrt(g_L*g_star_L)*d_L;


		  sigma = (p_star-p_R)/(u_star-u_R);
		  a_R = 2.0 - 0.5*(p_star-p_R)/(p_star+zeta*p_R);
		  b_R = -1.0*sigma/g_star_R/g_star_R - (a_R-1)/sigma;
		  L_rho = -1.0*(p_star-p_R)/2.0/rho_R;
		  L_u = sigma + rho_R*(u_star-u_R)/2.0 + g_R*g_R/sigma*(1.0+zeta*(p_star-p_R)/2.0/(p_star+zeta*p_R));
		  L_p = -1.0*(2+zeta/2.0*(p_star-p_R)/(p_star+zeta*p_R));
		  d_R = L_u*s_u_R + L_p*s_p_R + L_rho*s_rho_R;		 
	 }

  
  direvative[1] = (d_L*b_R-d_R*b_L)/(a_L*b_R-a_R*b_L);
  direvative[2] = (d_L*a_R-d_R*a_L)/(b_L*a_R-b_R*a_L);
  direvative[0] = 1.0/c_star_L/c_star_L*direvative[2];
  direvative[3] = 1.0/c_star_R/c_star_R*direvative[2];

}








