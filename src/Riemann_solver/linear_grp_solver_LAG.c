/**
 * @file  linear_GRP_solver_LAG.c
 * @brief This is a Lagrangian GRP solver for compressible inviscid flow in Ben-Artzi's paper.
 */

#include <math.h>
#include <stdio.h>

#include "../include/var_struc.h"
#include "../include/riemann_solver.h"


/**
 * @brief A Lagrangian GRP solver for unsteady compressible inviscid two-component flow in one space dimension.
 * @param[out] D: the temporal derivative of fluid variables. \n
 *                   [rho_L, u, p, rho_R]_t
 * @param[out] U:  the Riemann solutions. \n
 *                   [rho_star_L, u_star, p_star, rho_star_R]
 * @param[in] ifv_L: Left  States (rho_L, u_L, p_L, s_rho_L, s_u_L, s_p_L, gammaL).
 * @param[in] ifv_R: Right States (rho_R, u_R, p_R, s_rho_R, s_u_R, s_p_R, gammaR).
 *                   - s_rho, s_u, s_p: ξ-Lagrangian spatial derivatives.
 *                   - gamma: the constant of the perfect gas.
 * @param[in] eps: the largest value could be seen as zero.
 * @param[in] atc: Parameter that determines the solver type.
 *              - INFINITY: acoustic approximation
 *              - eps:      GRP solver(nonlinear + acoustic case)
 *              - -0.0:     GRP solver(only nonlinear case)
 * @par  Reference
 *       Theory is found in Reference [1]. \n
 *       [1] M. Ben-Artzi & J. Falcovitz, A second-order Godunov-type scheme for compressible fluid dynamics,
 *           Journal of Computational Physics, 55.1: 1-32, 1984
 */
void linear_GRP_solver_LAG(double * D, double * U, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R, const double eps, const double  atc)
{
  const double   rho_L = ifv_L->RHO,     rho_R = ifv_R->RHO;
  const double s_rho_L = ifv_L->t_rho, s_rho_R = ifv_R->t_rho;
  const double     u_L = ifv_L->U,         u_R = ifv_R->U;
  const double   s_u_L = ifv_L->t_u,     s_u_R = ifv_R->t_u;
  const double     p_L = ifv_L->P,         p_R = ifv_R->P;
  const double   s_p_L = ifv_L->t_p,     s_p_R = ifv_R->t_p;
  const double  gammaL = ifv_L->gamma,  gammaR = ifv_R->gamma;

  const double zetaL = (gammaL-1.0)/(gammaL+1.0);
  const double zetaR = (gammaR-1.0)/(gammaR+1.0);

  double dist; // Euclidean distance
  _Bool CRW[2];  // Centred Rarefaction Wave (CRW) Indicator

  double c_L, c_R, g_L, g_R; // g = rho * c
  c_L = sqrt(gammaL * p_L / rho_L);
  c_R = sqrt(gammaR * p_R / rho_R);
  g_L = rho_L*c_L;
  g_R = rho_R*c_R;
  double W_L, W_R; // Wave speed
  double c_star_L, c_star_R, g_star_L, g_star_R;
  double u_star, p_star, rho_star_L, rho_star_R;
  double beta_star;

  double a_L, b_L, d_L, a_R, b_R, d_R, L_rho, L_u, L_p, A, B;

  Riemann_solver_exact(&u_star, &p_star, gammaL, gammaR, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, eps, 500);

  if(CRW[0])
      {
	  rho_star_L = rho_L*pow(p_star/p_L, 1.0/gammaL);
	  c_star_L = c_L*pow(p_star/p_L, 0.5*(gammaL-1.0)/gammaL);
	  W_L = u_L - c_L;
      }
  else
      {
	  rho_star_L = rho_L*(p_star+zetaL*p_L)/(p_L+zetaL*p_star);
	  c_star_L = sqrt(gammaL * p_star / rho_star_L);
	  W_L = u_L - c_L*sqrt(0.5*((gammaL+1.0)*(p_star/p_L) + (gammaL-1.0))/gammaL);
      }
  if(CRW[1])
      {
	  rho_star_R = rho_R*pow(p_star/p_R,1.0/gammaR);
	  c_star_R = c_R*pow(p_star/p_R, 0.5*(gammaR-1.0)/gammaR);
	  W_R = u_R + c_R;
      }
  else
      {
	  rho_star_R = rho_R*(p_star+zetaR*p_R)/(p_R+zetaR*p_star);
	  c_star_R = sqrt(gammaR * p_star / rho_star_R);
	  W_R = u_R + c_R*sqrt(0.5*((gammaR+1.0)*(p_star/p_R) + (gammaR-1.0))/gammaR);
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
		  d_L = (s_u_L+s_p_L/g_L) + 1.0/g_L/(3.0*gammaL-1.0)*(c_L*c_L*s_rho_L-s_p_L)*(pow(beta_star,(3.0*gammaL-1.0)/2.0/(gammaL+1.0))-1.0); 
		  d_L = - 1.0 * sqrt(g_L*g_star_L)*d_L;
	      }
	  else //the 1-wave is a shock
	      {
		  W_L = (p_star-p_L) / (u_star-u_L);
		  A   = - 0.5/(p_star + zetaL * p_L);
		  a_L = 2.0 + A * (p_star-p_L);
		  b_L = - W_L/g_star_L/g_star_L - (a_L - 1.0)/W_L;
		  L_rho = (p_star-p_L)/2.0/rho_L;
		  B = 1.0/(p_star-p_L) - zetaL * A;
		  L_u = rho_L * (u_star-u_L) * (gammaL*p_L*B + 0.5) + W_L;
		  L_p = 1.0 + B * (p_star-p_L);
		  d_L = L_u*s_u_L - L_p*s_p_L - L_rho*s_rho_L;		  
	      }
	  //determine a_R, b_R and d_R
	  if(CRW[1]) //the 3-wave is a CRW
	      {
		  beta_star = g_star_R/g_R;
		  a_R = -1.0;
		  b_R = 1.0 / g_star_R;
		  d_R = (s_u_R-s_p_R/g_R) + 1.0/g_R/(3.0*gammaR-1.0)*(-c_L*c_L*s_rho_L+s_p_L)*(pow(beta_star,(3.0*gammaR-1.0)/2.0/(gammaR+1.0))-1.0);
		  d_R = - 1.0 * sqrt(g_R*g_star_R)*d_R;
	      }
	  else //the 3-wave is a shock
	      {
		  W_R = (p_star-p_R) / (u_star-u_R);
		  A   = - 0.5/(p_star + zetaR * p_R);
		  a_R = - 2.0 - A * (p_star-p_R);
		  b_R = W_R/g_star_R/g_star_R - (a_R + 1.0)/W_R;
		  L_rho = (p_star-p_R)/2.0/rho_R;
		  B = 1.0/(p_star-p_R) - zetaR * A;
		  L_u = rho_R * (u_R-u_star) * (gammaR*p_R*B + 0.5) - W_R;
		  L_p = 1.0 + B * (p_star-p_R);
		  d_R = L_u*s_u_R + L_p*s_p_R + L_rho*s_rho_R;
	      }
      }
  
  U[1] =   u_star;
  U[2] =   p_star;
  U[0] = rho_star_L;
  U[3] = rho_star_R;
  D[1] = (d_L*b_R-d_R*b_L)/(a_L*b_R-a_R*b_L);
  D[2] = (d_L*a_R-d_R*a_L)/(b_L*a_R-b_R*a_L);
  D[0] = 1.0/c_star_L/c_star_L*D[2];
  D[3] = 1.0/c_star_R/c_star_R*D[2];
}
