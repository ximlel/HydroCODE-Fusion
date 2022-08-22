#include <stdio.h>

#include "../include/var_struc.h"
#include "../include/Riemann_solver.h"


void flux_generator_y
(double h_y, int m, int n, double tau, double gamma,
 double *s_rho[], double *rho[m],
 double   *s_u[], double   *u[m],
 double   *s_v[], double   *v[m],
 double   *s_p[], double   *p[m],
 double *G1[m], double *G2[m], double *G3[m], double *G4[m],
 double *rhoy[m], double *uy[m], double *vy[m], double *py[m],
 double eps
 )
{
  int i = 0, j = 0;
  double rho_l, rho_r;
  double u_l, v_r;
  double v_l, u_r;
  double p_l, p_r;
  double s_rho_l, s_rho_r;
  double s_u_l, s_v_r;
  double s_v_l, s_u_r;
  double s_p_l, s_p_r;
  double dire[4], mid[4];


//===========================
  for(j = 0; j < m; ++j)
    for(i = 0; i <= n; ++i)
    {
      if(i)
      {
	rho_l = rho[j][i-1] + 0.5*h_y*s_rho[j][i-1];
          u_l =   u[j][i-1] + 0.5*h_y*  s_u[j][i-1];
          v_l =   v[j][i-1] + 0.5*h_y*  s_v[j][i-1];
          p_l =   p[j][i-1] + 0.5*h_y*  s_p[j][i-1];
      }
      else
      {
	rho_l = rho[j][i] + 0.5*h_y*s_rho[j][i];
          u_l =   u[j][i] + 0.5*h_y*  s_u[j][i];
          v_l =   v[j][i] + 0.5*h_y*  s_v[j][i];
          p_l =   p[j][i] + 0.5*h_y*  s_p[j][i];
      }
      if(i < n)
      {
	rho_r = rho[j][i] - 0.5*h_y*s_rho[j][i];
          u_r =   u[j][i] - 0.5*h_y*  s_u[j][i];
          v_r =   v[j][i] - 0.5*h_y*  s_v[j][i];
          p_r =   p[j][i] - 0.5*h_y*  s_p[j][i];
      }
      else
      {
	rho_r = rho[j][i-1] - 0.5*h_y*s_rho[j][i-1];
          u_r =   u[j][i-1] - 0.5*h_y*  s_u[j][i-1];
          v_r =   v[j][i-1] - 0.5*h_y*  s_v[j][i-1];
          p_r =   p[j][i-1] - 0.5*h_y*  s_p[j][i-1];
      }
//===========================
      if(i)
      {
	s_rho_l = s_rho[j][i-1];
          s_u_l =   s_u[j][i-1];
          s_v_l =   s_v[j][i-1];
          s_p_l =   s_p[j][i-1];
      }
      else
      {
	s_rho_l = s_rho[j][i];
          s_u_l =   s_u[j][i];
          s_v_l =   s_v[j][i];
          s_p_l =   s_p[j][i];
      }
      if(i < n)
      {
	s_rho_r = s_rho[j][i];
          s_u_r =   s_u[j][i];
          s_v_r =   s_v[j][i];
          s_p_r =   s_p[j][i];
      }
      else
      {
	s_rho_r = s_rho[j][i-1];
          s_u_r =   s_u[j][i-1];
          s_v_r =   s_v[j][i-1];
          s_p_r =   s_p[j][i-1];
      }
//===========================

      //rho, v, u, p
      linear_GRP_solver_Edir(dire, mid, rho_l, rho_r, s_rho_l, s_rho_r, v_l, v_r, s_v_l, s_v_r, u_l, u_r, s_u_l, s_u_r, p_l, p_r, s_p_l, s_p_r, gamma, eps);

      G1[j][i] = (mid[0]+0.5*tau*dire[0])*(mid[1]+0.5*tau*dire[1]);
      G2[j][i] = (mid[0]+0.5*tau*dire[0])*(mid[2]+0.5*tau*dire[2])*(mid[1]+0.5*tau*dire[1]);
      G3[j][i] = (mid[0]+0.5*tau*dire[0])*(mid[1]+0.5*tau*dire[1])*(mid[1]+0.5*tau*dire[1]) + (mid[3]+0.5*tau*dire[3]);
      G4[j][i] = 0.5*(mid[0]+0.5*tau*dire[0])*(mid[1]+0.5*tau*dire[1])*(mid[1]+0.5*tau*dire[1]);
      G4[j][i] = (gamma*(mid[3]+0.5*tau*dire[3])/(gamma-1.0) + G4[j][i])*(mid[1]+0.5*tau*dire[1]);

      rhoy[j][i] = mid[0] + tau*dire[0];
        uy[j][i] = mid[2] + tau*dire[2];
        vy[j][i] = mid[1] + tau*dire[1];
        py[j][i] = mid[3] + tau*dire[3];
    }
}
