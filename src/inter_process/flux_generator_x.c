#include <stdio.h>

#include "../include/var_struc.h"
#include "../include/Riemann_solver.h"


void flux_generator_x
(double h_x, int m, int n, double tau, double gamma,
 double *s_rho[], double *rho[m],
 double   *s_u[], double   *u[m],
 double   *s_v[], double   *v[m],
 double   *s_p[], double   *p[m],
 double *F1[m+1], double *F2[m+1], double *F3[m+1], double *F4[m+1],
 double *rhox[m+1], double *ux[m+1], double *vx[m+1], double *px[m+1],
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
  struct i_f_var ifv_L, ifv_R;
  /*
   * dire: the temporal derivative of fluid variables.
   *       \frac{\partial [rho, u, v, p]}{\partial t}
   * mid:  the Riemann solutions.
   *       [rho_star, u_star, v_star, p_star]
   */
  double dire[4], mid[4];
  double dire2[3], mid2[3];
  double f1, f2, f3;
  double errd = 0.0, errm = 0.0, errf1, errf2, errf3;

//===========================
  for(i = 0; i < n; ++i)
    for(j = 0; j <= m; ++j)
    {
      if(j)
      {
	rho_l = rho[j-1][i] + 0.5*h_x*s_rho[j-1][i];
          u_l =   u[j-1][i] + 0.5*h_x*  s_u[j-1][i];
          v_l =   v[j-1][i] + 0.5*h_x*  s_v[j-1][i];
          p_l =   p[j-1][i] + 0.5*h_x*  s_p[j-1][i];
      }
      else
      {
	rho_l = rho[j][i] + 0.5*h_x*s_rho[j][i];
          u_l =   u[j][i] + 0.5*h_x*  s_u[j][i];
          v_l =   v[j][i] + 0.5*h_x*  s_v[j][i];
          p_l =   p[j][i] + 0.5*h_x*  s_p[j][i];
      }
      if(j < m)
      {
	rho_r = rho[j][i] - 0.5*h_x*s_rho[j][i];
          u_r =   u[j][i] - 0.5*h_x*  s_u[j][i];
          v_r =   v[j][i] - 0.5*h_x*  s_v[j][i];
          p_r =   p[j][i] - 0.5*h_x*  s_p[j][i];
      }
      else
      {
	rho_r = rho[j-1][i] - 0.5*h_x*s_rho[j-1][i];
          u_r =   u[j-1][i] - 0.5*h_x*  s_u[j-1][i];
          v_r =   v[j-1][i] - 0.5*h_x*  s_v[j-1][i];
          p_r =   p[j-1][i] - 0.5*h_x*  s_p[j-1][i];
      }
//===========================
      if(j)
      {
	s_rho_l = s_rho[j-1][i];
          s_u_l =   s_u[j-1][i];
          s_v_l =   s_v[j-1][i];
          s_p_l =   s_p[j-1][i];
      }
      else
      {
	s_rho_l = s_rho[j][i];
          s_u_l =   s_u[j][i];
          s_v_l =   s_v[j][i];
          s_p_l =   s_p[j][i];
      }
      if(j < m)
      {
	s_rho_r = s_rho[j][i];
          s_u_r =   s_u[j][i];
          s_v_r =   s_v[j][i];
          s_p_r =   s_p[j][i];
      }
      else
      {
	s_rho_r = s_rho[j-1][i];
          s_u_r =   s_u[j-1][i];
          s_v_r =   s_v[j-1][i];
          s_p_r =   s_p[j-1][i];
      }
//===========================

      linear_GRP_solver_Edir(dire, mid, rho_l, rho_r, s_rho_l, s_rho_r, u_l, u_r, s_u_l, s_u_r, v_l, v_r, s_v_l, s_v_r, p_l, p_r, s_p_l, s_p_r, gamma, eps);


      F1[j][i] = (mid[0]+0.5*tau*dire[0])*(mid[1]+0.5*tau*dire[1]);
      F2[j][i] = (mid[0]+0.5*tau*dire[0])*(mid[1]+0.5*tau*dire[1])*(mid[1]+0.5*tau*dire[1]) + (mid[3]+0.5*tau*dire[3]);
      F3[j][i] = (mid[0]+0.5*tau*dire[0])*(mid[1]+0.5*tau*dire[1])*(mid[2]+0.5*tau*dire[2]);
      //F4[j][i] = 0.5*(mid[0]+0.5*tau*dire[0])*((mid[1]+0.5*tau*dire[1])*(mid[1]+0.5*tau*dire[1]) + (mid[2]+0.5*tau*dire[2])*(mid[2]+0.5*tau*dire[2]));
      F4[j][i] = 0.5*(mid[0]+0.5*tau*dire[0])*(mid[1]+0.5*tau*dire[1])*(mid[1]+0.5*tau*dire[1]);
      F4[j][i] = (gamma*(mid[3]+0.5*tau*dire[3])/(gamma-1.0) + F4[j][i])*(mid[1]+0.5*tau*dire[1]);

      /*
      f1 = (mid[0]+0.5*tau*dire[0])*(mid[1]+0.5*tau*dire[1]);
      f2 = f1*(mid[1]+0.5*tau*dire[1]) + mid[3]+0.5*tau*dire[3];
      f3 = (gamma/(gamma-1.0))*(mid[3]+0.5*tau*dire[3]) + 0.5*f1*(mid[1]+0.5*tau*dire[1]);
      f3 = f3*(mid[1]+0.5*tau*dire[1]);

      errf1 = (F1[j][i]-f1)*(F1[j][i]-f1);
      errf2 = (F2[j][i]-f2)*(F2[j][i]-f2);
      errf3 = (F4[j][i]-f3)*(F4[j][i]-f3);
      if((errf1 > eps))
	printf("sigularity f1\n");
      if((errf2 > eps))
	printf("sigularity f2\n");
      if((errf3 > eps))
	printf("sigularity f3\n");
      */

      rhox[j][i] = mid[0] + tau*dire[0];
        ux[j][i] = mid[1] + tau*dire[1];
        vx[j][i] = mid[2] + tau*dire[2];
        px[j][i] = mid[3] + tau*dire[3];
    }
}
