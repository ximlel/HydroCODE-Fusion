#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/var_struc.h"
#include "../include/Riemann_solver.h"
#include "../include/tools.h"

int Roe_GRP_solver_BN_1D
(double *Dt_U_all, double *U_all,
 const double   z_g_L, const double   z_g_R, const double   d_z_g_L, const double   d_z_g_R, const double   t_z_g_L, const double   t_z_g_R,
 const double rho_g_L, const double rho_g_R, const double d_rho_g_L, const double d_rho_g_R, const double t_rho_g_L, const double t_rho_g_R,
 const double   u_g_L, const double   u_g_R, const double   d_u_g_L, const double   d_u_g_R, const double   t_u_g_L, const double   t_u_g_R,
 const double   v_g_L, const double   v_g_R, const double   d_v_g_L, const double   d_v_g_R, const double   t_v_g_L, const double   t_v_g_R,
 const double   p_g_L, const double   p_g_R, const double   d_p_g_L, const double   d_p_g_R, const double   t_p_g_L, const double   t_p_g_R,
 const double rho_l_L, const double rho_l_R, const double d_rho_l_L, const double d_rho_l_R, const double t_rho_l_L, const double t_rho_l_R,
 const double   u_l_L, const double   u_l_R, const double   d_u_l_L, const double   d_u_l_R, const double   t_u_l_L, const double   t_u_l_R,
 const double   v_l_L, const double   v_l_R, const double   d_v_l_L, const double   d_v_l_R, const double   t_v_l_L, const double   t_v_l_R,
 const double   p_l_L, const double   p_l_R, const double   d_p_l_L, const double   d_p_l_R, const double   t_p_l_L, const double   t_p_l_R,
 const double gamma_g, const double gamma_l, const double  eps)
{
  double H_g_L, H_g_R, H_l_L, H_l_R;
  H_g_L = gamma_g/(gamma_g-1.0)*p_g_L/rho_g_L + 0.5*(u_g_L*u_g_L);
  H_g_R = gamma_g/(gamma_g-1.0)*p_g_R/rho_g_R + 0.5*(u_g_R*u_g_R);
  H_l_L = gamma_l/(gamma_l-1.0)*p_l_L/rho_l_L + 0.5*(u_l_L*u_l_L);
  H_l_R = gamma_l/(gamma_l-1.0)*p_l_R/rho_l_R + 0.5*(u_l_R*u_l_R);
  double z_l_L, z_l_R, d_z_l_L, d_z_l_R, t_z_l_L, t_z_l_R;
  z_l_L = 1.0-z_g_L;
  z_l_R = 1.0-z_g_R;
  d_z_l_L = -d_z_g_L;
  d_z_l_R = -d_z_g_R;
  t_z_l_L = -t_z_g_L;
  t_z_l_R = -t_z_g_R;
  double zrho_g_L, zrho_g_R, zrho_l_L, zrho_l_R;
  zrho_g_L = z_g_L*rho_g_L;
  zrho_g_R = z_g_R*rho_g_R;
  zrho_l_L = z_l_L*rho_l_L;
  zrho_l_R = z_l_R*rho_l_R;
  double zrho_g_S, u_g_S, H_g_S, c_g_S, zrho_l_S, u_l_S, H_l_S, c_l_S;
	zrho_g_S = sqrt(zrho_g_L*zrho_g_R);
	u_g_S = (u_g_L*sqrt(zrho_g_L)+u_g_R*sqrt(zrho_g_R)) / (sqrt(zrho_g_L)+sqrt(zrho_g_R));
	H_g_S = (H_g_L*sqrt(zrho_g_L)+H_g_R*sqrt(zrho_g_R)) / (sqrt(zrho_g_L)+sqrt(zrho_g_R));
	c_g_S = sqrt((gamma_g-1.0)*(H_g_S-0.5*u_g_S*u_g_S));
  	zrho_l_S = sqrt(zrho_l_L*zrho_l_R);
	u_l_S = (u_l_L*sqrt(zrho_l_L)+u_l_R*sqrt(zrho_l_R)) / (sqrt(zrho_l_L)+sqrt(zrho_l_R));
	H_l_S = (H_l_L*sqrt(zrho_l_L)+H_l_R*sqrt(zrho_l_R)) / (sqrt(zrho_l_L)+sqrt(zrho_l_R));
	c_l_S = sqrt((gamma_l-1.0)*(H_l_S-0.5*u_l_S*u_l_S));
  double U_INT, P_INT, b1, b2;
  U_INT = (zrho_g_S*u_g_S+zrho_l_S*u_l_S)/(zrho_g_S+zrho_l_S);
  P_INT = 0.5*(z_g_L*p_g_L+z_g_R*p_g_R)+0.5*(z_l_L*p_l_L+z_l_R*p_l_R);
//  P_INT = (gamma_g-1.0)*(H_g_S-0.5*u_g_S*u_g_S)/gamma_g*zrho_g_S + (gamma_l-1.0)*(H_l_S-0.5*u_l_S*u_l_S)/gamma_l*zrho_l_S;
  b1 =  gamma_g*P_INT/(c_g_S*c_g_S-(u_g_S-U_INT)*(u_g_S-U_INT));
  b2 = -gamma_l*P_INT/(c_l_S*c_l_S-(u_l_S-U_INT)*(u_l_S-U_INT));
  double R[7][7], R_inv[7][7];
  R[0][0] = 1.0;
  R[0][1] = 0.0;
  R[0][2] = 0.0;
  R[0][3] = 0.0;
  R[0][4] = 0.0;
  R[0][5] = 0.0;
  R[0][6] = 0.0;
  R[1][0] = b1;
  R[1][1] = 1.0;
  R[1][2] = 1.0;
  R[1][3] = 1.0;
  R[1][4] = 0.0;
  R[1][5] = 0.0;
  R[1][6] = 0.0;
  R[2][0] = U_INT*b1;
  R[2][1] = u_g_S - c_g_S;
  R[2][2] = u_g_S;
  R[2][3] = u_g_S + c_g_S;
  R[2][4] = 0.0;
  R[2][5] = 0.0;
  R[2][6] = 0.0;
  R[3][0] = (H_g_S-u_g_S*(u_g_S-U_INT))*b1-P_INT;
  R[3][1] = H_g_S - u_g_S*c_g_S;
  R[3][2] = 0.5*(u_g_S*u_g_S);
  R[3][3] = H_g_S + u_g_S*c_g_S;
  R[3][4] = 0.0;
  R[3][5] = 0.0;
  R[3][6] = 0.0;
  R[4][0] = b2;
  R[4][1] = 0.0;
  R[4][2] = 0.0;
  R[4][3] = 0.0;
  R[4][4] = 1.0;
  R[4][5] = 1.0;
  R[4][6] = 1.0;
  R[5][0] = U_INT*b2;
  R[5][1] = 0.0;
  R[5][2] = 0.0;
  R[5][3] = 0.0;
  R[5][4] = u_l_S - c_l_S;
  R[5][5] = u_l_S;
  R[5][6] = u_l_S + c_l_S;
  R[6][0] = (H_l_S-u_l_S*(u_l_S-U_INT))*b2+P_INT;
  R[6][1] = 0.0;
  R[6][2] = 0.0;
  R[6][3] = 0.0;
  R[6][4] = H_l_S - u_l_S*c_l_S;
  R[6][5] = 0.5*(u_l_S*u_l_S);
  R[6][6] = H_l_S + u_l_S*c_l_S;
	double lambda[7], alpha[7];
  	lambda[0] = U_INT;
  	lambda[1] = u_g_S - c_g_S;
	lambda[2] = u_g_S;
	lambda[3] = u_g_S + c_g_S;
  	lambda[4] = u_l_S - c_l_S;
	lambda[5] = u_l_S;
	lambda[6] = u_l_S + c_l_S;
  alpha[0] = z_g_R - z_g_L;
  alpha[1] = 0.5*((z_g_R*p_g_R-z_g_L*p_g_L)-zrho_g_S*c_g_S*(u_g_R-u_g_L)-(z_g_R-z_g_L)*(P_INT+gamma_g*P_INT/(c_g_S-(u_g_S-U_INT))*(u_g_S-U_INT)))/(c_g_S*c_g_S);
  alpha[2] = (zrho_g_R-zrho_g_L)-((z_g_R*p_g_R-z_g_L*p_g_L)+(gamma_g-1.0)*P_INT*(z_g_R-z_g_L))/(c_g_S*c_g_S);
  alpha[3] = 0.5*((z_g_R*p_g_R-z_g_L*p_g_L)+zrho_g_S*c_g_S*(u_g_R-u_g_L)-(z_g_R-z_g_L)*(P_INT-gamma_g*P_INT/(c_g_S+(u_g_S-U_INT))*(u_g_S-U_INT)))/(c_g_S*c_g_S);
  alpha[4] = 0.5*((z_l_R*p_l_R-z_l_L*p_l_L)-zrho_l_S*c_l_S*(u_l_R-u_l_L)-(z_l_R-z_l_L)*(P_INT+gamma_l*P_INT/(c_l_S-(u_l_S-U_INT))*(u_l_S-U_INT)))/(c_l_S*c_l_S);
  alpha[5] = (zrho_l_R-zrho_l_L)-((z_l_R*p_l_R-z_l_L*p_l_L)+(gamma_l-1.0)*P_INT*(z_l_R-z_l_L))/(c_l_S*c_l_S);
  alpha[6] = 0.5*((z_l_R*p_l_R-z_l_L*p_l_L)+zrho_l_S*c_l_S*(u_l_R-u_l_L)-(z_l_R-z_l_L)*(P_INT-gamma_l*P_INT/(c_l_S+(u_l_S-U_INT))*(u_l_S-U_INT)))/(c_l_S*c_l_S);
  double U[7];
	U[0] = z_g_L;
	U[1] = z_g_L*rho_g_L;
	U[2] = z_g_L*rho_g_L*u_g_L;
  	U[3] = z_g_L*(p_g_L/(gamma_g-1.0)+0.5*rho_g_L*u_g_L*u_g_L);
  	U[4] = z_l_L*rho_l_L;
	U[5] = z_l_L*rho_l_L*u_l_L;
  	U[6] = z_l_L*(p_l_L/(gamma_l-1.0)+0.5*rho_l_L*u_l_L*u_l_L);
  for(int j = 0; j < 7; j++)
      {
      for(int i = 0; i < 7; i++)
        {
          R_inv[i][j] = R[i][j];
          if(lambda[j]<=0)
              U[i] += alpha[j]*R[i][j];
        }
      }
 if (rinv(R_inv[0],7)==0)
	return 1;
 double v_g_S, v_l_S;
 if (u_g_S > 0.0)
  v_g_S = v_g_L;
 else
  v_g_S = v_g_R;
 if (u_l_S > 0.0)
   v_l_S = v_l_L;
 else
   v_l_S = v_l_R;

 U_all[0] = U[0];
 U_all[1] = U[1];
 U_all[2] = U[2];
 U_all[3] = U[1]*v_g_S;
 U_all[4] = U[3]+0.5*U[1]*v_g_S*v_g_S;
 U_all[5] = U[4];
 U_all[6] = U[5];
 U_all[7] = U[4]*v_l_S;
 U_all[8] = U[6]+0.5*U[4]*v_l_S*v_l_S;
/*
if (U_all[1]/U_all[0]>1.1)
{
	printf("%.8lf\n",U_all[1]);
}
*/

//derivative
 double D_U_L[7], D_U_R[7], T_U_L[7], T_U_R[7];
 D_U_L[0] = d_z_g_L;
 D_U_L[1] = z_g_L*d_rho_g_L+d_z_g_L*rho_g_L;
 D_U_L[2] = D_U_L[1]*u_g_L+z_g_L*rho_g_L*d_u_g_L;
 D_U_L[3] = (z_g_L*d_p_g_L+d_z_g_L*p_g_L)/(gamma_g-1)+0.5*D_U_L[2]*u_g_L+0.5*d_u_g_L*z_g_L*rho_g_L*u_g_L;
 D_U_L[4] = z_l_L*d_rho_l_L+d_z_l_L*rho_l_L;
 D_U_L[5] = D_U_L[4]*u_l_L+z_l_L*rho_l_L*d_u_l_L;
 D_U_L[6] = (z_l_L*d_p_l_L+d_z_l_L*p_l_L)/(gamma_l-1)+0.5*D_U_L[5]*u_l_L+0.5*d_u_l_L*z_l_L*rho_l_L*u_l_L;
 D_U_R[0] = d_z_g_R;
 D_U_R[1] = z_g_R*d_rho_g_R+d_z_g_R*rho_g_R;
 D_U_R[2] = D_U_R[1]*u_g_R+z_g_R*rho_g_R*d_u_g_R;
 D_U_R[3] = (z_g_R*d_p_g_R+d_z_g_R*p_g_R)/(gamma_g-1)+0.5*D_U_R[2]*u_g_R+0.5*d_u_g_R*z_g_R*rho_g_R*u_g_R;
 D_U_R[4] = z_l_R*d_rho_l_R+d_z_l_R*rho_l_R;
 D_U_R[5] = D_U_R[4]*u_l_R+z_l_R*rho_l_R*d_u_l_R;
 D_U_R[6] = (z_l_R*d_p_l_R+d_z_l_R*p_l_R)/(gamma_l-1)+0.5*D_U_R[5]*u_l_R+0.5*d_u_l_R*z_l_R*rho_l_R*u_l_R;
 T_U_L[0] = t_z_g_L;
 T_U_L[1] = z_g_L*t_rho_g_L+t_z_g_L*rho_g_L;
 T_U_L[2] = T_U_L[1]*u_g_L+z_g_L*rho_g_L*t_u_g_L;
 T_U_L[3] = (z_g_L*t_p_g_L+t_z_g_L*p_g_L)/(gamma_g-1)+0.5*T_U_L[2]*u_g_L+0.5*t_u_g_L*z_g_L*rho_g_L*u_g_L;
 T_U_L[4] = z_l_L*t_rho_l_L+t_z_l_L*rho_l_L;
 T_U_L[5] = T_U_L[4]*u_l_L+z_l_L*rho_l_L*t_u_l_L;
 T_U_L[6] = (z_l_L*t_p_l_L+t_z_l_L*p_l_L)/(gamma_l-1)+0.5*T_U_L[5]*u_l_L+0.5*t_u_l_L*z_l_L*rho_l_L*u_l_L;
 T_U_R[0] = t_z_g_R;
 T_U_R[1] = z_g_R*t_rho_g_R+t_z_g_R*rho_g_R;
 T_U_R[2] = T_U_R[1]*u_g_R+z_g_R*rho_g_R*t_u_g_R;
 T_U_R[3] = (z_g_R*t_p_g_R+t_z_g_R*p_g_R)/(gamma_g-1)+0.5*T_U_R[2]*u_g_R+0.5*t_u_g_R*z_g_R*rho_g_R*u_g_R;
 T_U_R[4] = z_l_R*t_rho_l_R+t_z_l_R*rho_l_R;
 T_U_R[5] = T_U_R[4]*u_l_R+z_l_R*rho_l_R*t_u_l_R;
 T_U_R[6] = (z_l_R*t_p_l_R+t_z_l_R*p_l_R)/(gamma_l-1)+0.5*T_U_R[5]*u_l_R+0.5*t_u_l_R*z_l_R*rho_l_R*u_l_R;

  u_g_S = U[2]/U[1];
  u_l_S = U[5]/U[4];
  U_INT = (U_all[2]+U_all[6])/(U_all[1]+U_all[5]);
  double V_INT;
  V_INT = (U_all[3]+U_all[7])/(U_all[1]+U_all[5]);
  double p_g_S, p_l_S;
  p_g_S = (U[3]-0.5*U[2]*U[2]/U[1])*(gamma_g-1.0)/U[0];
  p_l_S = (U[6]-0.5*U[5]*U[5]/U[4])*(gamma_l-1.0)/(1.0-U[0]);
  P_INT = p_g_S*U[0]+p_l_S*(1.0-U[0]);

 double D_W_L[7]={0}, D_W_R[7]={0}, T_W_L[7]={0}, T_W_R[7]={0};
 for(int i=0; i<7; i++)
  for(int j=0; j<7; j++)
  {
    D_W_L[i] += R_inv[i][j]*D_U_L[j];
    D_W_R[i] += R_inv[i][j]*D_U_R[j];
    T_W_L[i] += R_inv[i][j]*T_U_L[j];
    T_W_R[i] += R_inv[i][j]*T_U_R[j];
  }
 double D_W_S[7], T_W_S[7], D_U_S[7]={0}, T_U_S[7]={0};
 for(int i = 0; i < 7; i++)
   {
     if(lambda[i]>0.0)
        {
          D_W_S[i] = D_W_L[i];
          T_W_S[i] = T_W_L[i];
        }
     else
        {
          D_W_S[i] = D_W_R[i];
          T_W_S[i] = T_W_R[i];
        }
   }
  for(int i=0; i<7; i++)
   for(int j=0; j<7; j++)
    {
      D_U_S[i] += R[i][j]*D_W_S[j];
      T_U_S[i] += R[i][j]*T_W_S[j];
    }
 double d_v_g_S, d_v_l_S, t_v_g_S, t_v_l_S;
   if (u_g_S > 0.0)
   {
    d_v_g_S = d_v_g_L;
    t_v_g_S = t_v_g_L;
   }
   else
   {
    d_v_g_S = d_v_g_R;
    t_v_g_S = t_v_g_R;
   }
   if (u_l_S > 0.0)
   {
    d_v_l_S = d_v_l_L;
    t_v_l_S = t_v_l_L;
   }
   else
   {
    d_v_l_S = d_v_l_R;
    t_v_l_S = t_v_l_R;
   }
  double d_u_g_S, d_u_l_S, t_u_g_S, t_u_l_S;
  d_u_g_S = (D_U_S[2] - D_U_S[1]*u_g_S)/U[1];
  t_u_g_S = (T_U_S[2] - T_U_S[1]*u_g_S)/U[1];
  d_u_l_S = (D_U_S[5] - D_U_S[4]*u_l_S)/U[4];
  t_u_l_S = (T_U_S[5] - T_U_S[4]*u_l_S)/U[4];
  double d_p_g_S, d_p_l_S, t_p_g_S, t_p_l_S;
  d_p_g_S =  ((D_U_S[3] - 0.5*(D_U_S[2]*u_g_S + U[2]*d_u_g_S))*(gamma_g-1.0) - p_g_S*D_U_S[0])/U[0];
  t_p_g_S =  ((T_U_S[3] - 0.5*(T_U_S[2]*u_g_S + U[2]*t_u_g_S))*(gamma_g-1.0) - p_g_S*T_U_S[0])/U[0];
  d_p_l_S =  ((D_U_S[6] - 0.5*(D_U_S[5]*u_l_S + U[5]*d_u_l_S))*(gamma_l-1.0) + p_l_S*D_U_S[0])/(1.0-U[0]);
  t_p_l_S =  ((T_U_S[6] - 0.5*(T_U_S[5]*u_l_S + U[5]*t_u_l_S))*(gamma_l-1.0) + p_l_S*T_U_S[0])/(1.0-U[0]);

   Dt_U_all[0] = -U_INT*D_U_S[0];
   Dt_U_all[1] = -D_U_S[2];
   Dt_U_all[2] = -D_U_S[2]*u_g_S-U[2]*d_u_g_S-U[0]*d_p_g_S+(P_INT-p_g_S)*D_U_S[0];
   Dt_U_all[3] = -D_U_S[2]*v_g_S-U[2]*d_v_g_S;
   Dt_U_all[4] = -D_U_S[3]*u_g_S-U[3]*d_u_g_S-U[0]*(u_g_S*d_p_g_S+d_u_g_S*p_g_S)+(P_INT*U_INT-p_g_S*u_g_S)*D_U_S[0];
   Dt_U_all[5] = -D_U_S[5];;
   Dt_U_all[6] = -D_U_S[5]*u_l_S-U[5]*d_u_l_S-(1.0-U[0])*d_p_l_S-(P_INT-p_l_S)*D_U_S[0];
   Dt_U_all[7] = -D_U_S[5]*v_l_S-U[5]*d_v_l_S;
   Dt_U_all[8] = -D_U_S[6]*u_l_S-U[6]*d_u_l_S+(1.0-U[0])*(u_l_S*d_p_l_S+d_u_l_S*p_l_S)-(P_INT*U_INT-p_l_S*u_l_S)*D_U_S[0];

  return 0;
}
