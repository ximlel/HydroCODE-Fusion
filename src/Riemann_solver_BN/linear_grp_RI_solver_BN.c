#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/tools.h"
#include "../include/var_struc.h"
#include "../include/Riemann_solver.h"

int linear_GRP_RI_solver_BN
(struct RI_var *RI, const double D_z_s, const double z_s, const double *mid_g, const double *mid_s, 
 const struct GRP_LR_var GL, const struct GRP_LR_var GR,
 const double gamma_s, const double gamma_g, const double eps, const double tau, const int x_or_y)
{
    double z_g = 1.0-z_s;
    double rho_g = mid_g[0], u_g = mid_g[1], p_g = mid_g[3];
    double rho_s = mid_s[0], u_s = mid_s[1], p_s = mid_s[3];
    double c_s, c_g;
    c_s = sqrt(gamma_s * p_s / rho_s);
    c_g = sqrt(gamma_g * p_g / rho_g);
    double Lambda_v_p[7][7]={0.0}, Lambda_v_m[7][7]={0.0};
    Lambda_v_p[0][0]=fmax(u_s,0.0);
    Lambda_v_m[0][0]=fmin(u_s,0.0);
    Lambda_v_p[1][1]=fmax(u_s-c_s,0.0);
    Lambda_v_m[1][1]=fmin(u_s-c_s,0.0);
    Lambda_v_p[2][2]=fmax(u_s,0.0);
    Lambda_v_m[2][2]=fmin(u_s,0.0);
    Lambda_v_p[3][3]=fmax(u_s+c_s,0.0);
    Lambda_v_m[3][3]=fmin(u_s+c_s,0.0);
    Lambda_v_p[4][4]=fmax(u_g-c_g,0.0);
    Lambda_v_m[4][4]=fmin(u_g-c_g,0.0);
    Lambda_v_p[5][5]=fmax(u_g,0.0);
    Lambda_v_m[5][5]=fmin(u_g,0.0);
    Lambda_v_p[6][6]=fmax(u_g+c_g,0.0);
    Lambda_v_m[6][6]=fmin(u_g+c_g,0.0);
    double D_L[7], D_R[7];
    D_L[0] = D_z_s;
    D_R[0] = D_z_s;
    switch(x_or_y) {
    case 0:
	D_L[1] = GL.rho_sx;
	D_L[2] = GL.u_sx;
	D_L[3] = GL.Px;
	D_L[4] = GL.Qx;
	D_L[5] = GL.Hx;
	D_L[6] = GL.eta_gx;
	D_R[1] = GR.rho_sx;
	D_R[2] = GR.u_sx;
	D_R[3] = GR.Px;
	D_R[4] = GR.Qx;
	D_R[5] = GR.Hx;
	D_R[6] = GR.eta_gx;
	break;
    case 1:
	D_L[1] = GL.rho_sy;
	D_L[2] = GL.u_sy;
	D_L[3] = GL.Py;
	D_L[4] = GL.Qy;
	D_L[5] = GL.Hy;
	D_L[6] = GL.eta_gy;
	D_R[1] = GR.rho_sy;
	D_R[2] = GR.u_sy;
	D_R[3] = GR.Py;
	D_R[4] = GR.Qy;
	D_R[5] = GR.Hy;
	D_R[6] = GR.eta_gy;
	break;
    }
    double GAMMA_g = gamma_g-1.0;
    double V = u_g-u_s, T_g = pow(rho_g,GAMMA_g)/GAMMA_g;
    double R[7][7]={0.0};
    R[0][0] = 1.0;
    R[1][1] = 1.0/c_s;
    R[1][2] = 1.0;
    R[1][3] = 1.0/c_s;
    R[2][1] =-1.0/rho_s;
    R[2][2] = 0.0;
    R[2][3] = 1.0/rho_s;	
    R[3][1] = z_s*c_s + 2.0*z_g*rho_g*V/rho_s;
    R[3][3] = z_s*c_s - 2.0*z_g*rho_g*V/rho_s;
    R[3][4] = u_g-c_g-u_s;
    R[3][5] =-z_g*rho_g*T_g*GAMMA_g*V*V/c_g/c_g;
    R[3][6] = u_g+c_g-u_s;
    R[4][1] = z_g*rho_g/rho_s;
    R[4][3] =-z_g*rho_g/rho_s;
    R[4][4] = 1.0;
    R[4][5] =-z_g*rho_g*T_g*GAMMA_g*V/c_g/c_g;
    R[4][6] = 1.0;
    R[5][1] = V/rho_s;
    R[5][3] =-V/rho_s;
    R[5][4] =-c_g/z_g/rho_g;
    R[5][5] = T_g;
    R[5][6] = c_g/z_g/rho_g;
    R[6][5] = 1.0;	
    double BL[7][7], BR[7][7], W_tL[7], W_tR[7], D[7];
    mat_mul(R[0],Lambda_v_p[0],BL[0],7,7,7);
    mat_mul(R[0],Lambda_v_m[0],BR[0],7,7,7);
    if (rinv(R[0],7)==0)
    {
        exit(0);
        return 1;
    }
    mat_mul(BL[0],R[0],BL[0],7,7,7);
    mat_mul(BR[0],R[0],BR[0],7,7,7);
    mat_mul(BL[0],D_L,W_tL,7,1,7);
    mat_mul(BR[0],D_R,W_tR,7,1,7);
    mat_add(W_tL,W_tR,D,7,1);
    int i;
    for (i = 0; i < 7;i++)
	D[i] = -D[i];
    RI->z_s   = z_s   + 0.5*tau*D[0];
    RI->rho_s = rho_s + 0.5*tau*D[1];
    RI->u_s   = u_s   + 0.5*tau*D[2];
    RI->P     = z_g*rho_g*pow(u_g-u_s,2)+z_g*p_g+z_s*p_s;
    RI->P     = RI->P + 0.5*tau*D[3];
    RI->Q     = z_g*rho_g*(u_g-u_s);
    RI->Q     = RI->Q + 0.5*tau*D[4];
    RI->H     = 0.5*pow(u_g-u_s,2)+gamma_g/(gamma_g-1.0)*p_g/rho_g;
    RI->H     = RI->H + 0.5*tau*D[5];
    RI->eta_g = p_g/pow(rho_g,gamma_g);
    RI->eta_g = RI->eta_g + 0.5*tau*D[6];
    
    return 0;
}
