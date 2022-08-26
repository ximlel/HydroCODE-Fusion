#include <stdio.h>

#include "../include/var_struc.h"
#include "../include/flux_calc.h"


void flux_generator_x(const int m, const int n, const int nt, const double tau, struct cell_var_stru * CV,
		      struct b_f_var * bfv_L, struct b_f_var * bfv_R, const _Bool Transversal)
{
  double const h_x = config[10]; // the length of the initial x spatial grids
  struct i_f_var ifv_L = {.n_x = 1.0, .n_y = 0.0}, ifv_R = {.n_x = 1.0, .n_y = 0.0};
  int i, j;

//===========================
  for(i = 0; i < n; ++i)
    for(j = 0; j <= m; ++j)
    {
      if(j)
      {
          ifv_L.d_rho = CV->s_rho[j-1][i];
          ifv_L.d_u   =   CV->s_u[j-1][i];
          ifv_L.d_v   =   CV->s_v[j-1][i];
          ifv_L.d_p   =   CV->s_p[j-1][i];
          ifv_L.RHO  = CV[nt].RHO[j-1][i] + 0.5*h_x*CV->s_rho[j-1][i];
          ifv_L.U    =   CV[nt].U[j-1][i] + 0.5*h_x*  CV->s_u[j-1][i];
          ifv_L.V    =   CV[nt].V[j-1][i] + 0.5*h_x*  CV->s_v[j-1][i];
          ifv_L.P    =   CV[nt].P[j-1][i] + 0.5*h_x*  CV->s_p[j-1][i];
      }
      else
      {
          ifv_L.d_rho = bfv_L[i].SRHO;
          ifv_L.d_u   = bfv_L[i].SU;
          ifv_L.d_v   = bfv_L[i].SV;
          ifv_L.d_p   = bfv_L[i].SP;
          ifv_L.RHO   = bfv_L[i].RHO + 0.5*h_x*bfv_L[i].SRHO;
          ifv_L.U     = bfv_L[i].U   + 0.5*h_x*bfv_L[i].SU;
          ifv_L.V     = bfv_L[i].V   + 0.5*h_x*bfv_L[i].SV;
          ifv_L.P     = bfv_L[i].P   + 0.5*h_x*bfv_L[i].SP;
      }
      if(j < m)
      {
          ifv_R.d_rho = CV->s_rho[j][i];
          ifv_R.d_u   =   CV->s_u[j][i];
          ifv_R.d_v   =   CV->s_v[j][i];
          ifv_R.d_p   =   CV->s_p[j][i];
          ifv_R.RHO  = CV[nt].RHO[j][i] - 0.5*h_x*CV->s_rho[j][i];
          ifv_R.U    =   CV[nt].U[j][i] - 0.5*h_x*  CV->s_u[j][i];
          ifv_R.V    =   CV[nt].V[j][i] - 0.5*h_x*  CV->s_v[j][i];
          ifv_R.P    =   CV[nt].P[j][i] - 0.5*h_x*  CV->s_p[j][i];
      }
      else
      {
          ifv_R.d_rho = bfv_R[i].SRHO;
          ifv_R.d_u   = bfv_R[i].SU;
          ifv_R.d_v   = bfv_R[i].SV;
          ifv_R.d_p   = bfv_R[i].SP;
          ifv_R.RHO   = bfv_R[i].RHO - 0.5*h_x*bfv_R[i].SRHO;
          ifv_R.U     = bfv_R[i].U   - 0.5*h_x*bfv_R[i].SU;
          ifv_R.V     = bfv_R[i].V   - 0.5*h_x*bfv_R[i].SV;
          ifv_R.P     = bfv_R[i].P   - 0.5*h_x*bfv_R[i].SP;
      }
//===========================

      if (Transversal)
	  {
	      if(j)
		  {
		      ifv_L.t_rho = CV->t_rho[j-1][i];
		      ifv_L.t_u   =   CV->t_u[j-1][i];
		      ifv_L.t_v   =   CV->t_v[j-1][i];
		      ifv_L.t_p   =   CV->t_p[j-1][i];
		  }
	      else
		  {
		      ifv_L.t_rho = bfv_L[i].TRHO;
		      ifv_L.t_u   = bfv_L[i].TU;
		      ifv_L.t_v   = bfv_L[i].TV;
		      ifv_L.t_p   = bfv_L[i].TP;
		  }
	      if(j < m)
		  {
		      ifv_R.t_rho = CV->t_rho[j][i];
		      ifv_R.t_u   =   CV->t_u[j][i];
		      ifv_R.t_v   =   CV->t_v[j][i];
		      ifv_R.t_p   =   CV->t_p[j][i];
		  }
	      else
		  {
		      ifv_R.t_rho = bfv_R[i].TRHO;
		      ifv_R.t_u   = bfv_R[i].TU;
		      ifv_R.t_v   = bfv_R[i].TV;
		      ifv_R.t_p   = bfv_R[i].TP;
		  }
	  }
      else
	  {
	      ifv_L.t_rho = 0.0;
	      ifv_L.t_u   = 0.0;
	      ifv_L.t_v   = 0.0;
	      ifv_L.t_p   = 0.0;
	      ifv_R.t_rho = 0.0;
	      ifv_R.t_u   = 0.0;
	      ifv_R.t_v   = 0.0;
	      ifv_R.t_p   = 0.0;
	  }
//===========================

      GRP_2D_scheme(&ifv_L, &ifv_R, tau);

      CV->F_rho[j][i] = ifv_L.F_rho;
      CV->F_u[j][i]   = ifv_L.F_u;
      CV->F_v[j][i]   = ifv_L.F_v;
      CV->F_e[j][i]   = ifv_L.F_e;

      CV->rhoIx[j][i] = ifv_L.RHO_int;
      CV->uIx[j][i]   = ifv_L.U_int;
      CV->vIx[j][i]   = ifv_L.V_int;
      CV->pIx[j][i]   = ifv_L.P_int;
    }
}
