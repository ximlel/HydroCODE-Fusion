/**
 * @file flux_generator_x.c
 * @brief This file is a function which generates Eulerian fluxes in x-direction of 
 *        2-D Euler equations solved by 2-D GRP scheme.
 */
#include <stdio.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/inter_process.h"
#include "../include/flux_calc.h"


/**
 * @brief This function calculate Eulerian fluxes of 2-D Euler equations in x-direction by 2-D GRP solver.
 * @details Passes variable values on both sides of the interface to the structure variables b_f_var bfv_L and bfv_R,
 *          and use function GRP_2D_scheme() to calculate fluxes.
 * @param[in] m:      Number of the x-grids: n_x.
 * @param[in] n:      Number of the y-grids: n_y.
 * @param[in] nt:     Current plot time step for computing updates of conservative variables.
 * @param[in] tau:    The length of the time step.
 * @param[in,out] CV: Structure of cell variable data.
 * @param[in] bfv_L:  Structure pointer of fluid variables at left boundary.
 * @param[in] bfv_R:  Structure pointer of fluid variables at right boundary.
 * @param[in] Transversal: Whether the tangential effect is considered.
 * @return    miscalculation indicator.
 *   @retval  0: Successful calculation.
 *   @retval  1: Calculation error of left/right states.
 *   @retval  2: Calculation error of interfacial fluxes.
 */
int flux_generator_x(const int m, const int n, const int nt, const double tau, struct cell_var_stru * CV,
		      struct b_f_var * bfv_L, struct b_f_var * bfv_R, const _Bool Transversal)
{
  double const h_x = config[10]; // the length of the initial x spatial grids
  struct i_f_var ifv_L = {.n_x = 1.0, .n_y = 0.0, .gamma = config[6]};
  struct i_f_var ifv_R = ifv_L;
  int i, j, data_err, data_err_retval = 0;

//===========================
#ifdef _OPENMP
#pragma omp parallel for  firstprivate(ifv_L, ifv_R) collapse(2) schedule(dynamic, 8)
#elif defined _OPENACC
#pragma acc parallel loop firstprivate(ifv_L, ifv_R) collapse(2)
#endif
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
      if(ifvar_check(&ifv_L, &ifv_R, 2))
	  {
	      printf(" on [%d, %d, %d] (nt, x, y).\n", nt, j, i);
	      data_err_retval = 1;
	  }

//===========================

      data_err = GRP_2D_flux(&ifv_L, &ifv_R, tau);
      switch (data_err)
	  {
	  case 1:
	      printf("<0.0 error on [%d, %d, %d] (nt, x, y) - STAR_x\n", nt, j, i);
	      data_err_retval = 2;
	  case 2:
	      printf("NAN or INFinite error on [%d, %d, %d] (nt, x, y) - STA_x\n", nt, j, i); 
	      data_err_retval = 2;
	  case 3:
	      printf("NAN or INFinite error on [%d, %d, %d] (nt, x, y) - DIRE_x\n", nt, j, i); 
	      data_err_retval = 2;
	  }

      CV->F_rho[j][i] = ifv_L.F_rho;
      CV->F_u[j][i]   = ifv_L.F_u;
      CV->F_v[j][i]   = ifv_L.F_v;
      CV->F_e[j][i]   = ifv_L.F_e;

      CV->rhoIx[j][i] = ifv_L.RHO_int;
      CV->uIx[j][i]   = ifv_L.U_int;
      CV->vIx[j][i]   = ifv_L.V_int;
      CV->pIx[j][i]   = ifv_L.P_int;
    } // End of parallel region
  return data_err_retval;
}
