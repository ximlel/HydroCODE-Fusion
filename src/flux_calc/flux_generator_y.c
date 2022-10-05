/**
 * @file flux_generator_y.c
 * @brief This file is a function which generates Eulerian fluxes in y-direction of 
 *        2-D Euler equations solved by 2-D GRP scheme.
 */
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "../include/var_struc.h"
#include "../include/inter_process.h"
#include "../include/flux_calc.h"


/**
 * @brief This function calculate Eulerian fluxes of 2-D Euler equations in y-direction by 2-D GRP solver.
 * @details Passes variable values on both sides of the interface to the structure variables b_f_var bfv_L and bfv_R,
 *          and use function GRP_2D_scheme() to calculate fluxes.
 * @param[in] m:      Number of the x-grids: n_x.
 * @param[in] n:      Number of the y-grids: n_y.
 * @param[in] nt:     Current plot time step for computing updates of conservative variables.
 * @param[in] tau:    The length of the time step.
 * @param[in,out] CV: Structure of cell variable data.
 * @param[in] bfv_D:  Structure pointer of fluid variables at downside boundary.
 * @param[in] bfv_U:  Structure pointer of fluid variables at upper boundary.
 * @param[in] Transversal: Whether the tangential effect is considered.
 * @return    miscalculation indicator.
 *   @retval  0: Successful calculation.
 *   @retval  1: Calculation error of left/right states.
 *   @retval  2: Calculation error of interfacial fluxes.
 */
int flux_generator_y(const int m, const int n, const int nt, const double tau, struct cell_var_stru * CV,
		      struct b_f_var * bfv_D, struct b_f_var * bfv_U, const _Bool Transversal)
{
  double const h_y = config[11]; // the length of the initial y spatial grids
  struct i_f_var ifv_D = {.n_x = 1.0, .n_y = 0.0, .gamma = config[6]};
  struct i_f_var ifv_U = ifv_D;
  int i, j, data_err, data_err_retval = 0;

//===========================
#pragma omp parallel for firstprivate(ifv_U, ifv_D) collapse(2) schedule(dynamic)
  for(j = 0; j < m; ++j)
    for(i = 0; i <= n; ++i)
    {
      if(i)
      {
          ifv_D.d_rho = CV->t_rho[j][i-1];
          ifv_D.d_u   =   CV->t_u[j][i-1];
          ifv_D.d_v   =   CV->t_v[j][i-1];
          ifv_D.d_p   =   CV->t_p[j][i-1];
          ifv_D.RHO  = CV[nt].RHO[j][i-1] + 0.5*h_y*CV->t_rho[j][i-1];
          ifv_D.U    =   CV[nt].U[j][i-1] + 0.5*h_y*  CV->t_u[j][i-1];
          ifv_D.V    =   CV[nt].V[j][i-1] + 0.5*h_y*  CV->t_v[j][i-1];
          ifv_D.P    =   CV[nt].P[j][i-1] + 0.5*h_y*  CV->t_p[j][i-1];
      }
      else
      {
          ifv_D.d_rho = bfv_D[j].TRHO;
          ifv_D.d_u   = bfv_D[j].TU;
          ifv_D.d_v   = bfv_D[j].TV;
          ifv_D.d_p   = bfv_D[j].TP;
          ifv_D.RHO   = bfv_D[j].RHO + 0.5*h_y*bfv_D[j].TRHO;
          ifv_D.U     = bfv_D[j].U   + 0.5*h_y*bfv_D[j].TU;
          ifv_D.V     = bfv_D[j].V   + 0.5*h_y*bfv_D[j].TV;
          ifv_D.P     = bfv_D[j].P   + 0.5*h_y*bfv_D[j].TP;
      }
      if(i < n)
      {
          ifv_U.d_rho = CV->t_rho[j][i];
          ifv_U.d_u   =   CV->t_u[j][i];
          ifv_U.d_v   =   CV->t_v[j][i];
          ifv_U.d_p   =   CV->t_p[j][i];
          ifv_U.RHO  = CV[nt].RHO[j][i] - 0.5*h_y*CV->t_rho[j][i];
          ifv_U.U    =   CV[nt].U[j][i] - 0.5*h_y*  CV->t_u[j][i];
          ifv_U.V    =   CV[nt].V[j][i] - 0.5*h_y*  CV->t_v[j][i];
          ifv_U.P    =   CV[nt].P[j][i] - 0.5*h_y*  CV->t_p[j][i];
      }
      else
      {
          ifv_U.d_rho = bfv_U[j].TRHO;
          ifv_U.d_u   = bfv_U[j].TU;
          ifv_U.d_v   = bfv_U[j].TV;
          ifv_U.d_p   = bfv_U[j].TP;
          ifv_U.RHO   = bfv_U[j].RHO - 0.5*h_y*bfv_U[j].TRHO;
          ifv_U.U     = bfv_U[j].U   - 0.5*h_y*bfv_U[j].TU;
          ifv_U.V     = bfv_U[j].V   - 0.5*h_y*bfv_U[j].TV;
          ifv_U.P     = bfv_U[j].P   - 0.5*h_y*bfv_U[j].TP;
      }

//===========================
      if (Transversal)
	  {
	      if(i)
		  {
		      ifv_D.t_rho = -CV->s_rho[j][i-1];
		      ifv_D.t_u   = -  CV->s_u[j][i-1];
		      ifv_D.t_v   = -  CV->s_v[j][i-1];
		      ifv_D.t_p   = -  CV->s_p[j][i-1];
		  }
	      else
		  {
		      ifv_D.t_rho = -bfv_D[j].SRHO;
		      ifv_D.t_u   = -bfv_D[j].SU;
		      ifv_D.t_v   = -bfv_D[j].SV;
		      ifv_D.t_p   = -bfv_D[j].SP;
		  }
	      if(i < n)
		  {
		      ifv_U.t_rho = -CV->s_rho[j][i];
		      ifv_U.t_u   = -  CV->s_u[j][i];
		      ifv_U.t_v   = -  CV->s_v[j][i];
		      ifv_U.t_p   = -  CV->s_p[j][i];
		  }
	      else
		  {
		      ifv_U.t_rho = -bfv_U[j].SRHO;
		      ifv_U.t_u   = -bfv_U[j].SU;
		      ifv_U.t_v   = -bfv_U[j].SV;
		      ifv_U.t_p   = -bfv_U[j].SP;
		  }
	  }
      else
	  {
	      ifv_D.t_rho = -0.0;
	      ifv_D.t_u   = -0.0;
	      ifv_D.t_v   = -0.0;
	      ifv_D.t_p   = -0.0;
	      ifv_U.t_rho = -0.0;
	      ifv_U.t_u   = -0.0;
	      ifv_U.t_v   = -0.0;
	      ifv_U.t_p   = -0.0;
	  }
      if(ifvar_check(&ifv_D, &ifv_U, 2))
	  {
	      printf(" on [%d, %d, %d] (nt, x, y).\n", nt, j, i);
	      data_err_retval = 1;
	  }

//===========================

      data_err = GRP_2D_flux(&ifv_D, &ifv_U, tau);
      switch (data_err)
	  {
	  case 1:
	      printf("<0.0 error on [%d, %d, %d] (nt, x, y) - STAR_y\n", nt, j, i);
	      data_err_retval = 2;
	  case 2:
	      printf("NAN or INFinite error on [%d, %d, %d] (nt, x, y) - STAR_y\n", nt, j, i); 
	      data_err_retval = 2;
	  case 3:
	      printf("NAN or INFinite error on [%d, %d, %d] (nt, x, y) - DIRE_y\n", nt, j, i); 
	      data_err_retval = 2;
	  }

      CV->G_rho[j][i] = ifv_D.F_rho;
      CV->G_u[j][i]   = ifv_D.F_u;
      CV->G_v[j][i]   = ifv_D.F_v;
      CV->G_e[j][i]   = ifv_D.F_e;

      CV->rhoIy[j][i] = ifv_D.RHO_int;
      CV->uIy[j][i]   = ifv_D.U_int;
      CV->vIy[j][i]   = ifv_D.V_int;
      CV->pIy[j][i]   = ifv_D.P_int;
    }
  return data_err_retval;
}
