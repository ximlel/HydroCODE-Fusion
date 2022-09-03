/**
 * @file  godunov_solver_ALE_source.c
 * @brief This is an ALE Godunov scheme to solve 1-D Euler equations.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#include "../include/var_struc.h"
#include "../include/riemann_solver.h"
#include "../include/inter_process.h"
#include "../include/tools.h"


/**
 * @brief This function use Godunov scheme to solve 1-D Euler
 *        equations of motion on ALE coordinate.
 * @param[in] m:         Number of the grids.
 * @param[in,out] CV:     Structure of cell variable data.
 * @param[in,out] X[]:    Array of the coordinate data.
 * @param[out] cpu_time:  Array of the CPU time recording.
 * @param[out] time_plot: Array of the plotting time recording.
 * @todo All of the functionality of the ALE code has not yet been implemented.
 */
void Godunov_solver_ALE_source_Undone(const int m, struct cell_var_stru CV, double * X[], double * cpu_time, double * time_plot)
{
    /* 
     * j is a frequently used index for spatial variables.
     * k is a frequently used index for the time step.
     */
  int j, k;

  clock_t tic, toc;
  double cpu_time_sum = 0.0;

  double const t_all = config[1];        // the total time
  double const eps   = config[4];        // the largest value could be seen as zero
  int    const N     = (int)(config[5]); // the maximum number of time steps
  double const gamma = config[6];        // the constant of the perfect gas
  double const CFL   = config[7];        // the CFL number
  double const h     = config[10];       // the length of the initial spatial grids
  double       tau   = config[16];       // the length of the time step

  _Bool find_bound = false;

  double Mom, Ene;
  double c_L, c_R; // the speeds of sound
  double h_L, h_R; // length of spatial grids
  /*
   * mid:  the Riemann solutions.
   *       [rho_star, u_star, p_star]
   */
  double dire[3], mid[3];

  double ** RHO  = CV.RHO;
  double ** U    = CV.U;
  double ** P    = CV.P;
  double ** E    = CV.E;
  // the numerical flux at (x_{j-1/2}, t_{n}).
  double * F_rho = malloc((m+1) * sizeof(double));
  double * F_u   = malloc((m+1) * sizeof(double));
  double * F_e   = malloc((m+1) * sizeof(double));
  if(F_rho == NULL || F_u == NULL || F_e == NULL)
      {
	  printf("NOT enough memory! Flux\n");
	  goto return_NULL;
      }

  double nu;  // nu = tau/h
  double h_S_max; // h/S_max, S_max is the maximum wave speed
  double time_c = 0.0; // the current time
  int nt = 1; // the number of times storing plotting data

  struct b_f_var bfv_L = {.H = h, .SU = 0.0, .SP = 0.0, .SRHO = 0.0}; // Left  boundary condition
  struct b_f_var bfv_R = {.H = h, .SU = 0.0, .SP = 0.0, .SRHO = 0.0}; // Right boundary condition
  struct i_f_var ifv_L = {.gamma = gamma}, ifv_R = {.gamma = gamma};
  
//-----------------------THE MAIN LOOP--------------------------------
  for(k = 1; k <= N; ++k)
  {
      h_S_max = INFINITY; // h/S_max = INFINITY
      tic = clock();

      find_bound = bound_cond_slope_limiter(true, m, nt-1, CV, &bfv_L, &bfv_R, find_bound, false, time_c, X[nt-1]);
      if(!find_bound)
	  goto return_NULL;

      for(j = 0; j <= m; ++j)
	  { /*
	     *  j-1          j          j+1
	     * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	     *   o-----X-----o-----X-----o-----X--...
	     */
	      if(j) // Initialize the initial values.
		  {
		      h_L       =   X[nt-1][j] - X[nt-1][j-1];
		      ifv_L.RHO = RHO[nt-1][j-1];
		      ifv_L.U   =   U[nt-1][j-1];
		      ifv_L.P   =   P[nt-1][j-1];
		  }
	      else
		  {
		      h_L       = bfv_L.H;
		      ifv_L.RHO = bfv_L.RHO;
		      ifv_L.U   = bfv_L.U;
		      ifv_L.P   = bfv_L.P;
		  }
	      if(j < m)
		  {
		      h_R       =   X[nt-1][j+1] - X[nt-1][j];
		      ifv_R.RHO = RHO[nt-1][j];
		      ifv_R.U   =   U[nt-1][j];
		      ifv_R.P   =   P[nt-1][j];
		  }
	      else
		  {
		      h_R       = bfv_R.H;
		      ifv_R.RHO = bfv_R.RHO;
		      ifv_R.U   = bfv_R.U;
		      ifv_R.P   = bfv_R.P;
		  }

	      c_L = sqrt(gamma * ifv_L.P / ifv_L.RHO);
	      c_R = sqrt(gamma * ifv_R.P / ifv_R.RHO);
	      h_S_max = fmin(h_S_max, h_L/(fabs(ifv_L.U)+fabs(c_L)));
	      h_S_max = fmin(h_S_max, h_R/(fabs(ifv_R.U)+fabs(c_R)));

//========================Solve Riemann Problem========================
	      linear_GRP_solver_Edir(dire, mid, ifv_L, ifv_R, eps, INFINITY);

	      if(mid[2] < eps || mid[0] < eps)
		  {
		      printf("<0.0 error on [%d, %d] (t_n, x) - STAR\n", k, j);
		      time_c = t_all;
		  }
	      if(!isfinite(mid[1])|| !isfinite(mid[2])|| !isfinite(mid[0]))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) - STAR\n", k, j); 
		      time_c = t_all;
		  }

	      F_rho[j] = mid[0]*mid[1];
	      F_u[j] = F_rho[j]*mid[1] + mid[2];
	      F_e[j] = (gamma/(gamma-1.0))*mid[2] + 0.5*F_rho[j]*mid[1];
	      F_e[j] = F_e[j]*mid[1];
	  }

//====================Time step and grid fixed======================
    // If no total time, use fixed tau and time step N.
    if (isfinite(t_all) || !isfinite(config[16]) || config[16] <= 0.0)
	{
	    tau = CFL * h_S_max;
	    if ((time_c + tau) > (t_all - eps))
		tau = t_all - time_c;
	    else if(!isfinite(tau))
		{
		    printf("NAN or INFinite error on [%d, %g] (t_n, tau) - CFL\n", k, tau); 
		    tau = t_all - time_c;
		    goto return_NULL;
		}
	}
    nu = tau / h;

    for (j = 0; j <= m; ++j)
	X[nt][j] = X[nt-1][j];

//======================THE CORE ITERATION=========================(On Eulerian Coordinate)
    for(j = 0; j < m; ++j) // forward Euler
	{ /*
	   *  j-1          j          j+1
	   * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	   *   o-----X-----o-----X-----o-----X--...
	   */
	    RHO[nt][j] = RHO[nt-1][j]     - nu*(F_rho[j+1]-F_rho[j]);
	    Mom = RHO[nt-1][j]*U[nt-1][j] - nu*(F_u[j+1]  -F_u[j]);
	    Ene = RHO[nt-1][j]*E[nt-1][j] - nu*(F_e[j+1]  -F_e[j]);

	    U[nt][j] = Mom / RHO[nt][j];
	    E[nt][j] = Ene / RHO[nt][j];
	    P[nt][j] = (Ene - 0.5*Mom*U[nt][j])*(gamma-1.0);

	    if(P[nt][j] < eps || RHO[nt][j] < eps)
		{
		    printf("<0.0 error on [%d, %d] (t_n, x) - Update\n", k, j);
		    time_c = t_all;
		}
	}

//============================Time update=======================

    toc = clock();
    cpu_time[nt] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    cpu_time_sum += cpu_time[nt];

    time_c += tau;
    if (isfinite(t_all))
        DispPro(time_c*100.0/t_all, k);
    else
        DispPro(k*100.0/N, k);
    if(time_c > (t_all - eps) || isinf(time_c))
	{
	    config[5] = (double)k;
	    break;
	}

//===========================Fixed variable location=======================	
    for(j = 0; j < m; ++j)
	{
	    RHO[nt-1][j] = RHO[nt][j];
	    U[nt-1][j]   =   U[nt][j];
	    E[nt-1][j]   =   E[nt][j];  
	    P[nt-1][j]   =   P[nt][j];
	}
  }

  time_plot[0] = time_c - tau;
  time_plot[1] = time_c;
  printf("\nTime is up at time step %d.\n", k);
  printf("The cost of CPU time for 1D-Godunov Eulerian scheme for this problem is %g seconds.\n", cpu_time_sum);
//---------------------END OF THE MAIN LOOP----------------------

return_NULL:
  free(F_rho);
  free(F_u);
  free(F_e);
  F_rho = NULL;
  F_u   = NULL;
  F_e   = NULL;
}