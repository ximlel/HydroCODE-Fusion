/**
 * @file  godunov_solver_LAG_source.c
 * @brief This is a Lagrangian Godunov scheme to solve 1-D Euler equations.
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
 *        equations of motion on Lagrangian coordinate.
 * @param[in] m:          Number of the grids.
 * @param[in,out] CV:     Structure of cell variable data.
 * @param[in,out] X[]:    Array of the coordinate data.
 * @param[out] cpu_time:  Array of the CPU time recording.
 * @param[in,out] N_plot: Pointer to the number of time steps for plotting.
 * @param[in,out] time_plot: Array of the plotting time recording.
 */
void Godunov_solver_LAG_source(const int m, struct cell_var_stru CV, double * X[], double * cpu_time, int * N_plot, double time_plot[])
{
    /* 
     * j is a frequently used index for spatial variables.
     * k is a frequently used index for the time step.
     */
  int j, k = 0;

  clock_t tic, toc;
  double cpu_time_sum = 0.0;

  double const t_all = config[1];       // the total time
  double const eps   = config[4];       // the largest value could be seen as zero
  int    const N     = (int)config[5];  // the maximum number of time steps
  double const gamma = config[6];       // the constant of the perfect gas
  double const CFL   = config[7];       // the CFL number
  double const h     = config[10];      // the length of the initial spatial grids
  double       tau   = config[16];      // the length of the time step
  int    const bound = (int)config[17]; // the boundary condition in x-direction

  _Bool find_bound = false;

  double c_L, c_R; // the speeds of sound
  double h_L, h_R; // length of spatial grids
  _Bool CRW[2]; // Centred Rarefaction Wave (CRW) Indicator
  double u_star, p_star; // the Riemann solutions

  double h_S_max; // h/S_max, S_max is the maximum wave speed
  double time_c = 0.0; // the current time
  double C_m = 1.01; // a multiplicative coefficient allows the time step to increase.
  _Bool stop_t = false;
  int nt = 0; // the number of times storing plotting data

  struct b_f_var bfv_L = {.H = h}, bfv_R = bfv_L; // Left/Right boundary condition
  struct i_f_var ifv_L = {.gamma = gamma}, ifv_R = ifv_L;

  double ** RHO  = CV.RHO;
  double ** U    = CV.U;
  double ** P    = CV.P;
  double ** E    = CV.E;
  double * U_F  = (double*)malloc((m+1) * sizeof(double));
  double * P_F  = (double*)malloc((m+1) * sizeof(double));
  double * MASS = (double*)malloc(m * sizeof(double)); // Array of the mass data in computational cells.
  if(U_F == NULL || P_F == NULL || MASS == NULL)
      {
	  printf("NOT enough memory! Variables_F or MASS\n");
	  goto return_NULL;
      }
  for(k = 0; k < m; ++k) // Initialize the values of mass in computational cells
      MASS[k] = h * RHO[0][k];
  
//-----------------------THE MAIN LOOP--------------------------------
  for(k = 1; k <= N; ++k)
  {
      tic = clock();
      if (time_c >= time_plot[nt] && nt < (*N_plot-1))
	  {
	      for(j = 0; j < m; ++j)
		  {
		      RHO[nt+1][j] = RHO[nt][j];
		      U[nt+1][j]   =   U[nt][j];
		      E[nt+1][j]   =   E[nt][j];  
		      P[nt+1][j]   =   P[nt][j];
		      X[nt+1][j]   =   X[nt][j];
		  }
	      X[nt+1][m] = X[nt][m];
	      nt++;
	  }

      h_S_max = INFINITY; // h/S_max = INFINITY

      find_bound = bound_cond_slope_limiter(true, m, nt, &CV, &bfv_L, &bfv_R, find_bound, false, time_c, X[nt]);
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
		      h_L       =   X[nt][j] - X[nt][j-1];
		      ifv_L.RHO = RHO[nt][j-1];
		      ifv_L.U   =   U[nt][j-1];
		      ifv_L.P   =   P[nt][j-1];
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
		      h_R       =   X[nt][j+1] - X[nt][j];
		      ifv_R.RHO = RHO[nt][j];
		      ifv_R.U   =   U[nt][j];
		      ifv_R.P   =   P[nt][j];
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
	      h_S_max = fmin(h_S_max, h_L/c_L);
	      h_S_max = fmin(h_S_max, h_R/c_R);
	      if ((bound == -2 || bound == -24) && j == 0) // reflective boundary conditions
		  h_S_max = fmin(h_S_max, h_L/(fabs(ifv_L.U)+c_L));
	      if (bound == -2 && j == m)
		  h_S_max = fmin(h_S_max, h_R/(fabs(ifv_R.U)+c_R));

//========================Solve Riemann Problem========================

	      Riemann_solver_exact_single(&u_star, &p_star, gamma, ifv_L.U, ifv_R.U, ifv_L.P, ifv_R.P, c_L, c_R, CRW, eps, eps, 500);

	      if(p_star < eps)
		  {
		      printf("<0.0 error on [%d, %d] (t_n, x) - STAR\n", k, j);
		      stop_t = true;
		  }
	      if(!isfinite(p_star)|| !isfinite(u_star))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) - STAR\n", k, j); 
		      stop_t = true;
		  }

	      U_F[j] = u_star;
	      P_F[j] = p_star;
	  }

//====================Time step and grid movement======================
    // If no total time, use fixed tau and time step N.
    if(isfinite(t_all) || !isfinite(config[16]) || config[16] <= 0.0)
	{
	    tau = fmin(CFL * h_S_max, C_m * tau);
	    if(tau < eps)
		{
		    printf("\nThe length of the time step is so small on [%d, %g, %g] (t_n, time_c, tau)\n", k, time_c, tau);
		    stop_t = true;
		}
	    else if((time_c + tau) > (t_all - eps))
		tau = t_all - time_c;
	    else if(!isfinite(tau))
		{
		    printf("NAN or INFinite error on [%d, %g, %g] (t_n, time_c, tau) - CFL\n", k, time_c, tau); 
		    goto return_NULL;
		}
	}
    
    for(j = 0; j <= m; ++j)
	X[nt][j] += tau * U_F[j]; // motion along the contact discontinuity

//======================THE CORE ITERATION=========================(On Lagrangian Coordinate)
    for(j = 0; j < m; ++j) // forward Euler
	{ /*
	   *  j-1          j          j+1
	   * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	   *   o-----X-----o-----X-----o-----X--...
	   */
	    RHO[nt][j] = 1.0 / (1.0/RHO[nt][j] + tau/MASS[j]*(U_F[j+1] - U_F[j]));
	    U[nt][j]   = U[nt][j] - tau/MASS[j]*(P_F[j+1] - P_F[j]);
	    E[nt][j]   = E[nt][j] - tau/MASS[j]*(P_F[j+1]*U_F[j+1] - P_F[j]*U_F[j]);
	    P[nt][j]   = (E[nt][j] - 0.5 * U[nt][j]*U[nt][j]) * (gamma - 1.0) * RHO[nt][j];
	    if(P[nt][j] < eps || RHO[nt][j] < eps)
		{
		    printf("<0.0 error on [%d, %d] (t_n, x) - Update\n", k, j);
		    stop_t = true;
		}
	}

//============================Time update=======================

    time_c += tau;
    if(isfinite(t_all))
        DispPro(time_c*100.0/t_all, k);
    else
        DispPro(k*100.0/N, k);
    if(stop_t || time_c > (t_all - eps) || !isfinite(time_c))
	break;

//===========================Fixed variable location=======================

    toc = clock();
    cpu_time_sum += ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;
    cpu_time[nt]  = cpu_time_sum;
  }

  printf("\nTime is up at time step %d.\n", k);
  printf("The cost of CPU time for 1D-Godunov Lagrangian scheme for this problem is %g seconds.\n", cpu_time_sum);
//---------------------END OF THE MAIN LOOP----------------------

return_NULL:
  config[5] = (double)k;
  *N_plot = nt+1;
  if(isfinite(time_c))
      time_plot[nt] = time_c;
  else if(isfinite(t_all))
      time_plot[nt] = t_all;
  else if(isfinite(tau))
      time_plot[nt] = k*tau;

  free(U_F);
  free(P_F);
  U_F = NULL;
  P_F = NULL;
  free(MASS);
  MASS = NULL;
}
