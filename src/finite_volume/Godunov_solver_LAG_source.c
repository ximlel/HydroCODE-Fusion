/**
 * @file  Godunov_solver_LAG_source.c
 * @brief This is a Lagrangian Godunov scheme to solve 1-D Euler equations.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#include "../include/var_struc.h"
#include "../include/Riemann_solver.h"
#include "../include/inter_process.h"
#include "../include/tools.h"


/**
 * @brief This function use Godunov scheme to solve 1-D Euler
 *        equations of motion on Lagrangian coordinate.
 * @param[in]  m:        Number of the grids.
 * @param[in,out] CV:    Structural body of cell variable data.
 * @param[in,out] X[]:   Array of the coordinate data.
 * @param[out] cpu_time: Array of the CPU time recording.
 */
void Godunov_solver_LAG_source(const int m, struct cell_var_stru CV, double * X[], double * cpu_time)
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
  int    const bound = (int)(config[17]);// the boundary condition in x-direction

  _Bool find_bound = false;
  
  double u_L, p_L, rho_L;
  double u_R, p_R, rho_R;
  double c_L, c_R; // the speeds of sound
  double h_L, h_R; // length of spatial grids
  _Bool CRW[2]; // Centred Rarefaction Wave (CRW) Indicator
  double u_star, p_star; // the Riemann solutions

  double ** RHO  = CV.RHO;
  double ** U    = CV.U;
  double ** P    = CV.P;
  double ** E    = CV.E;
  double * u_mid = malloc((m + 1) * sizeof(double)); 
  double * p_mid = malloc((m + 1) * sizeof(double));
  double * MASS  = malloc(m * sizeof(double)); // Array of the mass data in computational cells.
  if(u_mid == NULL || p_mid == NULL || MASS == NULL)
      {
	  printf("NOT enough memory! Mid Variables or MASS\n");
	  goto return_NULL;
      }
  for(k = 0; k < m; ++k) // Initialize the values of mass in computational cells
      MASS[k] = h * RHO[0][k];

  double h_S_max; // h/S_max, S_max is the maximum wave speed
  double time_c = 0.0; // the current time
  double C_m = 1.01; // a multiplicative coefficient allows the time step to increase.
  int nt = 1; // the number of times storing plotting data

  struct b_f_var bfv_L = {.H = h}; // Left  boundary condition
  struct b_f_var bfv_R = {.H = h}; // Right boundary condition
  
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
		      h_L   =   X[nt-1][j] - X[nt-1][j-1];
		      rho_L = RHO[nt-1][j-1];
		      u_L   =   U[nt-1][j-1];
		      p_L   =   P[nt-1][j-1];
		  }
	      else
		  {
		      h_L   = bfv_L.H;
		      rho_L = bfv_L.RHO;
		      u_L   = bfv_L.U;
		      p_L   = bfv_L.P;
		  }
	      if(j < m)
		  {
		      h_R   =   X[nt-1][j+1] - X[nt-1][j];
		      rho_R = RHO[nt-1][j];
		      u_R   =   U[nt-1][j];
		      p_R   =   P[nt-1][j];
		  }
	      else
		  {
		      h_R   = bfv_R.H;
		      rho_R = bfv_R.RHO;
		      u_R   = bfv_R.U;
		      p_R   = bfv_R.P;
		  }

	      c_L = sqrt(gamma * p_L / rho_L);
	      c_R = sqrt(gamma * p_R / rho_R);
	      h_S_max = fmin(h_S_max, h_L/c_L);
	      h_S_max = fmin(h_S_max, h_R/c_R);
	      if ((bound == -2 || bound == -24) && j == 0) // reflective boundary conditions
		  h_S_max = fmin(h_S_max, h_L/(fabs(u_L)+c_L));
	      if (bound == -2 && j == m)
		  h_S_max = fmin(h_S_max, h_R/(fabs(u_R)+c_R));

//========================Solve Riemann Problem========================

	      Riemann_solver_exact_single(&u_star, &p_star, gamma, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, eps, 500);

	      if(p_star < eps)
		  {
		      printf("<0.0 error on [%d, %d] (t_n, x) - STAR\n", k, j);
		      time_c = t_all;
		  }
	      if(!isfinite(p_star)|| !isfinite(u_star))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) - STAR\n", k, j); 
		      time_c = t_all;
		  }

	      u_mid[j] = u_star;
	      p_mid[j] = p_star;
	  }

//====================Time step and grid movement======================
    // If no total time, use fixed tau and time step N.
    if (isfinite(t_all) || !isfinite(config[16]) || config[16] <= 0.0)
	{
	    tau = fmin(CFL * h_S_max, C_m * tau);
	    if ((time_c + tau) > (t_all - eps))
		tau = t_all - time_c;
	}
    
    for(j = 0; j <= m; ++j)
	X[nt][j] = X[nt-1][j] + tau * u_mid[j]; // motion along the contact discontinuity

//======================THE CORE ITERATION=========================(On Lagrangian Coordinate)
    for(j = 0; j < m; ++j) // forward Euler
	{ /*
	   *  j-1          j          j+1
	   * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	   *   o-----X-----o-----X-----o-----X--...
	   */
	    RHO[nt][j] = 1.0 / (1.0/RHO[nt-1][j] + tau/MASS[j]*(u_mid[j+1] - u_mid[j]));
	    U[nt][j]   = U[nt-1][j] - tau/MASS[j]*(p_mid[j+1] - p_mid[j]);
	    E[nt][j]   = E[nt-1][j] - tau/MASS[j]*(p_mid[j+1]*u_mid[j+1] - p_mid[j]*u_mid[j]);
	    P[nt][j]   = (E[nt][j] - 0.5 * U[nt][j]*U[nt][j]) * (gamma - 1.0) * RHO[nt][j];
	    if(P[nt][j] < eps || RHO[nt][j] < eps)
		{
		    printf("<0.0 error on [%d, %d] (t_n, x) - Update\n", k, j);
		    time_c = t_all;
		}
	    if(!isfinite(P[nt][j])|| !isfinite(U[nt][j])|| !isfinite(RHO[nt][j]))
		{
		    printf("NAN or INFinite error on [%d, %d] (t_n, x) - Update\n", k, j); 
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
    for(j = 0; j <= m; ++j)
	X[nt-1][j] = X[nt][j];
    for(j = 0; j < m; ++j)
	{
	    RHO[nt-1][j] = RHO[nt][j];
	    U[nt-1][j]   =   U[nt][j];
	    E[nt-1][j]   =   E[nt][j];  
	    P[nt-1][j]   =   P[nt][j];
	}
  }

  printf("\nTime is up at time step %d.\n", k);
  printf("The cost of CPU time for 1D-Godunov Lagrangian scheme for this problem is %g seconds.\n", cpu_time_sum);
//---------------------END OF THE MAIN LOOP----------------------

return_NULL:
  free(u_mid);
  free(p_mid);
  u_mid = NULL;
  p_mid = NULL;
  free(MASS);
  MASS = NULL;
}