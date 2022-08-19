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
#include "../include/tools.h"


/**
 * @brief This function use Godunov scheme to solve 1-D Euler
 *        equations of motion on Lagrangian coordinate.
 * @param[in]  m:        Number of the grids.
 * @param[in,out] CV:    Structural body of cell variable data.
 * @param[in,out] X[]:   Array of the coordinate data.
 * @param[out] cpu_time: Array of the CPU time recording.
 */
void Godunov_solver_LAG_source
(const int m, struct cell_var CV, double * X[], double * cpu_time)
{
    double ** RHO = CV.RHO;
    double ** U   = CV.U;
    double ** P   = CV.P;
    double ** E   = CV.E;
    /* 
     * j is a frequently used index for spatial variables.
     * k is a frequently used index for the time step.
     */
  int j, k;

  clock_t tic, toc;
  double cpu_time_sum = 0.0;

  double const t_all = config[1]<0.0 ? INFINITY : config[1]; // the total time
  double const eps   = config[4];        // the largest value could be seen as zero
  int    const N     = (int)(config[5]); // the maximum number of time steps
  double const gamma = config[6];        // the constant of the perfect gas
  double const CFL   = config[7];        // the CFL number
  double const h     = config[10];       // the length of the initial spatial grids
  double       tau   = config[16];       // the length of the time step
  int    const bound = (int)(config[17]);// the boundary condition

  _Bool find_bound = false;
  
  double u_L, p_L, rho_L;
  double u_R, p_R, rho_R;
  double c_L, c_R; // the speeds of sound
  double h_L, h_R; // length of spatial grids
  _Bool CRW[2]; // Centred Rarefaction Wave (CRW) Indicator
  double u_star, p_star; // the Riemann solutions
  double *u_mid = malloc((m + 1) * sizeof(double)); 
  double *p_mid = malloc((m + 1) * sizeof(double));
  double * MASS; // Array of the mass data in computational cells.
  if(u_mid == NULL || p_mid == NULL)
      {
	  printf("NOT enough memory! Mid Variables\n");
	  goto return_NULL;
      }
  MASS = malloc(m * sizeof(double));
  if(MASS == NULL)
      {
	  printf("NOT enough memory! MASS\n");
	  goto return_NULL;
      }
  for(k = 0; k < m; ++k) // Initialize the values of mass in computational cells
      MASS[k] = h * RHO[0][k];

  double h_S_max; // h/S_max, S_max is the maximum wave speed
  double time_c = 0.0; // the current time
  double C_m = 1.01; // a multiplicative coefficient allows the time step to increase.
  int n = 1; // the number of times storing plotting data

  double UL, PL, RHOL, CL, HL; // Left  boundary condition
  double UR, PR, RHOR, CR, HR; // Right boundary condition

//-----------------------THE MAIN LOOP--------------------------------
  for(k = 1; k <= N; ++k)
  {
      h_S_max = INFINITY; // h/S_max = INFINITY
      tic = clock();

      switch (bound)
	  {
	  case -1: // initial boudary conditions
	      if(find_bound)
		  break;
	      else
		  printf("Initial boudary conditions.\n");		  
	      find_bound = true;
	      UL   =   U[0][0]; UR   =   U[0][m-1];
	      PL   =   P[0][0]; PR   =   P[0][m-1];
	      RHOL = RHO[0][0]; RHOR = RHO[0][m-1];
	      HL  = h; HR = h;
	      break;
	  case -2: // reflective boundary conditions
	      if(!find_bound)
		  printf("Reflective boudary conditions.\n");
	      find_bound = true;
	      UL   = - U[n-1][0]; UR   = - U[n-1][m-1];
	      PL   =   P[n-1][0]; PR   =   P[n-1][m-1];
	      RHOL = RHO[n-1][0]; RHOR = RHO[n-1][m-1];
	      HL = X[n-1][1] - X[n-1][0];
	      HR = X[n-1][m] - X[n-1][m-1];
	      CL = sqrt(gamma * PL / RHOL);
	      CR = sqrt(gamma * PR / RHOR);
	      h_S_max = fmin(HL/(fabs(UL)+CL), HR/(fabs(UR)+CR));
	      break;
	  case -4: // free boundary conditions
	      if(!find_bound)
		  printf("Free boudary conditions.\n");
	      find_bound = true;
	      UL   =   U[n-1][0]; UR   =   U[n-1][m-1];
	      PL   =   P[n-1][0]; PR   =   P[n-1][m-1];
	      RHOL = RHO[n-1][0]; RHOR = RHO[n-1][m-1];
	      HL = X[n-1][1] - X[n-1][0];
	      HR = X[n-1][m] - X[n-1][m-1];
	      break;
	  case -5: // periodic boundary conditions
	      if(!find_bound)
		  printf("Periodic boudary conditions.\n");
	      find_bound = true;
	      UL   =   U[n-1][m-1]; UR   =   U[n-1][0];
	      PL   =   P[n-1][m-1]; PR   =   P[n-1][0];
	      RHOL = RHO[n-1][m-1]; RHOR = RHO[n-1][0];
	      HL = X[n-1][m] - X[n-1][m-1];
	      HR = X[n-1][1] - X[n-1][0];
	      break;
	  case -24: // reflective + free boundary conditions
	      if(!find_bound)
		  printf("Reflective + Free boudary conditions.\n");
	      find_bound = true;
	      UL   = - U[n-1][0]; UR   =   U[n-1][m-1];
	      PL   =   P[n-1][0]; PR   =   P[n-1][m-1];
	      RHOL = RHO[n-1][0]; RHOR = RHO[n-1][m-1];
	      HL = X[n-1][1] - X[n-1][0];
	      HR = X[n-1][m] - X[n-1][m-1];
	      CL = sqrt(gamma * PL / RHOL);
	      h_S_max = HL/(fabs(UL)+CL);
	      break;
	  default:
	      printf("No suitable boundary coditions!\n");
	      goto return_NULL;
	  }
      
      for(j = 0; j <= m; ++j)
	  { /*
	     *  j-1          j          j+1
	     * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	     *   o-----X-----o-----X-----o-----X--...
	     */
	      if(j) // Initialize the initial values.
		  {
		      h_L   =   X[n-1][j] - X[n-1][j-1];
		      rho_L = RHO[n-1][j-1];
		      u_L   =   U[n-1][j-1];
		      p_L   =   P[n-1][j-1];
		  }
	      else
		  {
		      h_L   =   HL;
		      rho_L = RHOL;
		      u_L   =   UL;
		      p_L   =   PL;
		  }
	      if(j < m)
		  {
		      h_R   =   X[n-1][j+1] - X[n-1][j];
		      rho_R = RHO[n-1][j];
		      u_R   =   U[n-1][j];
		      p_R   =   P[n-1][j];
		  }
	      else
		  {
		      h_R   =   HR;
		      rho_R = RHOR;
		      u_R   =   UR;
		      p_R   =   PR;
		  }

	      c_L = sqrt(gamma * p_L / rho_L);
	      c_R = sqrt(gamma * p_R / rho_R);
	      h_S_max = fmin(h_S_max, h_L/c_L);
	      h_S_max = fmin(h_S_max, h_R/c_R);

//========================Solve Riemann Problem========================

	      Riemann_solver_exact(&u_star, &p_star, gamma, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, eps, 500);

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
	X[n][j] = X[n-1][j] + tau * u_mid[j]; // motion along the contact discontinuity

//======================THE CORE ITERATION=========================(On Lagrangian Coordinate)
    for(j = 0; j < m; ++j) // forward Euler
	{ /*
	   *  j-1          j          j+1
	   * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	   *   o-----X-----o-----X-----o-----X--...
	   */
	    RHO[n][j] = 1.0 / (1.0/RHO[n-1][j] + tau/MASS[j]*(u_mid[j+1] - u_mid[j]));
	    U[n][j]   = U[n-1][j] - tau/MASS[j]*(p_mid[j+1] - p_mid[j]);
	    E[n][j]   = E[n-1][j] - tau/MASS[j]*(p_mid[j+1]*u_mid[j+1] - p_mid[j]*u_mid[j]);
	    P[n][j]   = (E[n][j] - 0.5 * U[n][j]*U[n][j]) * (gamma - 1.0) * RHO[n][j];
	    if(P[n][j] < eps || RHO[n][j] < eps)
		{
		    printf("<0.0 error on [%d, %d] (t_n, x) - Update\n", k, j);
		    time_c = t_all;
		}
	    if(!isfinite(P[n][j])|| !isfinite(U[n][j])|| !isfinite(RHO[n][j]))
		{
		    printf("NAN or INFinite error on [%d, %d] (t_n, x) - Update\n", k, j); 
		    time_c = t_all;
		}
	}

//============================Time update=======================

    toc = clock();
    cpu_time[n] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    cpu_time_sum += cpu_time[n];

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
	X[n-1][j] = X[n][j];
    for(j = 0; j < m; ++j)
	{
	    RHO[n-1][j] = RHO[n][j];
	    U[n-1][j]   =   U[n][j];
	    E[n-1][j]   =   E[n][j];  
	    P[n-1][j]   =   P[n][j];
	}
  }

  printf("\nTime is up in time step %d.\n", k);
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
