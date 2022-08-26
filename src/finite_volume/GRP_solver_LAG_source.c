/**
 * @file  GRP_solver_LAG_source.c
 * @brief This is a Lagrangian GRP scheme to solve 1-D Euler equations.
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
 * @brief This function use GRP scheme to solve 1-D Euler
 *        equations of motion on Lagrangian coordinate.
 * @param[in]  m:        Number of the grids.
 * @param[in,out] CV:    Structural body of cell variable data.
 * @param[in,out] X[]:   Array of the coordinate data.
 * @param[out] cpu_time: Array of the CPU time recording.
 */
void GRP_solver_LAG_source(const int m, struct cell_var_stru CV, double * X[], double * cpu_time)
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

  double t_u_L, t_p_L, t_rho_L; // spatial derivatives in Lagrangian coordinate ξ
  double t_u_R, t_p_R, t_rho_R;
  /*
   * dire: the temporal derivative of fluid variables.
   *       \frac{\partial [rho_L, u, p, rho_R]}{\partial t}
   * mid:  the Riemann solutions.
   *       [rho_star_L, u_star, p_star, rho_star_R]
   */
  double dire[4], mid[4];

  double ** RHO  = CV.RHO;
  double ** U    = CV.U;
  double ** P    = CV.P;
  double ** E    = CV.E;
  // the slopes of variable values
  double * s_rho = calloc(m, sizeof(double));
  double * s_u   = calloc(m, sizeof(double));
  double * s_p   = calloc(m, sizeof(double));
  CV.d_rho = s_rho;
  CV.d_u   = s_u;
  CV.d_p   = s_p;
  // the variable values at (x_{j-1/2}, t_{n+1}).
  double * U_next     = malloc((m+1) * sizeof(double));
  double * P_next     = malloc((m+1) * sizeof(double));
  double * RHO_next_L = malloc((m+1) * sizeof(double));
  double * RHO_next_R = malloc((m+1) * sizeof(double));
  // the temporal derivatives at (x_{j-1/2}, t_{n}).
  double * U_t     = malloc((m+1) * sizeof(double));
  double * P_t     = malloc((m+1) * sizeof(double));
  double * RHO_t_L = malloc((m+1) * sizeof(double));
  double * RHO_t_R = malloc((m+1) * sizeof(double));
  // the numerical flux at (x_{j-1/2}, t_{n+1/2}).
  double * U_F  = malloc((m+1) * sizeof(double));
  double * P_F  = malloc((m+1) * sizeof(double));
  double * MASS = malloc(m * sizeof(double));; // Array of the mass data in computational cells.
  if(s_rho == NULL || s_u == NULL || s_p == NULL)
      {
	  printf("NOT enough memory! Slope\n");
	  goto return_NULL;
      }
  if(U_next == NULL || P_next == NULL || RHO_next_L == NULL || RHO_next_R == NULL)
      {
	  printf("NOT enough memory! Variables_next\n");
	  goto return_NULL;
      }
  if(U_t == NULL || P_t == NULL || RHO_t_L == NULL || RHO_t_R == NULL)
      {
	  printf("NOT enough memory! Temproal derivative\n");
	  goto return_NULL;
      }
  if(U_F == NULL || P_F == NULL || MASS == NULL)
      {
	  printf("NOT enough memory! Variables_F or MASS\n");
	  goto return_NULL;
      }
  for(k = 0; k < m; ++k) // Initialize the values of mass in computational cells
      MASS[k] = h * RHO[0][k];

  double h_S_max; // h/S_max, S_max is the maximum wave speed
  double time_c = 0.0; // the current time
  double C_m = 1.01; // a multiplicative coefficient allows the time step to increase.
  int nt = 1; // the number of times storing plotting data

  struct b_f_var bfv_L = {.H = h, .SU = 0.0, .SP = 0.0, .SRHO = 0.0}; // Left  boundary condition
  struct b_f_var bfv_R = {.H = h, .SU = 0.0, .SP = 0.0, .SRHO = 0.0}; // Right boundary condition
  
//-----------------------THE MAIN LOOP--------------------------------
  for(k = 1; k <= N; ++k)
  {
      h_S_max = INFINITY; // h/S_max = INFINITY
      tic = clock();

      find_bound = bound_cond_slope_limiter(true, m, nt-1, CV, &bfv_L, &bfv_R, find_bound, true, time_c, X[nt-1]);
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
		      rho_L = RHO[nt-1][j-1] + 0.5*h_L*s_rho[j-1];
		      u_L   =   U[nt-1][j-1] + 0.5*h_L*s_u[j-1];
		      p_L   =   P[nt-1][j-1] + 0.5*h_L*s_p[j-1];
		  }
	      else
		  {
		      h_L   = bfv_L.H;
		      rho_L = bfv_L.RHO + 0.5*h_L*bfv_L.SRHO;
		      u_L   = bfv_L.U   + 0.5*h_L*bfv_L.SU;
		      p_L   = bfv_L.P   + 0.5*h_L*bfv_L.SP;
		  }
	      if(j < m)
		  {
		      h_R   =   X[nt-1][j+1] - X[nt-1][j];
		      rho_R = RHO[nt-1][j] - 0.5*h_R*s_rho[j];
		      u_R   =   U[nt-1][j] - 0.5*h_R*s_u[j];
		      p_R   =   P[nt-1][j] - 0.5*h_R*s_p[j];
		  }
	      else
		  {
		      h_R   = bfv_R.H;
		      rho_R = bfv_R.RHO + 0.5*h_R*bfv_R.SRHO;
		      u_R   = bfv_R.U   + 0.5*h_R*bfv_R.SU;
		      p_R   = bfv_R.P   + 0.5*h_R*bfv_R.SP;
		  }
	      if(p_L < eps || p_R < eps || rho_L < eps || rho_R < eps)
		  {
		      printf("<0.0 error on [%d, %d] (t_n, x) - Reconstruction\n", k, j);
		      goto return_NULL;
		  }
	      if(!isfinite(p_L)|| !isfinite(p_R)|| !isfinite(u_L)|| !isfinite(u_R)|| !isfinite(rho_L)|| !isfinite(rho_R))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) - Reconstruction\n", k, j); 
		      goto return_NULL;
		  }

	      c_L = sqrt(gamma * p_L / rho_L);
	      c_R = sqrt(gamma * p_R / rho_R);
	      h_S_max = fmin(h_S_max, h_L/c_L);
	      h_S_max = fmin(h_S_max, h_R/c_R);
	      if ((bound == -2 || bound == -24) && j == 0) // reflective boundary conditions
		  h_S_max = fmin(h_S_max, h_L/(fabs(u_L)+c_L));
	      if (bound == -2 && j == m)
		  h_S_max = fmin(h_S_max, h_R/(fabs(u_R)+c_R));

	      if(j) //calculate the material derivatives
		  {
		      t_u_L   =   s_u[j-1]/rho_L;
		      t_p_L   =   s_p[j-1]/rho_L;
		      t_rho_L = s_rho[j-1]/rho_L;
		  }
	      else
		  {
		      t_rho_L = bfv_L.SRHO/rho_L;
		      t_u_L   = bfv_L.SU  /rho_L;
		      t_p_L   = bfv_L.SP  /rho_L;
		  }
	      if(j < m)
		  {
		      t_u_R   =   s_u[j]/rho_R;
		      t_p_R   =   s_p[j]/rho_R;
		      t_rho_R = s_rho[j]/rho_R;
		  }
	      else
		  {
		      t_rho_R = bfv_R.SRHO/rho_R;
		      t_u_R   = bfv_R.SU  /rho_R;
		      t_p_R   = bfv_R.SP  /rho_R;
		  }

//========================Solve GRP========================
	      linear_GRP_solver_LAG(dire, mid, rho_L, rho_R, t_rho_L, t_rho_R, u_L, u_R, t_u_L, t_u_R, p_L, p_R, t_p_L, t_p_R, gamma, eps, eps);

	      if(mid[2] < eps)
		  {
		      printf("<0.0 error on [%d, %d] (t_n, x) - STAR\n", k, j);
		      time_c = t_all;
		  }
	      if(!isfinite(mid[1])|| !isfinite(mid[2]))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) - STAR\n", k, j); 
		      time_c = t_all;
		  }

	      RHO_next_L[j] = mid[0];
	      RHO_next_R[j] = mid[3];
	      U_next[j]     = mid[1];
	      P_next[j]     = mid[2];
	      RHO_t_L[j] = dire[0];
	      RHO_t_R[j] = dire[3];
	      U_t[j]     = dire[1];
	      P_t[j]     = dire[2];
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
	{
	    U_F[j] = U_next[j] + 0.5 * tau * U_t[j];
	    P_F[j] = P_next[j] + 0.5 * tau * P_t[j];

	    RHO_next_L[j] += tau * RHO_t_L[j];
	    RHO_next_R[j] += tau * RHO_t_R[j];
	    U_next[j]     += tau * U_t[j];
	    P_next[j]     += tau * P_t[j];

	    X[nt][j] = X[nt-1][j] + tau * U_F[j]; // motion along the contact discontinuity
	}

//======================THE CORE ITERATION=========================(On Lagrangian Coordinate)
    for(j = 0; j < m; ++j) // forward Euler
	{ /*
	   *  j-1          j          j+1
	   * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	   *   o-----X-----o-----X-----o-----X--...
	   */
	    RHO[nt][j] = 1.0 / (1.0/RHO[nt-1][j] + tau/MASS[j]*(U_F[j+1] - U_F[j]));
	    U[nt][j]   = U[nt-1][j] - tau/MASS[j]*(P_F[j+1] - P_F[j]);
	    E[nt][j]   = E[nt-1][j] - tau/MASS[j]*(P_F[j+1]*U_F[j+1] - P_F[j]*U_F[j]);
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
	    
//============================compute the slopes============================
	    s_u[j]   = (    U_next[j+1] -     U_next[j])/(X[nt][j+1]-X[nt][j]);
	    s_p[j]   = (    P_next[j+1] -     P_next[j])/(X[nt][j+1]-X[nt][j]);
	    s_rho[j] = (RHO_next_L[j+1] - RHO_next_R[j])/(X[nt][j+1]-X[nt][j]);
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
  printf("The cost of CPU time for 1D-GRP Lagrangian scheme for this problem is %g seconds.\n", cpu_time_sum);
//---------------------END OF THE MAIN LOOP----------------------

return_NULL:
  free(s_u);
  free(s_p);
  free(s_rho);
  s_u   = NULL;
  s_p   = NULL;
  s_rho = NULL;
  free(U_next);
  free(P_next);
  free(RHO_next_L);
  free(RHO_next_R);
  U_next     = NULL;
  P_next     = NULL;
  RHO_next_L = NULL;
  RHO_next_R = NULL;
  free(U_t);
  free(P_t);
  free(RHO_t_L);
  free(RHO_t_R);
  U_t     = NULL;
  P_t     = NULL;
  RHO_t_L = NULL;
  RHO_t_R = NULL;
  free(U_F);
  free(P_F);
  U_F = NULL;
  P_F = NULL;
  free(MASS);
  MASS = NULL;
}