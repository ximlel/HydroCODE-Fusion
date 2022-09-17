/**
 * @file  grp_solver_EUL_source.c
 * @brief This is an Eulerian GRP scheme to solve 1-D Euler equations.
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
 * @brief This function use GRP scheme to solve 1-D Euler
 *        equations of motion on Eulerian coordinate.
 * @param[in] m:          Number of the grids.
 * @param[in,out] CV:     Structure of cell variable data.
 * @param[out] cpu_time:  Array of the CPU time recording.
 * @param[out] time_plot: Array of the plotting time recording.
 */
void GRP_solver_EUL_source(const int m, struct cell_var_stru CV, double * cpu_time, const int N_plot, double time_plot[])
{
    /* 
     * j is a frequently used index for spatial variables.
     * k is a frequently used index for the time step.
     */
  int j, k = 0;

  clock_t tic, toc;
  double cpu_time_sum = 0.0;

  double const t_all = config[1];      // the total time
  double const eps   = config[4];      // the largest value could be seen as zero
  int    const N     = (int)config[5]; // the maximum number of time steps
  double const gamma = config[6];      // the constant of the perfect gas
  double const CFL   = config[7];      // the CFL number
  double const h     = config[10];     // the length of the initial spatial grids
  double       tau   = config[16];     // the length of the time step

  _Bool find_bound = false;

  double Mom, Ene;
  double c_L, c_R; // the speeds of sound
  /*
   * dire: the temporal derivative of fluid variables.
   *       \frac{\partial [rho, u, p]}{\partial t}
   * mid:  the Riemann solutions.
   *       [rho_star, u_star, p_star]
   */
  double dire[3], mid[3];

  double nu;  // nu = tau/h
  double h_S_max; // h/S_max, S_max is the maximum wave speed
  double time_c = 0.0; // the current time
  int nt = 1; // the number of times storing plotting data

  struct b_f_var bfv_L = {.SRHO = 0.0, .SP = 0.0, .SU = 0.0}, bfv_R = bfv_L; // Left/Right boundary condition
  struct i_f_var ifv_L = {.gamma = gamma}, ifv_R = ifv_L;

  double ** RHO  = CV.RHO;
  double ** U    = CV.U;
  double ** P    = CV.P;
  double ** E    = CV.E;
  // the slopes of variable values
  double * s_rho = (double*)calloc(m, sizeof(double));
  double * s_u   = (double*)calloc(m, sizeof(double));
  double * s_p   = (double*)calloc(m, sizeof(double));
  CV.d_rho = s_rho;
  CV.d_u   = s_u;
  CV.d_p   = s_p;
  // the variable values at (x_{j-1/2}, t_{n+1}).
  double * U_next   = (double*)malloc((m+1) * sizeof(double));
  double * P_next   = (double*)malloc((m+1) * sizeof(double));
  double * RHO_next = (double*)malloc((m+1) * sizeof(double));
  // the temporal derivatives at (x_{j-1/2}, t_{n}).
  double * U_t   = (double*)malloc((m+1) * sizeof(double));
  double * P_t   = (double*)malloc((m+1) * sizeof(double));
  double * RHO_t = (double*)malloc((m+1) * sizeof(double));
  // the numerical flux at (x_{j-1/2}, t_{n}).
  double * F_rho = (double*)malloc((m+1) * sizeof(double));
  double * F_u   = (double*)malloc((m+1) * sizeof(double));
  double * F_e   = (double*)malloc((m+1) * sizeof(double));
  if(s_rho == NULL || s_u == NULL || s_p == NULL)
      {
	  printf("NOT enough memory! Slope\n");
	  goto return_NULL;
      }
  if(U_next == NULL || P_next == NULL || RHO_next == NULL)
      {
	  printf("NOT enough memory! Variables_next\n");
	  goto return_NULL;
      }
  if(U_t == NULL || P_t == NULL || RHO_t == NULL)
      {
	  printf("NOT enough memory! Temproal derivative\n");
	  goto return_NULL;
      }
  if(F_rho == NULL || F_u == NULL || F_e == NULL)
      {
	  printf("NOT enough memory! Flux\n");
	  goto return_NULL;
      }
  
//-----------------------THE MAIN LOOP--------------------------------
  for(k = 1; k <= N; ++k)
  {
      tic = clock();
      if (time_c > time_plot[nt] && nt < (N_plot-1))
	  nt++;

      h_S_max = INFINITY; // h/S_max = INFINITY

      find_bound = bound_cond_slope_limiter(false, m, nt-1, &CV, &bfv_L, &bfv_R, find_bound, true, time_c);
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
		      ifv_L.RHO = RHO[nt-1][j-1] + 0.5*h*s_rho[j-1];
		      ifv_L.U   =   U[nt-1][j-1] + 0.5*h*s_u[j-1];
		      ifv_L.P   =   P[nt-1][j-1] + 0.5*h*s_p[j-1];
		  }
	      else
		  {
		      ifv_L.RHO = bfv_L.RHO + 0.5*h*bfv_L.SRHO;
		      ifv_L.U   = bfv_L.U   + 0.5*h*bfv_L.SU;
		      ifv_L.P   = bfv_L.P   + 0.5*h*bfv_L.SP;
		  }
	      if(j < m)
		  {
		      ifv_R.RHO = RHO[nt-1][j] - 0.5*h*s_rho[j];
		      ifv_R.U   =   U[nt-1][j] - 0.5*h*s_u[j];
		      ifv_R.P   =   P[nt-1][j] - 0.5*h*s_p[j];
		  }
	      else
		  {
		      ifv_R.RHO = bfv_R.RHO + 0.5*h*bfv_R.SRHO;
		      ifv_R.U   = bfv_R.U   + 0.5*h*bfv_R.SU;
		      ifv_R.P   = bfv_R.P   + 0.5*h*bfv_R.SP;
		  }
	      if(ifv_L.P < eps || ifv_R.P < eps || ifv_L.RHO < eps || ifv_R.RHO < eps)
		  {
		      printf("<0.0 error on [%d, %d] (t_n, x) - Reconstruction\n", k, j);
		      goto return_NULL;
		  }

	      c_L = sqrt(gamma * ifv_L.P / ifv_L.RHO);
	      c_R = sqrt(gamma * ifv_R.P / ifv_R.RHO);
	      h_S_max = fmin(h_S_max, h/(fabs(ifv_L.U)+fabs(c_L)));
	      h_S_max = fmin(h_S_max, h/(fabs(ifv_R.U)+fabs(c_R)));

	      if(j) //calculate the material derivatives
		  {
		      ifv_L.d_u   =   s_u[j-1];
		      ifv_L.d_p   =   s_p[j-1];
		      ifv_L.d_rho = s_rho[j-1];
		  }
	      else
		  {
		      ifv_L.d_rho = bfv_L.SRHO;
		      ifv_L.d_u   = bfv_L.SU;
		      ifv_L.d_p   = bfv_L.SP;
		  }
	      if(j < m)
		  {
		      ifv_R.d_u   =   s_u[j];
		      ifv_R.d_p   =   s_p[j];
		      ifv_R.d_rho = s_rho[j];
		  }
	      else
		  {
		      ifv_R.d_rho = bfv_R.SRHO;
		      ifv_R.d_u   = bfv_R.SU;
		      ifv_R.d_p   = bfv_R.SP;
		  }
	      if(!isfinite(ifv_L.d_p)|| !isfinite(ifv_R.d_p)|| !isfinite(ifv_L.d_u)|| !isfinite(ifv_R.d_u)|| !isfinite(ifv_L.d_rho)|| !isfinite(ifv_R.d_rho))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) - Slope\n", k, j); 
		      goto return_NULL;
		  }

//========================Solve GRP========================
	      linear_GRP_solver_Edir(dire, mid, &ifv_L, &ifv_R, eps, eps);

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
	      if(!isfinite(dire[1])|| !isfinite(dire[2])|| !isfinite(dire[0]))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) - DIRE\n", k, j); 
		      time_c = t_all;
		  }

	      RHO_next[j] = mid[0];
	      U_next[j]   = mid[1];
	      P_next[j]   = mid[2];
	      RHO_t[j] = dire[0];
	      U_t[j]   = dire[1];
	      P_t[j]   = dire[2];
	  }

//====================Time step and grid fixed======================
    // If no total time, use fixed tau and time step N.
    if(isfinite(t_all) || !isfinite(config[16]) || config[16] <= 0.0)
	{
	    tau = CFL * h_S_max;
	    if(tau < eps)
		{
		    printf("\nThe length of the time step is so small on [%d, %g, %g] (t_n, time_c, tau)\n", k, time_c, tau);
		    time_c = t_all;
		}
	    else if((time_c + tau) > (t_all - eps))
		tau = t_all - time_c;
	    else if(!isfinite(tau))
		{
		    printf("NAN or INFinite error on [%d, %g, %g] (t_n, time_c, tau) - CFL\n", k, time_c, tau); 
		    tau = t_all - time_c;
		    goto return_NULL;
		}
	}
    nu = tau / h;
    
    for(j = 0; j <= m; ++j)
	{
	    RHO_next[j] += 0.5 * tau * RHO_t[j];;
	    U_next[j]   += 0.5 * tau * U_t[j];
	    P_next[j]   += 0.5 * tau * P_t[j];

	    F_rho[j] = RHO_next[j]*U_next[j];
	    F_u[j] = F_rho[j]*U_next[j] + P_next[j];
	    F_e[j] = (gamma/(gamma-1.0))*P_next[j] + 0.5*F_rho[j]*U_next[j];
	    F_e[j] = F_e[j]*U_next[j];

	    RHO_next[j] += 0.5 * tau * RHO_t[j];;
	    U_next[j]   += 0.5 * tau * U_t[j];
	    P_next[j]   += 0.5 * tau * P_t[j];
	}

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
	    
//============================compute the slopes============================
	    s_u[j]   = (  U_next[j+1] -   U_next[j])/h;
	    s_p[j]   = (  P_next[j+1] -   P_next[j])/h;
	    s_rho[j] = (RHO_next[j+1] - RHO_next[j])/h;
	}

//============================Time update=======================

    time_c += tau;
    if(isfinite(t_all))
        DispPro(time_c*100.0/t_all, k);
    else
        DispPro(k*100.0/N, k);
    if(time_c > (t_all - eps) || isinf(time_c))
	break;

//===========================Fixed variable location=======================	
    for(j = 0; j < m; ++j)
	{
	    RHO[nt-1][j] = RHO[nt][j];
	    U[nt-1][j]   =   U[nt][j];
	    E[nt-1][j]   =   E[nt][j];  
	    P[nt-1][j]   =   P[nt][j];
	}

    toc = clock();
    cpu_time[nt] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    cpu_time_sum += cpu_time[nt];
  }

  printf("\nTime is up at time step %d.\n", k);
  printf("The cost of CPU time for 1D-GRP Eulerian scheme for this problem is %g seconds.\n", cpu_time_sum);
//---------------------END OF THE MAIN LOOP----------------------

return_NULL:
  config[5] = (double)k;
  if(fabs(time_plot[1]) < eps)
      {
	  time_plot[1] = time_c;
	  if(isfinite(time_c))
	      {
		  time_plot[N_plot-2] = time_c - tau;
		  time_plot[N_plot-1] = time_c;
	      }
	  else
	      {
		  time_plot[N_plot-2] = N*tau - tau;
		  time_plot[N_plot-1] = N*tau;
	      }
      }

  free(s_u);
  free(s_p);
  free(s_rho);
  s_u   = NULL;
  s_p   = NULL;
  s_rho = NULL;
  free(U_next);
  free(P_next);
  free(RHO_next);
  U_next   = NULL;
  P_next   = NULL;
  RHO_next = NULL;
  free(U_t);
  free(P_t);
  free(RHO_t);
  U_t   = NULL;
  P_t   = NULL;
  RHO_t = NULL;
  free(F_rho);
  free(F_u);
  free(F_e);
  F_rho = NULL;
  F_u   = NULL;
  F_e   = NULL;
}
