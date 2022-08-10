/**
 * @file  GRP_solver_source.c
 * @brief This is a Lagrangian/Eulerian GRP scheme to solve 1-D Euler equations.
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "../include/Riemann_solver.h"

#ifdef _WIN32
#define ISNAN(a) _isnan((a))
#define ISERR(a) !_finite((a))
#elif __linux__
#define ISNAN(a)    isnan((a))
#define ISERR(a) !isfinite((a))
#endif

/**
 * @brief Minmod limiter of two variables.
 */
static inline double minmod2(double s_L, double s_R)
{
    if(s_L * s_R < 0.0)
	return 0.0;
    else if(fabs(s_R) < fabs(s_L))
	return s_R;
    else
	return s_L;
}

/**
 * @brief Minmod limiter of three variables.
 */
static inline double minmod3(double s_L, double s_R, double s_m)
{
    if(s_L * s_m < 0.0 || s_R * s_m < 0.0)
	return 0.0;
    else if(fabs(s_m) < fabs(s_L) && fabs(s_m) < fabs(s_R))
	return s_m;
    else if(fabs(s_R) < fabs(s_L))
	return s_R;
    else
	return s_L;
}

/**
 * @brief This function use GRP scheme to solve 1-D Euler
 *        equations of motion on Lagrangian coordinate.
 * @param[in]  config:   The array of configuration data.
 * @param[in]  m:        The number of the grids.
 * @param[in,out] RHO,U,P,E,X[]: Array of the density/velocity/pressure/energy/coordinate data.
 * @param[out] cpu_time: Array of the CPU time recording.
 */
void GRP_solver_LAG_source
(double * config, const int m, double * RHO[], double * U[], double * P[],
 double * E[], double * X[], double * cpu_time)
{
    /* 
     * j is a frequently used index for spatial variables.
     * k is a frequently used index for the time step.
     */
  int j, k;

  clock_t tic, toc;
  double cpu_time_sum = 0.0;

  double const gamma = config[0];        // the constant of the perfect gas
  double       tau   = config[1];        // the length of the time step
  double const h     = config[2];        // the length of the initial spatial grids
  double const eps   = config[3];        // the largest value could be seen as zero
  int    const N     = (int)(config[4]); // the number of time steps
  double const t_all = config[5];        // the total time
  double const CFL   = config[6];        // the CFL number
  int    const bound = (int)(config[7]); // the boundary condition

  double u_L, p_L, rho_L;
  double u_R, p_R, rho_R;
  double c_L, c_R; // the speeds of sound
  double h_L, h_R; // length of spatial grids

  double s_u_L, s_p_L, s_rho_L; // spatial derivatives in coordinate x (slopes)
  double s_u_R, s_p_R, s_rho_R;
  double t_u_L, t_p_L, t_rho_L; // spatial derivatives in Lagrangian coordinate ξ
  double t_u_R, t_p_R, t_rho_R;
  /*
   * dire: the temporal derivative of fluid variables.
   *       \frac{\partial [rho_L, u, p, rho_R]}{\partial t}
   * mid:  the Riemann solutions.
   *       [rho_star_L, u_star, p_star, rho_star_R]
   */
  double dire[4], mid[4];
  // the paramater in slope limiters.
  double alpha = 1.5;

  // the slopes of variable values
  double *s_rho, *s_u, *s_p;
  // the variable values at (x_{j-1/2}, t_{n+1}).
  double *U_next, *P_next, *RHO_next_L, *RHO_next_R;
  // the temporal derivatives at (x_{j-1/2}, t_{n}).
  double *U_t, *P_t, *RHO_t_L, *RHO_t_R;
  double *U_F, *P_F; // the numerical flux at (x_{j-1/2}, t_{n+1/2}).
  double * MASS; // Array of the mass data in computational cells.
  s_rho = calloc(m, sizeof(double));
  s_u   = calloc(m, sizeof(double));
  s_p   = calloc(m, sizeof(double));
  if(s_rho == NULL || s_u == NULL || s_p == NULL)
      {
	  printf("NOT enough memory! Slope\n");
	  goto _END_;
      }
  U_next     = calloc(m+1, sizeof(double));
  P_next     = calloc(m+1, sizeof(double));
  RHO_next_L = calloc(m+1, sizeof(double));
  RHO_next_R = calloc(m+1, sizeof(double));
  if(U_next == NULL || P_next == NULL || RHO_next_L == NULL || RHO_next_R == NULL)
      {
	  printf("NOT enough memory! Variables_next\n");
	  goto _END_;
      }
  U_t     = calloc(m+1, sizeof(double));
  P_t     = calloc(m+1, sizeof(double));
  RHO_t_L = calloc(m+1, sizeof(double));
  RHO_t_R = calloc(m+1, sizeof(double));
  if(U_t == NULL || P_t == NULL || RHO_t_L == NULL || RHO_t_R == NULL)
      {
	  printf("NOT enough memory! Temproal derivative\n");
	  goto _END_;
      }
  U_F = calloc(m+1, sizeof(double));
  P_F = calloc(m+1, sizeof(double));
  if(U_F == NULL || P_F == NULL)
      {
	  printf("NOT enough memory! Variables_F\n");
	  goto _END_;
      }
  MASS = malloc(m * sizeof(double));
  if(MASS == NULL)
      {
	  printf("NOT enough memory! MASS\n");
	  goto _END_;
      }
  for(k = 0; k < m; ++k) // Initialize the values of mass in computational cells
      MASS[k] = h * RHO[0][k];

  double h_S_max; // h/S_max, S_max is the maximum wave speed
  double time_c = 0.0; // the current time
  int n = 1; // the number of times storing plotting data

  double UL, PL, RHOL, HL, SUL, SPL, SRHOL; // Left  boundary condition
  double UR, PR, RHOR, HR, SUR, SPR, SRHOR; // Right boundary condition
  if (bound == -1) // initial boudary conditions
      {
	  UL   =   U[0][0]; UR   =   U[0][m-1];
	  PL   =   P[0][0]; PR   =   P[0][m-1];
	  RHOL = RHO[0][0]; RHOR = RHO[0][m-1];
	  HL  = h; HR = h;
      }
  SUL   = 0.0;   SUR = 0.0;
  SPL   = 0.0;   SPR = 0.0;
  SRHOL = 0.0; SRHOR = 0.0;

//-----------------------THE MAIN LOOP--------------------------------
  for(k = 1; k <= N; ++k)
  {
      h_S_max = INFINITY; // h/S_max = INF
      tic = clock();

      if (bound == -2) // reflective boundary conditions
	  {
	      UL   = - U[n-1][0]; UR   = - U[n-1][m-1];
	      PL   =   P[n-1][0]; PR   =   P[n-1][m-1];
	      RHOL = RHO[n-1][0]; RHOR = RHO[n-1][m-1];
	      HL = X[n-1][1] - X[n-1][0];
	      HR = X[n-1][m] - X[n-1][m-1];
	}
      if (bound == -4) // free boundary conditions
	  {
	      UL   =   U[n-1][0]; UR   =   U[n-1][m-1];
	      PL   =   P[n-1][0]; PR   =   P[n-1][m-1];
	      RHOL = RHO[n-1][0]; RHOR = RHO[n-1][m-1];
	      HL = X[n-1][1] - X[n-1][0];
	      HR = X[n-1][m] - X[n-1][m-1];
	  }
      if (bound == -5) // periodic boundary conditions
	  {
	      UL   =   U[n-1][m-1]; UR   =   U[n-1][0];
	      PL   =   P[n-1][m-1]; PR   =   P[n-1][0];
	      RHOL = RHO[n-1][m-1]; RHOR = RHO[n-1][0];
	      HL = X[n-1][m] - X[n-1][m-1];
	      HR = X[n-1][1] - X[n-1][0];
	  }
      
//=================Initialize slopes=====================
      for(j = 0; j < m; ++j) // Reconstruct slopes
	  { /*
	     *  j-1          j          j+1
	     * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	     *   o-----X-----o-----X-----o-----X--...
	     */
	      if(j)
		  {
		      h_L     = 0.5 * (X[n-1][j+1] - X[n-1][j-1]);
		      s_u_L   = (U[n-1][j]   -   U[n-1][j-1]) / h_L;
		      s_p_L   = (P[n-1][j]   -   P[n-1][j-1]) / h_L;
		      s_rho_L = (RHO[n-1][j] - RHO[n-1][j-1]) / h_L;
		  }
	      else
		  {
		      h_L     = 0.5 * (X[n-1][j+1] - X[n-1][j] + HL);
		      s_u_L   = (U[n-1][j]   -   UL) / h_L;
		      s_p_L   = (P[n-1][j]   -   PL) / h_L;
		      s_rho_L = (RHO[n-1][j] - RHOL) / h_L;
		  }
	      if(j < m-1)
		  {
		      h_R     = 0.5 * (X[n-1][j+2] - X[n-1][j]);
		      s_u_R   = (U[n-1][j+1]   -   U[n-1][j]) / h_R;
		      s_p_R   = (P[n-1][j+1]   -   P[n-1][j]) / h_R;
		      s_rho_R = (RHO[n-1][j+1] - RHO[n-1][j]) / h_R;
		  }
	      else
		  {
		      h_R     = 0.5 * (X[n-1][j+1] - X[n-1][j] + HR);
		      s_u_R   = (UR   -   U[n-1][j]) / h_R;
		      s_p_R   = (PR   -   P[n-1][j]) / h_R;
		      s_rho_R = (RHOR - RHO[n-1][j]) / h_R;
		  }
	      if (k == 1)
		  {
		      s_u[j]   = minmod2(s_u_L,   s_u_R);
		      s_p[j]   = minmod2(s_p_L,   s_p_R);
		      s_rho[j] = minmod2(s_rho_L, s_rho_R);
		  }
	      else
		  {
		      s_u[j]   = minmod3(alpha*s_u_L,   alpha*s_u_R,   s_u[j]);
		      s_p[j]   = minmod3(alpha*s_p_L,   alpha*s_p_R,   s_p[j]);
		      s_rho[j] = minmod3(alpha*s_rho_L, alpha*s_rho_R, s_rho[j]);
		  }
	  }
      if (bound == -2) // reflective boundary conditions
	  {
	      SUL = - s_u[0]; SUR = - s_u[m-1];
	  }
      if (bound == -5) // periodic boundary conditions
	  {
	      SUL   =   s_u[m-1]; SUR   =   s_u[0];
	      SPL   =   s_p[m-1]; SPR   =   s_p[0];
	      SRHOL = s_rho[m-1]; SRHOR = s_rho[0];
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
		      rho_L = RHO[n-1][j-1] + 0.5*h_L*s_rho[j-1];
		      u_L   =   U[n-1][j-1] + 0.5*h_L*s_u[j-1];
		      p_L   =   P[n-1][j-1] + 0.5*h_L*s_p[j-1];
		  }
	      else
		  {
		      h_L   =   HL;
		      rho_L = RHOL + 0.5*h_L*SRHOL;
		      u_L   =   UL + 0.5*h_L*SUL;
		      p_L   =   PL + 0.5*h_L*SPL;
		  }
	      if(j < m)
		  {
		      h_R   =   X[n-1][j+1] - X[n-1][j];
		      rho_R = RHO[n-1][j] - 0.5*h_R*s_rho[j];
		      u_R   =   U[n-1][j] - 0.5*h_R*s_u[j];
		      p_R   =   P[n-1][j] - 0.5*h_R*s_p[j];
		  }
	      else
		  {
		      h_R   =   HR;
		      rho_R = RHOR + 0.5*h_R*SRHOR;
		      u_R   =   UR + 0.5*h_R*SUR;
		      p_R   =   PR + 0.5*h_R*SPR;
		  }
	      if(p_L < eps || p_R < eps || rho_L < eps || rho_R < eps)
		  {
		      printf("<0.0 error on [%d, %d] (t_n, x) - Reconstruction\n", k, j);
		      goto _END_;
		  }
	      if(ISERR(p_L)||ISERR(p_R)||ISERR(u_L)||ISERR(u_R)||ISERR(rho_L)||ISERR(rho_R))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) - Reconstruction\n", k, j); 
		      goto _END_;
		  }

	      c_L = sqrt(gamma * p_L / rho_L);
	      c_R = sqrt(gamma * p_R / rho_R);
	      h_S_max = fmin(h_S_max, h_L/c_L);
	      h_S_max = fmin(h_S_max, h_R/c_R);

	      if(j) //calculate the material derivatives
		  {
		      t_u_L   =   s_u[j-1]/rho_L;
		      t_p_L   =   s_p[j-1]/rho_L;
		      t_rho_L = s_rho[j-1]/rho_L;
		  }
	      else
		  {
		      t_rho_L = SRHOL/rho_L;
		      t_u_L   =   SUL/rho_L;
		      t_p_L   =   SPL/rho_L;
		  }
	      if(j < m)
		  {
		      t_u_R   =   s_u[j]/rho_R;
		      t_p_R   =   s_p[j]/rho_R;
		      t_rho_R = s_rho[j]/rho_R;
		  }
	      else
		  {
		      t_rho_R = SRHOR/rho_R;
		      t_u_R   =   SUR/rho_R;
		      t_p_R   =   SPR/rho_R;
		  }

//========================Solve GRP========================
	      linear_GRP_solver_LAG(dire, mid, rho_L, rho_R, t_rho_L, t_rho_R, u_L, u_R, t_u_L, t_u_R, p_L, p_R, t_p_L, t_p_R, gamma, eps, eps);

	      if(mid[2] < eps)
		  {
		      printf("<0.0 error on [%d, %d] (t_n, x) STAR\n", k, j);
		      time_c = t_all;
		  }
	      if(ISERR(mid[1])||ISERR(mid[2]))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) STAR\n", k, j); 
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
    tau = CFL * h_S_max;
    if ((time_c + tau) > (t_all - eps))
        tau = t_all - time_c;
    
    for(j = 0; j <= m; ++j)
	{
	    U_F[j] = U_next[j] + 0.5 * tau * U_t[j];
	    P_F[j] = P_next[j] + 0.5 * tau * P_t[j];         

	    RHO_next_L[j] += tau * RHO_t_L[j];
	    RHO_next_R[j] += tau * RHO_t_R[j];
	    U_next[j]     += tau * U_t[j];
	    P_next[j]     += tau * P_t[j];

	    X[n][j] = X[n-1][j] + tau * U_F[j]; // motion along the contact discontinuity
	}

//======================THE CORE ITERATION=========================(On Lagrangian Coordinate)
    for(j = 0; j < m; ++j) // forward Euler
	{ /*
	   *  j-1          j          j+1
	   * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	   *   o-----X-----o-----X-----o-----X--...
	   */
	    RHO[n][j] = 1 / (1/RHO[n-1][j] + tau/MASS[j]*(U_F[j+1] - U_F[j]));
	    U[n][j]   = U[n-1][j] - tau/MASS[j]*(P_F[j+1] - P_F[j]);
	    E[n][j]   = E[n-1][j] - tau/MASS[j]*(P_F[j+1]*U_F[j+1] - P_F[j]*U_F[j]);
	    P[n][j]   = (E[n][j] - 0.5 * U[n][j]*U[n][j]) * (gamma - 1.0) * RHO[n][j];
	    if(P[n][j] < eps || RHO[n][j] < eps)
		{
		    printf("<0.0 error on [%d, %d] (t_n, x)\n", k, j);
		    time_c = t_all;
		}
	    if(ISERR(P[n][j])||ISERR(U[n][j])||ISERR(RHO[n][j]))
		{
		    printf("NAN or INFinite error on [%d, %d] (t_n, x)\n", k, j); 
		    time_c = t_all;
		}
	    
//============================compute the slopes============================
	    s_u[j]   = (    U_next[j+1] -     U_next[j])/(X[n][j+1]-X[n][j]);
	    s_p[j]   = (    P_next[j+1] -     P_next[j])/(X[n][j+1]-X[n][j]);
	    s_rho[j] = (RHO_next_L[j+1] - RHO_next_R[j])/(X[n][j+1]-X[n][j]);
	}

//============================Time update=======================

    toc = clock();
    cpu_time[n] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    cpu_time_sum += cpu_time[n];

    time_c += tau;
    if(time_c > (t_all - eps))
	{
	    printf("Time is up in time step %d.\n", k);
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

  printf("The cost of CPU time for 1D-GRP scheme for this problem is %g seconds.\n", cpu_time_sum);
//---------------------END OF THE MAIN LOOP----------------------

_END_:
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

/**
 * @brief This function use GRP scheme to solve 1-D Euler
 *        equations of motion on Eulerian coordinate.
 * @param[in]  config:   The array of configuration data.
 * @param[in]  m:        The number of the grids.
 * @param[in,out] RHO,U,P,E,X[]: Array of the density/velocity/pressure/energy/coordinate data.
 * @param[out] cpu_time: Array of the CPU time recording.
 */
void GRP_solver_EUL_source
(double * config, const int m, double * RHO[], double * U[], double * P[],
 double * E[], double * X[], double * cpu_time)
{
    /* 
     * j is a frequently used index for spatial variables.
     * k is a frequently used index for the time step.
     */
  int j, k;

  clock_t tic, toc;
  double cpu_time_sum = 0.0;

  double const gamma = config[0];        // the constant of the perfect gas
  double       tau   = config[1];        // the length of the time step
  double const h     = config[2];        // the length of the initial spatial grids
  double const eps   = config[3];        // the largest value could be seen as zero
  int    const N     = (int)(config[4]); // the number of time steps
  double const t_all = config[5];        // the total time
  double const CFL   = config[6];        // the CFL number
  int    const bound = (int)(config[7]); // the boundary condition

  double Mom, Ene;
  double u_L, p_L, rho_L;
  double u_R, p_R, rho_R;
  double c_L, c_R; // the speeds of sound
  double h_L, h_R; // length of spatial grids

  double s_u_L, s_p_L, s_rho_L; // spatial derivatives in coordinate x (slopes)
  double s_u_R, s_p_R, s_rho_R;
  /*
   * dire: the temporal derivative of fluid variables.
   *       \frac{\partial [rho, u, p]}{\partial t}
   * mid:  the Riemann solutions.
   *       [rho_star, u_star, p_star]
   */
  double dire[3], mid[3];
  // the paramater in slope limiters.
  double alpha = 1.5;

  // the slopes of variable values
  double *s_rho, *s_u, *s_p;
  // the variable values at (x_{j-1/2}, t_{n+1}).
  double *U_next, *P_next, *RHO_next;
  // the temporal derivatives at (x_{j-1/2}, t_{n}).
  double *U_t, *P_t, *RHO_t;
  double nu;  // nu = tau/h
  double *F1, *F2, *F3; // the numerical flux at (x_{j-1/2}, t_{n}).
  double * MASS; // Array of the mass data in computational cells.
  s_rho = calloc(m, sizeof(double));
  s_u   = calloc(m, sizeof(double));
  s_p   = calloc(m, sizeof(double));
  if(s_rho == NULL || s_u == NULL || s_p == NULL)
      {
	  printf("NOT enough memory! Slope\n");
	  goto _END_;
      }
  U_next   = calloc(m+1, sizeof(double));
  P_next   = calloc(m+1, sizeof(double));
  RHO_next = calloc(m+1, sizeof(double));
  if(U_next == NULL || P_next == NULL || RHO_next == NULL)
      {
	  printf("NOT enough memory! Variables_next\n");
	  goto _END_;
      }
  U_t   = calloc(m+1, sizeof(double));
  P_t   = calloc(m+1, sizeof(double));
  RHO_t = calloc(m+1, sizeof(double));
  if(U_t == NULL || P_t == NULL || RHO_t == NULL)
      {
	  printf("NOT enough memory! Temproal derivative\n");
	  goto _END_;
      }
  F1 = malloc((m+1) * sizeof(double));
  F2 = malloc((m+1) * sizeof(double));
  F3 = malloc((m+1) * sizeof(double));
  if(F1 == NULL || F2 == NULL || F3 == NULL)
      {
	  printf("NOT enough memory! Flux\n");
	  goto _END_;
      }
  MASS = malloc(m * sizeof(double));
  if(MASS == NULL)
      {
	  printf("NOT enough memory! MASS\n");
	  goto _END_;
      }
  for(k = 0; k < m; ++k) // Initialize the values of mass in computational cells
      MASS[k] = h * RHO[0][k];

  double h_S_max; // h/S_max, S_max is the maximum wave speed
  double time_c = 0.0; // the current time
  int n = 1; // the number of times storing plotting data

  double UL, PL, RHOL, HL, SUL, SPL, SRHOL; // Left  boundary condition
  double UR, PR, RHOR, HR, SUR, SPR, SRHOR; // Right boundary condition
  if (bound == -1) // initial boudary conditions
      {
	  UL   =   U[0][0]; UR   =   U[0][m-1];
	  PL   =   P[0][0]; PR   =   P[0][m-1];
	  RHOL = RHO[0][0]; RHOR = RHO[0][m-1];
	  HL  = h; HR = h;
      }
  SUL   = 0.0;   SUR = 0.0;
  SPL   = 0.0;   SPR = 0.0;
  SRHOL = 0.0; SRHOR = 0.0;

//-----------------------THE MAIN LOOP--------------------------------
  for(k = 1; k <= N; ++k)
  {
      h_S_max = INFINITY; // h/S_max = INF
      tic = clock();

      if (bound == -2) // reflective boundary conditions
	  {
	      UL   = - U[n-1][0]; UR   = - U[n-1][m-1];
	      PL   =   P[n-1][0]; PR   =   P[n-1][m-1];
	      RHOL = RHO[n-1][0]; RHOR = RHO[n-1][m-1];
	      HL = X[n-1][1] - X[n-1][0];
	      HR = X[n-1][m] - X[n-1][m-1];
	}
      if (bound == -4) // free boundary conditions
	  {
	      UL   =   U[n-1][0]; UR   =   U[n-1][m-1];
	      PL   =   P[n-1][0]; PR   =   P[n-1][m-1];
	      RHOL = RHO[n-1][0]; RHOR = RHO[n-1][m-1];
	      HL = X[n-1][1] - X[n-1][0];
	      HR = X[n-1][m] - X[n-1][m-1];
	  }
      if (bound == -5) // periodic boundary conditions
	  {
	      UL   =   U[n-1][m-1]; UR   =   U[n-1][0];
	      PL   =   P[n-1][m-1]; PR   =   P[n-1][0];
	      RHOL = RHO[n-1][m-1]; RHOR = RHO[n-1][0];
	      HL = X[n-1][m] - X[n-1][m-1];
	      HR = X[n-1][1] - X[n-1][0];
	  }
      
//=================Initialize slopes=====================
      for(j = 0; j < m; ++j) // Reconstruct slopes
	  { /*
	     *  j-1          j          j+1
	     * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	     *   o-----X-----o-----X-----o-----X--...
	     */
	      if(j)
		  {
		      h_L     = 0.5 * (X[n-1][j+1] - X[n-1][j-1]);
		      s_u_L   = (U[n-1][j]   -   U[n-1][j-1]) / h_L;
		      s_p_L   = (P[n-1][j]   -   P[n-1][j-1]) / h_L;
		      s_rho_L = (RHO[n-1][j] - RHO[n-1][j-1]) / h_L;
		  }
	      else
		  {
		      h_L     = 0.5 * (X[n-1][j+1] - X[n-1][j] + HL);
		      s_u_L   = (U[n-1][j]   -   UL) / h_L;
		      s_p_L   = (P[n-1][j]   -   PL) / h_L;
		      s_rho_L = (RHO[n-1][j] - RHOL) / h_L;
		  }
	      if(j < m-1)
		  {
		      h_R     = 0.5 * (X[n-1][j+2] - X[n-1][j]);
		      s_u_R   = (U[n-1][j+1]   -   U[n-1][j]) / h_R;
		      s_p_R   = (P[n-1][j+1]   -   P[n-1][j]) / h_R;
		      s_rho_R = (RHO[n-1][j+1] - RHO[n-1][j]) / h_R;
		  }
	      else
		  {
		      h_R     = 0.5 * (X[n-1][j+1] - X[n-1][j] + HR);
		      s_u_R   = (UR   -   U[n-1][j]) / h_R;
		      s_p_R   = (PR   -   P[n-1][j]) / h_R;
		      s_rho_R = (RHOR - RHO[n-1][j]) / h_R;
		  }
	      if (k == 1)
		  {
		      s_u[j]   = minmod2(s_u_L,   s_u_R);
		      s_p[j]   = minmod2(s_p_L,   s_p_R);
		      s_rho[j] = minmod2(s_rho_L, s_rho_R);
		  }
	      else
		  {
		      s_u[j]   = minmod3(alpha*s_u_L,   alpha*s_u_R,   s_u[j]);
		      s_p[j]   = minmod3(alpha*s_p_L,   alpha*s_p_R,   s_p[j]);
		      s_rho[j] = minmod3(alpha*s_rho_L, alpha*s_rho_R, s_rho[j]);
		  }
	  }
      if (bound == -2) // reflective boundary conditions
	  {
	      SUL = - s_u[0]; SUR = - s_u[m-1];
	  }
      if (bound == -5) // periodic boundary conditions
	  {
	      SUL   =   s_u[m-1]; SUR   =   s_u[0];
	      SPL   =   s_p[m-1]; SPR   =   s_p[0];
	      SRHOL = s_rho[m-1]; SRHOR = s_rho[0];
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
		      rho_L = RHO[n-1][j-1] + 0.5*h_L*s_rho[j-1];
		      u_L   =   U[n-1][j-1] + 0.5*h_L*s_u[j-1];
		      p_L   =   P[n-1][j-1] + 0.5*h_L*s_p[j-1];
		  }
	      else
		  {
		      h_L   =   HL;
		      rho_L = RHOL + 0.5*h_L*SRHOL;
		      u_L   =   UL + 0.5*h_L*SUL;
		      p_L   =   PL + 0.5*h_L*SPL;
		  }
	      if(j < m)
		  {
		      h_R   =   X[n-1][j+1] - X[n-1][j];
		      rho_R = RHO[n-1][j] - 0.5*h_R*s_rho[j];
		      u_R   =   U[n-1][j] - 0.5*h_R*s_u[j];
		      p_R   =   P[n-1][j] - 0.5*h_R*s_p[j];
		  }
	      else
		  {
		      h_R   =   HR;
		      rho_R = RHOR + 0.5*h_R*SRHOR;
		      u_R   =   UR + 0.5*h_R*SUR;
		      p_R   =   PR + 0.5*h_R*SPR;
		  }
	      if(p_L < eps || p_R < eps || rho_L < eps || rho_R < eps)
		  {
		      printf("<0.0 error on [%d, %d] (t_n, x) - Reconstruction\n", k, j);
		      goto _END_;
		  }
	      if(ISERR(p_L)||ISERR(p_R)||ISERR(u_L)||ISERR(u_R)||ISERR(rho_L)||ISERR(rho_R))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) - Reconstruction\n", k, j); 
		      goto _END_;
		  }

	      c_L = sqrt(gamma * p_L / rho_L);
	      c_R = sqrt(gamma * p_R / rho_R);
	      h_S_max = fmin(h_S_max, h_L/(fabs(u_L)+fabs(c_L)));
	      h_S_max = fmin(h_S_max, h_R/(fabs(u_R)+fabs(c_R)));

	      if(j) //calculate the material derivatives
		  {
		      s_u_L   =   s_u[j-1];
		      s_p_L   =   s_p[j-1];
		      s_rho_L = s_rho[j-1];
		  }
	      else
		  {
		      s_rho_L = SRHOL;
		      s_u_L   =   SUL;
		      s_p_L   =   SPL;
		  }
	      if(j < m)
		  {
		      s_u_R   =   s_u[j];
		      s_p_R   =   s_p[j];
		      s_rho_R = s_rho[j];
		  }
	      else
		  {
		      s_rho_R = SRHOR;
		      s_u_R   =   SUR;
		      s_p_R   =   SPR;
		  }

//========================Solve GRP========================
	      linear_GRP_solver_Edir(dire, mid, rho_L, rho_R, s_rho_L, s_rho_R, u_L, u_R, s_u_L, s_u_R, p_L, p_R, s_p_L, s_p_R, gamma, eps);

	      if(mid[2] < eps)
		  {
		      printf("<0.0 error on [%d, %d] (t_n, x) STAR\n", k, j);
		      time_c = t_all;
		  }
	      if(ISERR(mid[1])||ISERR(mid[2]))
		  {
		      printf("NAN or INFinite error on [%d, %d] (t_n, x) STAR\n", k, j); 
		      time_c = t_all;
		  }
	      
	      RHO_next[j] = mid[0];
	      U_next[j]   = mid[1];
	      P_next[j]   = mid[2];
	      RHO_t[j] = dire[0];
	      U_t[j]   = dire[1];
	      P_t[j]   = dire[2];
	  }

//====================Time step and grid movement======================
    tau = CFL * h_S_max;
    if ((time_c + tau) > (t_all - eps))
        tau = t_all - time_c;
    nu = tau / h;
    
    for(j = 0; j <= m; ++j)
	{
	    RHO_next[j] += 0.5 * tau * RHO_t[j];;
	    U_next[j]   += 0.5 * tau * U_t[j];
	    P_next[j]   += 0.5 * tau * P_t[j];

	    F1[j] = RHO_next[j]*U_next[j];
	    F2[j] = F1[j]*U_next[j] + P_next[j];
	    F3[j] = (gamma/(gamma-1.0))*P_next[j] + 0.5*F1[j]*U_next[j];
	    F3[j] = F3[j]*U_next[j];

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
	    RHO[n][j] = RHO[n-1][j] - nu*(F1[j+1]-F1[j]);
	    Mom = RHO[n-1][j]*U[n-1][j] - nu*(F2[j+1]-F2[j]);
	    Ene = RHO[n-1][j]*E[n-1][j] - nu*(F3[j+1]-F3[j]);
	    
	    U[n][j] = Mom / RHO[n][j];
	    E[n][j] = Ene / RHO[n][j];
	    P[n][j] = (Ene - 0.5*Mom*U[n][j])*(gamma-1.0);

	    if(P[n][j] < eps || RHO[n][j] < eps)
		{
		    printf("<0.0 error on [%d, %d] (t_n, x)\n", k, j);
		    time_c = t_all;
		}
	    if(ISERR(P[n][j])||ISERR(U[n][j])||ISERR(RHO[n][j]))
		{
		    printf("NAN or INFinite error on [%d, %d] (t_n, x)\n", k, j); 
		    time_c = t_all;
		}
	    
//============================compute the slopes============================
	    s_u[j]   = (  U_next[j+1] -   U_next[j])/(X[n][j+1]-X[n][j]);
	    s_p[j]   = (  P_next[j+1] -   P_next[j])/(X[n][j+1]-X[n][j]);
	    s_rho[j] = (RHO_next[j+1] - RHO_next[j])/(X[n][j+1]-X[n][j]);
	}

//============================Time update=======================

    toc = clock();
    cpu_time[n] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    cpu_time_sum += cpu_time[n];

    time_c += tau;
    if(time_c > (t_all - eps))
	{
	    printf("Time is up in time step %d.\n", k);
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

  printf("The cost of CPU time for 1D-GRP scheme for this problem is %g seconds.\n", cpu_time_sum);
//---------------------END OF THE MAIN LOOP----------------------

_END_:
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
  free(F1);
  free(F2);
  free(F3);
  F1 = NULL;
  F2 = NULL;
  F3 = NULL;
  free(MASS);
  MASS = NULL;
}
