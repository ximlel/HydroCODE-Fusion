/**
 * @file  GRP_solver_2D_EUL_source.c
 * @brief This is an Eulerian GRP scheme to solve 2-D Euler equations.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#include "../include/var_struc.h"
#include "../include/Riemann_solver.h"
#include "../include/tools.h"


#define _2D_INIT_MEM(v, M, N)						\
    do {								\
	v = (double **)malloc((M) * sizeof(double *));			\
	if(v == NULL)							\
	    {								\
		printf("NOT enough memory! %s\n", #v);			\
		goto return_NULL;					\
	    }								\
	for(j = 0; j < (M); ++j)					\
	    {								\
		v[j] = (double *)malloc((N) * sizeof(double));		\
		if(v[j] == NULL)					\
		    {							\
			printf("NOT enough memory! %s[%d]\n", #v, j);	\
			goto return_NULL;				\
		    }							\
	    }								\
    } while (0)

/**
 * @brief This function use GRP scheme to solve 2-D Euler
 *        equations of motion on Eulerian coordinate.
 * @param[in]  m:        Number of the x-grids: n_x.
 * @param[in]  n:        Number of the y-grids: n_y.
 * @param[in,out] CV:    Structural body of cell variable data.
 * @param[out] cpu_time: Array of the CPU time recording.
 */
void GRP_solver_2D_EUL_source(const int m, const int n, struct cell_var_stru * CV, double * cpu_time)
{
    /* 
     * i is a frequently used index for spatial variables.
     * j is a frequently used index for x - spatial variables.
     * k is a frequently used index for the time step.
     */
  int i, j, k;

  clock_t tic, toc;
  double cpu_time_sum = 0.0;

  double const t_all   = config[1];        // the total time
  double const eps     = config[4];        // the largest value could be seen as zero
  int    const N       = (int)(config[5]); // the maximum number of time steps
  double const gamma   = config[6];        // the constant of the perfect gas
  double const CFL     = config[7];        // the CFL number
  double const h_x     = config[10];       // the length of the initial x-spatial grids
  double const h_y     = config[11];       // the length of the initial y-spatial grids
  double       tau     = config[16];       // the length of the time step
  int    const bound_x = (int)(config[17]);// the boundary condition in x-direction
  int    const bound_y = (int)(config[18]);// the boundary condition in y-direction
  double const alpha   = config[41];       // the paramater in slope limiters.

  double mom_x, mom_y, ene;
  double c;

  double dire[4], mid[4];

  double ** s_rho, ** s_u, ** s_v, ** s_p;
  double ** t_rho, ** t_u, ** t_v, ** t_p;
  double ** rhox, ** ux, ** vx, ** px;
  double ** F1, ** F2, ** F3, ** F4;
  double ** rhoy, ** uy, ** vy, ** py;
  double ** G1, ** G2, ** G3, ** G4;
  _2D_INIT_MEM(s_rho, m, n);
  _2D_INIT_MEM(s_u,   m, n);
  _2D_INIT_MEM(s_v,   m, n);
  _2D_INIT_MEM(s_p,   m, n);
  _2D_INIT_MEM(t_rho, m, n);
  _2D_INIT_MEM(t_u,   m, n);
  _2D_INIT_MEM(t_v,   m, n);
  _2D_INIT_MEM(t_p,   m, n);
  _2D_INIT_MEM(rhox, m+1, n);
  _2D_INIT_MEM(ux,   m+1, n);
  _2D_INIT_MEM(vx,   m+1, n);
  _2D_INIT_MEM(px,   m+1, n);
  _2D_INIT_MEM(F1,   m+1, n);
  _2D_INIT_MEM(F2,   m+1, n);
  _2D_INIT_MEM(F3,   m+1, n);
  _2D_INIT_MEM(F4,   m+1, n); 
  _2D_INIT_MEM(rhoy, m, n+1);
  _2D_INIT_MEM(uy,   m, n+1);
  _2D_INIT_MEM(vy,   m, n+1);
  _2D_INIT_MEM(py,   m, n+1);
  _2D_INIT_MEM(G1,   m, n+1);
  _2D_INIT_MEM(G2,   m, n+1);
  _2D_INIT_MEM(G3,   m, n+1);
  _2D_INIT_MEM(G4,   m, n+1);

  double sigma = 0.0, speed_max = 0.0;  // speed_max denote the largest character speed at each
  // time step
  double half_tau, half_nu, mu, nu;  // tau is the length of the time step, nu = tau/h

  double h_S_max; // h/S_max, S_max is the maximum wave speed
  double time_c = 0.0; // the current time
  int nt = 1; // the number of times storing plotting data

  double UL, VL, PL, RHOL, SUL = 0.0, SVL = 0.0, SPL = 0.0, SRHOL = 0.0; // Left     boundary condition
  double UR, VR, PR, RHOR, SUR = 0.0, SVR = 0.0, SPR = 0.0, SRHOR = 0.0; // Right    boundary condition
  double UU, VU, PU, RHOU, SUU = 0.0, SVU = 0.0, SPU = 0.0, SRHOU = 0.0; // Upper    boundary condition
  double UD, VD, PD, RHOD, SUD = 0.0, SVD = 0.0, SPD = 0.0, SRHOD = 0.0; // Downside boundary condition

  for(j = 0; j < m; ++j)
    for(i = 0; i < n; ++i)
    {
      c = sqrt(gamma * CV->P[j][i] / CV->RHO[j][i]);
      sigma = fabs(c) + fabs(CV->U[j][i] + fabs(CV->V[j][i]));
      speed_max = ((speed_max < sigma) ? sigma : speed_max);
    }
  tau = (CFL * ((h_x<h_y)?h_x:h_y)) / speed_max;
  half_tau = tau * 0.5;
  half_nu = half_tau / h_x;
  mu = tau / h_y;


  for(i = 0; i < n; ++i)
  {
    for(j = 1; j < m; ++j)
    {
      s_rho[j][i] = alpha*(CV->RHO[j][i] - CV->RHO[j-1][i])/h_x;
        s_u[j][i] = alpha*(  CV->U[j][i] -   CV->U[j-1][i])/h_x;
        s_v[j][i] = alpha*(  CV->V[j][i] -   CV->V[j-1][i])/h_x;
        s_p[j][i] = alpha*(  CV->P[j][i] -   CV->P[j-1][i])/h_x;
    }
      s_rho[0][i] = alpha*(CV->RHO[1][i] - CV->RHO[0][i])/h_x;
        s_u[0][i] = alpha*(  CV->U[1][i] -   CV->U[0][i])/h_x;
        s_v[0][i] = alpha*(  CV->V[1][i] -   CV->V[0][i])/h_x;
        s_p[0][i] = alpha*(  CV->P[1][i] -   CV->P[0][i])/h_x;
  }
  for(j = 0; j < m; ++j)
  {
    for(i = 1; i < n; ++i)
    {
      t_rho[j][i] = alpha*(CV->RHO[j][i] - CV->RHO[j][i-1])/h_y;
        t_u[j][i] = alpha*(  CV->U[j][i] -   CV->U[j][i-1])/h_y;
        t_v[j][i] = alpha*(  CV->V[j][i] -   CV->V[j][i-1])/h_y;
        t_p[j][i] = alpha*(  CV->P[j][i] -   CV->P[j][i-1])/h_y;
    }
    t_rho[j][0] = alpha*(CV->RHO[j][1] - CV->RHO[j][0])/h_y;
      t_u[j][0] = alpha*(  CV->U[j][1] -   CV->U[j][0])/h_y;
      t_v[j][0] = alpha*(  CV->V[j][1] -   CV->V[j][0])/h_y;
      t_p[j][0] = alpha*(  CV->P[j][1] -   CV->P[j][0])/h_y;
  }


//------------THE MAIN LOOP-------------
  for(k = 1; k <= N; ++k)
  {
    /* evaluate f and a at some grid points for the iteration
     * and evaluate the character speed to decide the length
     * of the time step by (tau * speed_max)/h = CFL
     */
    h_S_max = INFINITY; // h/S_max = INFINITY
    tic = clock();

    slope_limiter_x(h_x, m, n, s_rho, (CV+nt-1)->RHO, s_u, (CV+nt-1)->U, s_v, (CV+nt-1)->V, s_p, (CV+nt-1)->P, alpha);
    flux_generator_x(h_x, m, n, half_tau, gamma, s_rho, (CV+nt-1)->RHO, s_u, (CV+nt-1)->U, s_v, (CV+nt-1)->V, s_p, (CV+nt-1)->P, F1, F2, F3, F4, rhox, ux, vx, px, eps);
//===============THE CORE ITERATION=================
    for(i = 0; i < n; ++i)
      for(j = 0; j < m; ++j)
      { /*
	 *  j-1          j          j+1
	 * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	 *   o-----X-----o-----X-----o-----X--...
	 */
	(CV+nt)->RHO[j][i] = (CV+nt-1)->RHO[j][i]       - half_nu*(F1[j+1][i]-F1[j][i]);
	mom_x = (CV+nt-1)->RHO[j][i]*(CV+nt-1)->U[j][i] - half_nu*(F2[j+1][i]-F2[j][i]);
	mom_y = (CV+nt-1)->RHO[j][i]*(CV+nt-1)->V[j][i] - half_nu*(F3[j+1][i]-F3[j][i]);
	ene = (CV+nt-1)->P[j][i]/(gamma-1.0) + 0.5*(CV+nt-1)->RHO[j][i]*((CV+nt-1)->U[j][i]*(CV+nt-1)->U[j][i]);
	ene = ene - half_nu*(F4[j+1][i]-F4[j][i]);

	(CV+nt)->U[j][i] = mom_x / (CV+nt)->RHO[j][i];
	(CV+nt)->V[j][i] = mom_y / (CV+nt)->RHO[j][i];
	(CV+nt)->P[j][i] = (ene - 0.5*mom_x*(CV+nt)->U[j][i])*(gamma-1.0);

	s_rho[j][i] = (rhox[j+1][i] - rhox[j][i])/h_x;
	  s_u[j][i] = (  ux[j+1][i] -   ux[j][i])/h_x;
	  s_v[j][i] = (  vx[j+1][i] -   vx[j][i])/h_x;
	  s_p[j][i] = (  px[j+1][i] -   px[j][i])/h_x;
      }
    //*/
//==================================================

    slope_limiter_y(h_y, m, n, t_rho, (CV+nt)->RHO, t_u, (CV+nt)->U, t_v, (CV+nt)->V, t_p, (CV+nt)->P, alpha);
    flux_generator_y(h_y, m, n, tau, gamma, t_rho, (CV+nt)->RHO, t_u, (CV+nt)->U, t_v, (CV+nt)->V, t_p, (CV+nt)->P, G1, G2, G3, G4, rhoy, uy, vy, py, eps);
    //printf("Hello world!\n");
//===============THE CORE ITERATION=================
    for(j = 0; j < m; ++j)
      for(i = 0; i < n; ++i)
      { /*
	 *  j-1          j          j+1
	 * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	 *   o-----X-----o-----X-----o-----X--...
	 */
	mom_x = (CV+nt)->RHO[j][i]*(CV+nt)->U[j][i] - mu*(G2[j][i+1]-G2[j][i]);
	mom_y = (CV+nt)->RHO[j][i]*(CV+nt)->V[j][i] - mu*(G3[j][i+1]-G3[j][i]);
	ene = (CV+nt)->P[j][i]/(gamma-1.0) + 0.5*(CV+nt)->RHO[j][i]*(CV+nt)->V[j][i]*(CV+nt)->V[j][i];
	ene = ene - mu*(G4[j][i+1]-G4[j][i]);
    (CV+nt)->RHO[j][i] = (CV+nt)->RHO[j][i]     - mu*(G1[j][i+1]-G1[j][i]);

	(CV+nt)->U[j][i] = mom_x / (CV+nt)->RHO[j][i];
	(CV+nt)->V[j][i] = mom_y / (CV+nt)->RHO[j][i];
	(CV+nt)->P[j][i] = (ene - 0.5*mom_y*(CV+nt)->V[j][i])*(gamma-1.0);

	s_rho[j][i] = (rhoy[j][i+1] - rhoy[j][i])/h_y;
	  s_u[j][i] = (  uy[j][i+1] -   uy[j][i])/h_y;
	  s_v[j][i] = (  vy[j][i+1] -   vy[j][i])/h_y;
	  s_p[j][i] = (  py[j][i+1] -   py[j][i])/h_y;
      }
//==================================================

    slope_limiter_x(h_x, m, n, s_rho, (CV+nt)->RHO, s_u, (CV+nt)->U, s_v, (CV+nt)->V, s_p, (CV+nt)->P, alpha);
    flux_generator_x(h_x, m, n, half_tau, gamma, s_rho, (CV+nt)->RHO, s_u, (CV+nt)->U, s_v, (CV+nt)->V, s_p, (CV+nt)->P, F1, F2, F3, F4, rhox, ux, vx, px, eps);
//===============THE CORE ITERATION=================
    for(i = 0; i < n; ++i)
      for(j = 0; j < m; ++j)
      { /*
	 *  j-1          j          j+1
	 * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	 *   o-----X-----o-----X-----o-----X--...
	 */
	mom_x = (CV+nt)->RHO[j][i]*(CV+nt)->U[j][i] - half_nu*(F2[j+1][i]-F2[j][i]);
	mom_y = (CV+nt)->RHO[j][i]*(CV+nt)->V[j][i] - half_nu*(F3[j+1][i]-F3[j][i]);
	ene = (CV+nt)->P[j][i]/(gamma-1.0) + 0.5*(CV+nt)->RHO[j][i]*((CV+nt)->U[j][i]*(CV+nt)->U[j][i]);
	ene = ene - half_nu*(F4[j+1][i]-F4[j][i]);
    (CV+nt)->RHO[j][i] = (CV+nt)->RHO[j][i]     - half_nu*(F1[j+1][i]-F1[j][i]);

	(CV+nt)->U[j][i] = mom_x / (CV+nt)->RHO[j][i];
	(CV+nt)->V[j][i] = mom_y / (CV+nt)->RHO[j][i];
	(CV+nt)->P[j][i] = (ene - 0.5*mom_x*(CV+nt)->U[j][i])*(gamma-1.0);

	s_rho[j][i] = (rhox[j+1][i] - rhox[j][i])/h_x;
	  s_u[j][i] = (  ux[j+1][i] -   ux[j][i])/h_x;
	  s_v[j][i] = (  vx[j+1][i] -   vx[j][i])/h_x;
	  s_p[j][i] = (  px[j+1][i] -   px[j][i])/h_x;
      }
//==================================================

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
    for(i = 0; i < n; ++i)
	{
	    (CV+nt-1)->RHO[j][i] = (CV+nt)->RHO[j][i];
	    (CV+nt-1)->U[j][i]   =   (CV+nt)->U[j][i];
	    (CV+nt-1)->V[j][i]   =   (CV+nt)->V[j][i];
	    (CV+nt-1)->E[j][i]   =   (CV+nt)->E[j][i];  
	    (CV+nt-1)->P[j][i]   =   (CV+nt)->P[j][i];
	}
  }

  printf("\nTime is up at time step %d.\n", k);
  printf("The cost of CPU time for 1D-GRP Eulerian scheme for this problem is %g seconds.\n", cpu_time_sum);
//------------END OF THE MAIN LOOP-------------

return_NULL:
  for(j = 0; j < m+1; ++j)
  {
    free(rhox[j]); free(ux[j]); free(vx[j]); free(px[j]);
    free(F1[j]);   free(F2[j]); free(F3[j]); free(F4[j]);
  }
  for(j = 0; j < m; ++j)
  {
    free(rhoy[j]); free(uy[j]); free(vy[j]); free(py[j]);
    free(G1[j]);   free(G2[j]); free(G3[j]); free(G4[j]);
    free(s_rho[j]); free(s_u[j]); free(s_v[j]); free(s_p[j]);
    free(t_rho[j]); free(t_u[j]); free(t_v[j]); free(t_p[j]);
  }
    free(rhox); free(ux); free(vx); free(px);
    free(F1);   free(F2); free(F3); free(F4);
    free(rhoy); free(uy); free(vy); free(py); 
    free(G1);   free(G2); free(G3); free(G4);
    free(s_rho); free(s_u); free(s_v); free(s_p);
    free(t_rho); free(t_u); free(t_v); free(t_p);
}
