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
#include "../include/flux_calc.h"
#include "../include/inter_process.h"
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

#define _1D_BC_INIT_MEM(v, M)				\
    do {						\
	v = (struct b_f_var *)calloc((M), sizeof(struct b_f_var));	\
	if(v == NULL)					\
	    {						\
		printf("NOT enough memory! %s\n", #v);	\
		goto return_NULL;			\
	    }						\
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
     * i is a frequently used index for y-spatial variables.
     * j is a frequently used index for x-spatial variables.
     * k is a frequently used index for the time step.
     */
  int i, j, k;

  clock_t tic, toc;
  double cpu_time_sum = 0.0;

  double const t_all     = config[1];        // the total time
  double const eps       = config[4];        // the largest value could be seen as zero
  int    const N         = (int)(config[5]); // the maximum number of time steps
  double const gamma     = config[6];        // the constant of the perfect gas
  double const CFL       = config[7];        // the CFL number
  double const h_x       = config[10];       // the length of the initial x-spatial grids
  double const h_y       = config[11];       // the length of the initial y-spatial grids
  double       tau       = config[16];       // the length of the time step
  _Bool  const dim_split = (_Bool)config[33];// Dimensional splitting?

  _Bool find_bound_x = false, find_bound_y = false;

  double mom_x, mom_y, ene;
  double c; // the speeds of sound
  
  // Left/Right/Upper/Downside boundary condition
  struct b_f_var * bfv_L = NULL, * bfv_R = NULL, * bfv_U = NULL, * bfv_D = NULL; 
  // the slopes of variable values.
  _2D_INIT_MEM(CV->s_rho, m, n); _2D_INIT_MEM(CV->t_rho, m, n);
  _2D_INIT_MEM(CV->s_u,   m, n); _2D_INIT_MEM(CV->t_u,   m, n);
  _2D_INIT_MEM(CV->s_v,   m, n); _2D_INIT_MEM(CV->t_v,   m, n);
  _2D_INIT_MEM(CV->s_p,   m, n); _2D_INIT_MEM(CV->t_p,   m, n);
  // the variable values at (x_{j-1/2}, t_{n+1}).
  _2D_INIT_MEM(CV->rhoIx, m+1, n);
  _2D_INIT_MEM(CV->uIx,   m+1, n);
  _2D_INIT_MEM(CV->vIx,   m+1, n);
  _2D_INIT_MEM(CV->pIx,   m+1, n);
  _2D_INIT_MEM(CV->F_rho, m+1, n);
  _2D_INIT_MEM(CV->F_u,   m+1, n);
  _2D_INIT_MEM(CV->F_v,   m+1, n);
  _2D_INIT_MEM(CV->F_e,   m+1, n); 
  // the variable values at (y_{j-1/2}, t_{n+1}).
  _2D_INIT_MEM(CV->rhoIy, m, n+1);
  _2D_INIT_MEM(CV->uIy,   m, n+1);
  _2D_INIT_MEM(CV->vIy,   m, n+1);
  _2D_INIT_MEM(CV->pIy,   m, n+1);
  _2D_INIT_MEM(CV->G_rho, m, n+1);
  _2D_INIT_MEM(CV->G_u,   m, n+1);
  _2D_INIT_MEM(CV->G_v,   m, n+1);
  _2D_INIT_MEM(CV->G_e,   m, n+1);
  // boundary condition
  _1D_BC_INIT_MEM(bfv_L, n); _1D_BC_INIT_MEM(bfv_R, n);
  _1D_BC_INIT_MEM(bfv_D, m); _1D_BC_INIT_MEM(bfv_U, m);
  
  double half_tau, half_nu, mu, nu;  // nu = tau/h_x, mu = tau/h_y.

  double h_S_max, sigma; // h/S_max, S_max is the maximum character speed, sigma is the character speed
  double time_c = 0.0; // the current time
  int nt = 1; // the number of times storing plotting data

//------------THE MAIN LOOP-------------
  for(k = 1; k <= N; ++k)
  {
    /* evaluate f and a at some grid points for the iteration
     * and evaluate the character speed to decide the length
     * of the time step by (tau * speed_max)/h = CFL
     */
      h_S_max = INFINITY; // h/S_max = INFINITY
      tic = clock();

    for(j = 0; j < m; ++j)
	for(i = 0; i < n; ++i)
	    {
		c = sqrt(gamma * CV->P[j][i] / CV->RHO[j][i]);
		sigma = fabs(c) + fabs(CV->U[j][i]) + fabs(CV->V[j][i]);
		h_S_max = fmin(h_S_max, fmin(h_x,h_y) / sigma);
	    }
    // If no total time, use fixed tau and time step N.
    if (isfinite(t_all) || !isfinite(config[16]) || config[16] <= 0.0)
	{
	    tau = CFL * h_S_max;
	    if ((time_c + tau) > (t_all - eps))
		tau = t_all - time_c;
	}
    mu = tau / h_y;
    nu = tau / h_x;
    if(dim_split)
	{
	    half_tau = tau * 0.5;
	    half_nu = half_tau / h_x;
	}
    else // NO dimensional splitting!
	{
	    half_tau = tau;
	    half_nu  = nu;
	}

    find_bound_x = bound_cond_slope_limiter_x(m, n, nt-1, CV, bfv_L, bfv_R, find_bound_x);
    if(find_bound_x)
        goto return_NULL;

    flux_generator_x(m, n, nt-1, half_tau, CV, bfv_L, bfv_R);

//===============THE CORE ITERATION=================
    for(i = 0; i < n; ++i)
      for(j = 0; j < m; ++j)
      { /*
	 *  j-1          j          j+1
	 * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	 *   o-----X-----o-----X-----o-----X--...
	 */
	  (CV+nt)->RHO[j][i] = (CV+nt-1)->RHO[j][i]       - half_nu*(CV->F_rho[j+1][i]-CV->F_rho[j][i]);
	  mom_x = (CV+nt-1)->RHO[j][i]*(CV+nt-1)->U[j][i] - half_nu*(CV->F_u[j+1][i]  -CV->F_u[j][i]);
	  mom_y = (CV+nt-1)->RHO[j][i]*(CV+nt-1)->V[j][i] - half_nu*(CV->F_v[j+1][i]  -CV->F_v[j][i]);
	  ene   = (CV+nt-1)->RHO[j][i]*(CV+nt-1)->E[j][i] - half_nu*(CV->F_e[j+1][i]  -CV->F_e[j][i]);
	  
	  (CV+nt)->U[j][i] = mom_x / (CV+nt)->RHO[j][i];
	  (CV+nt)->V[j][i] = mom_y / (CV+nt)->RHO[j][i];
	  (CV+nt)->E[j][i] = ene   / (CV+nt)->RHO[j][i];
	  (CV+nt)->P[j][i] = (ene - 0.5*mom_x*(CV+nt)->U[j][i])*(gamma-1.0);
	  
	  CV->s_rho[j][i] = (CV->rhoIx[j+1][i] - CV->rhoIx[j][i])/h_x;
	  CV->s_u[j][i]   = (  CV->uIx[j+1][i] -   CV->uIx[j][i])/h_x;
	  CV->s_v[j][i]   = (  CV->vIx[j+1][i] -   CV->vIx[j][i])/h_x;
	  CV->s_p[j][i]   = (  CV->pIx[j+1][i] -   CV->pIx[j][i])/h_x;
      }

//==================================================

    find_bound_y = bound_cond_slope_limiter_y(m, n, nt, CV, bfv_D, bfv_U, find_bound_y);
    if(find_bound_y)
        goto return_NULL;
    flux_generator_y(m, n, nt, tau, CV, bfv_D, bfv_U);

//===============THE CORE ITERATION=================
    for(j = 0; j < m; ++j)
      for(i = 0; i < n; ++i)
      { /*
	 *  j-1          j          j+1
	 * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	 *   o-----X-----o-----X-----o-----X--...
	 */
	  mom_x = (CV+nt)->RHO[j][i]*(CV+nt)->U[j][i] - mu*(CV->G_u[j][i+1]  -CV->G_u[j][i]);
	  mom_y = (CV+nt)->RHO[j][i]*(CV+nt)->V[j][i] - mu*(CV->G_v[j][i+1]  -CV->G_v[j][i]);
	  ene   = (CV+nt)->RHO[j][i]*(CV+nt)->E[j][i] - mu*(CV->G_e[j][i+1]  -CV->G_e[j][i]);
	  (CV+nt)->RHO[j][i] = (CV+nt)->RHO[j][i]     - mu*(CV->G_rho[j][i+1]-CV->G_rho[j][i]);
	  
	  (CV+nt)->U[j][i] = mom_x / (CV+nt)->RHO[j][i];
	  (CV+nt)->V[j][i] = mom_y / (CV+nt)->RHO[j][i];
	  (CV+nt)->E[j][i] = ene   / (CV+nt)->RHO[j][i];
	  (CV+nt)->P[j][i] = (ene - 0.5*mom_y*(CV+nt)->V[j][i])*(gamma-1.0);
	  
	  CV->t_rho[j][i] = (CV->rhoIy[j][i+1] - CV->rhoIy[j][i])/h_y;
	  CV->t_u[j][i]   = (  CV->uIy[j][i+1] -   CV->uIy[j][i])/h_y;
	  CV->t_v[j][i]   = (  CV->vIy[j][i+1] -   CV->vIy[j][i])/h_y;
	  CV->t_p[j][i]   = (  CV->pIy[j][i+1] -   CV->pIy[j][i])/h_y;
      }
//==================================================

    if(dim_split)
    {
    bound_cond_slope_limiter_x(m, n, nt-1, CV, bfv_L, bfv_R, find_bound_x);
    flux_generator_x(m, n, nt-1, half_tau, CV, bfv_L, bfv_R);

//===============THE CORE ITERATION=================
    for(i = 0; i < n; ++i)
      for(j = 0; j < m; ++j)
      { /*
	 *  j-1          j          j+1
	 * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	 *   o-----X-----o-----X-----o-----X--...
	 */
	  mom_x = (CV+nt)->RHO[j][i]*(CV+nt)->U[j][i] - half_nu*(CV->F_u[j+1][i]  -CV->F_u[j][i]);
	  mom_y = (CV+nt)->RHO[j][i]*(CV+nt)->V[j][i] - half_nu*(CV->F_v[j+1][i]  -CV->F_v[j][i]);
	  ene   = (CV+nt)->RHO[j][i]*(CV+nt)->E[j][i] - half_nu*(CV->F_e[j+1][i]  -CV->F_e[j][i]);
	  (CV+nt)->RHO[j][i] = (CV+nt)->RHO[j][i]     - half_nu*(CV->F_rho[j+1][i]-CV->F_rho[j][i]);
	  
	  (CV+nt)->U[j][i] = mom_x / (CV+nt)->RHO[j][i];
	  (CV+nt)->V[j][i] = mom_y / (CV+nt)->RHO[j][i];
	  (CV+nt)->E[j][i] = ene   / (CV+nt)->RHO[j][i];
	  (CV+nt)->P[j][i] = (ene - 0.5*mom_x*(CV+nt)->U[j][i])*(gamma-1.0);
	  
	  CV->s_rho[j][i] = (CV->rhoIx[j+1][i] - CV->rhoIx[j][i])/h_x;
	  CV->s_u[j][i]   = (  CV->uIx[j+1][i] -   CV->uIx[j][i])/h_x;
	  CV->s_v[j][i]   = (  CV->vIx[j+1][i] -   CV->vIx[j][i])/h_x;
	  CV->s_p[j][i]   = (  CV->pIx[j+1][i] -   CV->pIx[j][i])/h_x;
      }
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
    free(CV->F_rho[j]); free(CV->F_u[j]); free(CV->F_v[j]); free(CV->F_e[j]);
    free(CV->rhoIx[j]); free(CV->uIx[j]); free(CV->vIx[j]); free(CV->pIx[j]);
    CV->F_rho[j]= NULL; CV->F_u[j]= NULL; CV->F_v[j]= NULL; CV->F_e[j]= NULL;
    CV->rhoIx[j]= NULL; CV->uIx[j]= NULL; CV->vIx[j]= NULL; CV->pIx[j]= NULL;
  }
  for(j = 0; j < m; ++j)
  {
    free(CV->G_rho[j]); free(CV->G_u[j]); free(CV->G_v[j]); free(CV->G_e[j]);
    free(CV->rhoIy[j]); free(CV->uIy[j]); free(CV->vIy[j]); free(CV->pIy[j]);
    free(CV->s_rho[j]); free(CV->s_u[j]); free(CV->s_v[j]); free(CV->s_p[j]);
    free(CV->t_rho[j]); free(CV->t_u[j]); free(CV->t_v[j]); free(CV->t_p[j]);

    CV->G_rho[j]= NULL; CV->G_u[j]= NULL; CV->G_v[j]= NULL; CV->G_e[j]= NULL; 
    CV->rhoIy[j]= NULL; CV->uIy[j]= NULL; CV->vIy[j]= NULL; CV->pIy[j]= NULL; 
    CV->s_rho[j]= NULL; CV->s_u[j]= NULL; CV->s_v[j]= NULL; CV->s_p[j]= NULL; 
    CV->t_rho[j]= NULL; CV->t_u[j]= NULL; CV->t_v[j]= NULL; CV->t_p[j]= NULL; 
  }
    free(CV->F_rho); free(CV->F_u); free(CV->F_v); free(CV->F_e);
    free(CV->rhoIx); free(CV->uIx); free(CV->vIx); free(CV->pIx);
    free(CV->G_rho); free(CV->G_u); free(CV->G_v); free(CV->G_e);
    free(CV->rhoIy); free(CV->uIy); free(CV->vIy); free(CV->pIy); 
    free(CV->s_rho); free(CV->s_u); free(CV->s_v); free(CV->s_p);
    free(CV->t_rho); free(CV->t_u); free(CV->t_v); free(CV->t_p);
    free(bfv_L); free(bfv_R);
    free(bfv_D); free(bfv_U);
    
    CV->F_rho= NULL; CV->F_u= NULL; CV->F_v= NULL; CV->F_e= NULL;
    CV->rhoIx= NULL; CV->uIx= NULL; CV->vIx= NULL; CV->pIx= NULL;
    CV->G_rho= NULL; CV->G_u= NULL; CV->G_v= NULL; CV->G_e= NULL;
    CV->rhoIy= NULL; CV->uIy= NULL; CV->vIy= NULL; CV->pIy= NULL; 
    CV->s_rho= NULL; CV->s_u= NULL; CV->s_v= NULL; CV->s_p= NULL;
    CV->t_rho= NULL; CV->t_u= NULL; CV->t_v= NULL; CV->t_p= NULL;
    bfv_L= NULL; bfv_R= NULL;
    bfv_D= NULL; bfv_U= NULL;
}
