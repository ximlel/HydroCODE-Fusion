/**
 * @file  grp_solver_2D_EUL_source.c
 * @brief This is an Eulerian GRP scheme to solve 2-D Euler equations without dimension splitting.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#ifdef _OPENACC
#include <omp.h>
#include <openacc.h>
#elif defined _OPENMP
#include <omp.h>
#endif

#include "../include/var_struc.h"
#include "../include/riemann_solver.h"
#include "../include/flux_calc.h"
#include "../include/inter_process.h"
#include "../include/tools.h"


/**
 * @brief M*N memory allocations to the variable 'v' in the structure cell_var_stru.
 */
#define INIT_MEM_2D(v, M, N)						\
    do {								\
	CV->v = (double **)malloc((M) * sizeof(double *));			\
	if(CV->v == NULL)							\
	    {								\
		printf("NOT enough memory! %s\n", #v);			\
		goto return_NULL;					\
	    }								\
	for(j = 0; j < (M); ++j)					\
	    {								\
		CV->v[j] = (double *)malloc((N) * sizeof(double));		\
		if(CV->v[j] == NULL)					\
		    {							\
			printf("NOT enough memory! %s[%d]\n", #v, j);	\
			goto return_NULL;				\
		    }							\
	    }								\
    } while (0)

/**
 * @brief M memory allocations to the structure variable b_f_var 'bfv'.
 */
#define BC_INIT_MEM_1D(bfv, M)				\
    do {						\
	bfv = (struct b_f_var *)calloc((M), sizeof(struct b_f_var));	\
	if(bfv == NULL)					\
	    {						\
		printf("NOT enough memory! %s\n", #bfv);	\
		goto return_NULL;			\
	    }						\
    } while (0)

/**
 * @brief This function use GRP scheme to solve 2-D Euler
 *        equations of motion on Eulerian coordinate without dimension splitting.
 * @param[in] m:          Number of the x-grids: n_x.
 * @param[in] n:          Number of the y-grids: n_y.
 * @param[in,out] CV:     Structure of cell variable data.
 * @param[out] cpu_time:  Array of the CPU time recording.
 * @param[out] time_plot: Array of the plotting time recording.
 */
void GRP_solver_2D_EUL_source(const int m, const int n, struct cell_var_stru * CV, double * cpu_time, int * N_plot, double time_plot[])
{
#ifdef _OPENMP
  printf("@@ Number of threads for OpenMP: %d\n", omp_get_max_threads());
#endif
#ifdef _OPENACC
  printf("@@ Number of CPU devices for OpenACC: %d\n", acc_get_num_devices(acc_device_host));
  printf("@@ Number of GPU devices for OpenACC: %d\n", acc_get_num_devices(acc_device_not_host));
#endif
    /* 
     * i is a frequently used index for y-spatial variables.
     * j is a frequently used index for x-spatial variables.
     * k is a frequently used index for the time step.
     */
  int i, j, k = 0;

  double tic, toc;
  double cpu_time_sum = 0.0;

  double const t_all     = config[1];      // the total time
  double const eps       = config[4];      // the largest value could be seen as zero
  int    const N         = (int)config[5]; // the maximum number of time steps
  double const gamma     = config[6];      // the constant of the perfect gas
  double const CFL       = config[7];      // the CFL number
  double const h_x       = config[10];     // the length of the initial x-spatial grids
  double const h_y       = config[11];     // the length of the initial y-spatial grids
  double       tau       = config[16];     // the length of the time step

  _Bool find_bound_x = false, find_bound_y = false;
  int flux_err;

  double mom_x, mom_y, ene;
  double c; // the speeds of sound

  double mu, nu;  // nu = tau/h_x, mu = tau/h_y.
  double h_S_max, sigma; // h/S_max, S_max is the maximum character speed, sigma is the character speed
  double time_c = 0.0; // the current time
  _Bool stop_t = false;
  int nt = 0; // the number of times storing plotting data
  
  // Left/Right/Upper/Downside boundary condition
  struct b_f_var * bfv_L = NULL, * bfv_R = NULL, * bfv_U = NULL, * bfv_D = NULL;
  // the slopes of variable values.
  INIT_MEM_2D(s_rho, m, n); INIT_MEM_2D(t_rho, m, n);
  INIT_MEM_2D(s_u,   m, n); INIT_MEM_2D(t_u,   m, n);
  INIT_MEM_2D(s_v,   m, n); INIT_MEM_2D(t_v,   m, n);
  INIT_MEM_2D(s_p,   m, n); INIT_MEM_2D(t_p,   m, n);
  // the variable values at (x_{j-1/2}, t_{n+1}).
  INIT_MEM_2D(rhoIx, m+1, n);
  INIT_MEM_2D(uIx,   m+1, n);
  INIT_MEM_2D(vIx,   m+1, n);
  INIT_MEM_2D(pIx,   m+1, n);
  INIT_MEM_2D(F_rho, m+1, n);
  INIT_MEM_2D(F_u,   m+1, n);
  INIT_MEM_2D(F_v,   m+1, n);
  INIT_MEM_2D(F_e,   m+1, n); 
  // the variable values at (y_{j-1/2}, t_{n+1}).
  INIT_MEM_2D(rhoIy, m, n+1);
  INIT_MEM_2D(uIy,   m, n+1);
  INIT_MEM_2D(vIy,   m, n+1);
  INIT_MEM_2D(pIy,   m, n+1);
  INIT_MEM_2D(G_rho, m, n+1);
  INIT_MEM_2D(G_u,   m, n+1);
  INIT_MEM_2D(G_v,   m, n+1);
  INIT_MEM_2D(G_e,   m, n+1);
  // boundary condition
  BC_INIT_MEM_1D(bfv_L, n); BC_INIT_MEM_1D(bfv_R, n);
  BC_INIT_MEM_1D(bfv_D, m); BC_INIT_MEM_1D(bfv_U, m);

//------------THE MAIN LOOP-------------
  for(k = 1; k <= N; ++k)
  {
#ifdef _OPENMP
    tic = omp_get_wtime();
#elif defined _OPENACC
    tic = omp_get_wtime();
#else
    tic = (double)clock() / (double)CLOCKS_PER_SEC;
#endif
    if (time_c >= time_plot[nt] && nt < (*N_plot-1))
	{
	    for(j = 0; j < m; ++j)
		for(i = 0; i < n; ++i)
		    {
			CV[nt+1].RHO[j][i] = CV[nt].RHO[j][i];
			CV[nt+1].U[j][i]   =   CV[nt].U[j][i];
			CV[nt+1].V[j][i]   =   CV[nt].V[j][i];
			CV[nt+1].E[j][i]   =   CV[nt].E[j][i];  
			CV[nt+1].P[j][i]   =   CV[nt].P[j][i];
		    }
	}

    /* evaluate f and a at some grid points for the iteration
     * and evaluate the character speed to decide the length
     * of the time step by (tau * speed_max)/h = CFL
     */
    h_S_max = INFINITY; // h/S_max = INFINITY

    for(j = 0; j < m; ++j)
	for(i = 0; i < n; ++i)
	    {
		c = sqrt(gamma * CV->P[j][i] / CV->RHO[j][i]);
		sigma = fabs(c) + fabs(CV->U[j][i]) + fabs(CV->V[j][i]);
		h_S_max = fmin(h_S_max, fmin(h_x,h_y) / sigma);
	    }
    // If no total time, use fixed tau and time step N.
    if(isfinite(t_all) || !isfinite(config[16]) || config[16] <= 0.0)
	{
	    tau = CFL * h_S_max;
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
    nu = tau / h_x;
    mu = tau / h_y;

    find_bound_x = bound_cond_slope_limiter_x(m, n, nt, CV, bfv_L, bfv_R, bfv_D, bfv_U, find_bound_x, true, time_c);
    if(!find_bound_x)
        goto return_NULL;
    find_bound_y = bound_cond_slope_limiter_y(m, n, nt, CV, bfv_L, bfv_R, bfv_D, bfv_U, find_bound_y, true, time_c);
    if(!find_bound_y)
        goto return_NULL;

    flux_err = flux_generator_x(m, n, nt, tau, CV, bfv_L, bfv_R, true);
    if(flux_err == 1)
        goto return_NULL;
    else if(flux_err == 2)
	stop_t = true;
    flux_err = flux_generator_y(m, n, nt, tau, CV, bfv_D, bfv_U, true);
    if(flux_err == 1)
        goto return_NULL;
    else if(flux_err == 2)
	stop_t = true;

//===============THE CORE ITERATION=================
#ifdef _OPENMP
#pragma omp parallel for  private(mom_x, mom_y, ene) collapse(2) schedule(dynamic, 8)
#elif defined _OPENACC
#pragma acc parallel loop private(mom_x, mom_y, ene) collapse(2) worker
#endif
    for(i = 0; i < n; ++i)
      for(j = 0; j < m; ++j)
      { /*
	 *  j-1          j          j+1
	 * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	 *   o-----X-----o-----X-----o-----X--...
	 */
	  mom_x = CV[nt].RHO[j][i]*CV[nt].U[j][i] - nu*(CV->F_u[j+1][i]  -CV->F_u[j][i])   - mu*(CV->G_u[j][i+1]  -CV->G_u[j][i]);
	  mom_y = CV[nt].RHO[j][i]*CV[nt].V[j][i] - nu*(CV->F_v[j+1][i]  -CV->F_v[j][i])   - mu*(CV->G_v[j][i+1]  -CV->G_v[j][i]);
	  ene   = CV[nt].RHO[j][i]*CV[nt].E[j][i] - nu*(CV->F_e[j+1][i]  -CV->F_e[j][i])   - mu*(CV->G_e[j][i+1]  -CV->G_e[j][i]);
	  CV[nt].RHO[j][i]   =   CV[nt].RHO[j][i] - nu*(CV->F_rho[j+1][i]-CV->F_rho[j][i]) - mu*(CV->G_rho[j][i+1]-CV->G_rho[j][i]);

	  CV[nt].U[j][i] = mom_x / CV[nt].RHO[j][i];
	  CV[nt].V[j][i] = mom_y / CV[nt].RHO[j][i];
	  CV[nt].E[j][i] = ene   / CV[nt].RHO[j][i];
	  CV[nt].P[j][i] = (ene - 0.5*mom_x*CV[nt].U[j][i] - 0.5*mom_y*CV[nt].V[j][i])*(gamma-1.0);
	  if(CV[nt].P[j][i] < eps || CV[nt].RHO[j][i] < eps)
	      {
		  printf("<0.0 error on [%d, %d, %d] (t_n, x, y) - Update\n", k, j, i);
		  stop_t = true;
	      }

	  CV->s_rho[j][i] = (CV->rhoIx[j+1][i] - CV->rhoIx[j][i])/h_x;
	  CV->s_u[j][i]   = (  CV->uIx[j+1][i] -   CV->uIx[j][i])/h_x;
	  CV->s_v[j][i]   = (  CV->vIx[j+1][i] -   CV->vIx[j][i])/h_x;
	  CV->s_p[j][i]   = (  CV->pIx[j+1][i] -   CV->pIx[j][i])/h_x;
	  CV->t_rho[j][i] = (CV->rhoIy[j][i+1] - CV->rhoIy[j][i])/h_y;
	  CV->t_u[j][i]   = (  CV->uIy[j][i+1] -   CV->uIy[j][i])/h_y;
	  CV->t_v[j][i]   = (  CV->vIy[j][i+1] -   CV->vIy[j][i])/h_y;
	  CV->t_p[j][i]   = (  CV->pIy[j][i+1] -   CV->pIy[j][i])/h_y;
      } // End of parallel region

//==================================================
    
    time_c += tau;
    if(isfinite(t_all))
        DispPro(time_c*100.0/t_all, k);
    else
        DispPro(k*100.0/N, k);
    if(stop_t || time_c > (t_all - eps) || !isfinite(time_c))
	break;

    //===========================Fixed variable location=======================

#ifdef _OPENMP
    toc = omp_get_wtime();
#elif defined _OPENACC
    toc = omp_get_wtime();
#else
    toc = (double)clock() / (double)CLOCKS_PER_SEC;
#endif
    cpu_time_sum += toc - tic;
    cpu_time[nt]  = cpu_time_sum;
  }

  printf("\nTime is up at time step %d.\n", k);
  printf("The cost of CPU time for genuinely 2D-GRP Eulerian scheme without dimension splitting for this problem is %g seconds.\n", cpu_time_sum);
  //------------END OF THE MAIN LOOP-------------
  
return_NULL:
  config[5] = (double)k;
  *N_plot = nt+1;
  if(isfinite(time_c))
      time_plot[nt] = time_c;
  else if(isfinite(t_all))
      time_plot[nt] = t_all;
  else if(isfinite(tau))
      time_plot[nt] = k*tau;

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
