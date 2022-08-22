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

  double const t_all   = config[1];        // the total time
  double const eps     = config[4];        // the largest value could be seen as zero
  int    const N       = (int)(config[5]); // the maximum number of time steps
  double const gamma   = config[6];        // the constant of the perfect gas
  double const CFL     = config[7];        // the CFL number
  double const h_x     = config[10];       // the length of the initial x-spatial grids
  double const h_y     = config[11];       // the length of the initial y-spatial grids
  double       tau     = config[16];       // the length of the time step
  double const alpha   = config[41];       // the paramater in slope limiters.

  _Bool find_bound_x = false, find_bound_y = false;

  double mom_x, mom_y, ene;
  double c; // the speeds of sound
  struct i_f_var ifv_L, ifv_R;
  /*
   * dire: the temporal derivative of fluid variables.
   *       \frac{\partial [rho, u, v, p]}{\partial t}
   * mid:  the Riemann solutions.
   *       [rho_star, u_star, v_star, p_star]
   */
  double dire[4], mid[4];

  // the numerical flux at (x_{j-1/2}, t_{n}).
  double ** F1, ** F2, ** F3, ** F4;
  // the numerical flux at (y_{j-1/2}, t_{n}).
  double ** G1, ** G2, ** G3, ** G4;
  struct b_f_var * bfv_L, * bfv_R, * bfv_U, * bfv_D; // Left/Right/Upper/Downside boundary condition
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
  _2D_INIT_MEM(F1, m+1, n); _2D_INIT_MEM(F2, m+1, n);
  _2D_INIT_MEM(F3, m+1, n); _2D_INIT_MEM(F4, m+1, n); 
  // the variable values at (y_{j-1/2}, t_{n+1}).
  _2D_INIT_MEM(CV->rhoIy, m, n+1);
  _2D_INIT_MEM(CV->uIy,   m, n+1);
  _2D_INIT_MEM(CV->vIy,   m, n+1);
  _2D_INIT_MEM(CV->pIy,   m, n+1);
  _2D_INIT_MEM(G1, m, n+1); _2D_INIT_MEM(G2, m, n+1);
  _2D_INIT_MEM(G3, m, n+1); _2D_INIT_MEM(G4, m, n+1);
  // boundary condition
  _1D_BC_INIT_MEM(bfv_L, n); _1D_BC_INIT_MEM(bfv_R, n);
  _1D_BC_INIT_MEM(bfv_U, m); _1D_BC_INIT_MEM(bfv_D, m);
  
  double half_tau, half_nu, mu, nu;  // nu = tau/h_x, mu = tau/h_y.

  double h_S_max, sigma; // h/S_max, S_max is the maximum character speed, sigma is the character speed
  double time_c = 0.0; // the current time
  int nt = 1; // the number of times storing plotting data

  for(i = 0; i < n; ++i)
  {
    for(j = 1; j < m; ++j)
    {
	CV->s_rho[j][i] = alpha*(CV->RHO[j][i] - CV->RHO[j-1][i])/h_x;
        CV->s_u[j][i] = alpha*(  CV->U[j][i] -   CV->U[j-1][i])/h_x;
        CV->s_v[j][i] = alpha*(  CV->V[j][i] -   CV->V[j-1][i])/h_x;
        CV->s_p[j][i] = alpha*(  CV->P[j][i] -   CV->P[j-1][i])/h_x;
    }
    CV->s_rho[0][i] = alpha*(CV->RHO[1][i] - CV->RHO[0][i])/h_x;
    CV->s_u[0][i] = alpha*(  CV->U[1][i] -   CV->U[0][i])/h_x;
    CV->s_v[0][i] = alpha*(  CV->V[1][i] -   CV->V[0][i])/h_x;
    CV->s_p[0][i] = alpha*(  CV->P[1][i] -   CV->P[0][i])/h_x;
  }
  for(j = 0; j < m; ++j)
  {
    for(i = 1; i < n; ++i)
    {
	CV->t_rho[j][i] = alpha*(CV->RHO[j][i] - CV->RHO[j][i-1])/h_y;
        CV->t_u[j][i] = alpha*(  CV->U[j][i] -   CV->U[j][i-1])/h_y;
        CV->t_v[j][i] = alpha*(  CV->V[j][i] -   CV->V[j][i-1])/h_y;
        CV->t_p[j][i] = alpha*(  CV->P[j][i] -   CV->P[j][i-1])/h_y;
    }
    CV->t_rho[j][0] = alpha*(CV->RHO[j][1] - CV->RHO[j][0])/h_y;
    CV->t_u[j][0] = alpha*(  CV->U[j][1] -   CV->U[j][0])/h_y;
    CV->t_v[j][0] = alpha*(  CV->V[j][1] -   CV->V[j][0])/h_y;
    CV->t_p[j][0] = alpha*(  CV->P[j][1] -   CV->P[j][0])/h_y;
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
    half_tau = tau * 0.5;
    half_nu = half_tau / h_x;
    mu = tau / h_y;

    if(find_bound_x = bound_cond_slope_limiter_x(m, n, nt, CV, bfv_L, bfv_R, find_bound_x))
        goto return_NULL;

    flux_generator_x(m, n, half_tau, CV, F1, F2, F3, F4);

//===============THE CORE ITERATION=================
    for(i = 0; i < n; ++i)
      for(j = 0; j < m; ++j)
      { /*
	 *  j-1          j          j+1
	 * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	 *   o-----X-----o-----X-----o-----X--...
	 */
	  mom_x = (CV+nt-1)->RHO[j][i]*(CV+nt-1)->U[j][i] - half_nu*(F2[j+1][i]-F2[j][i]);
	  mom_y = (CV+nt-1)->RHO[j][i]*(CV+nt-1)->V[j][i] - half_nu*(F3[j+1][i]-F3[j][i]);
	  ene   = (CV+nt-1)->RHO[j][i]*(CV+nt-1)->E[j][i] - half_nu*(F4[j+1][i]-F4[j][i]);
	  (CV+nt)->RHO[j][i] = (CV+nt-1)->RHO[j][i]       - half_nu*(F1[j+1][i]-F1[j][i]);
	  
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

    if(find_bound_y = bound_cond_slope_limiter_y(m, n, nt, CV, bfv_U, bfv_D, find_bound_y))
        goto return_NULL;
    flux_generator_y(m, n, tau, CV, G1, G2, G3, G4);

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
	  ene   = (CV+nt)->RHO[j][i]*(CV+nt)->E[j][i] - mu*(G4[j][i+1]-G4[j][i]);
	  (CV+nt)->RHO[j][i] = (CV+nt)->RHO[j][i]     - mu*(G1[j][i+1]-G1[j][i]);
	  
	  (CV+nt)->U[j][i] = mom_x / (CV+nt)->RHO[j][i];
	  (CV+nt)->V[j][i] = mom_y / (CV+nt)->RHO[j][i];
	  (CV+nt)->E[j][i] = ene   / (CV+nt)->RHO[j][i];
	  (CV+nt)->P[j][i] = (ene - 0.5*mom_y*(CV+nt)->V[j][i])*(gamma-1.0);
	  
	  CV->s_rho[j][i] = (CV->rhoIy[j][i+1] - CV->rhoIy[j][i])/h_y;
	  CV->s_u[j][i]   = (  CV->uIy[j][i+1] -   CV->uIy[j][i])/h_y;
	  CV->s_v[j][i]   = (  CV->vIy[j][i+1] -   CV->vIy[j][i])/h_y;
	  CV->s_p[j][i]   = (  CV->pIy[j][i+1] -   CV->pIy[j][i])/h_y;
      }
//==================================================

    bound_cond_slope_limiter_x(m, n, nt, CV, bfv_L, bfv_R, find_bound_x);
    flux_generator_x(m, n, half_tau, CV, F1, F2, F3, F4);

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
	  ene   = (CV+nt)->RHO[j][i]*(CV+nt)->E[j][i] - half_nu*(F4[j+1][i]-F4[j][i]);
	  (CV+nt)->RHO[j][i] = (CV+nt)->RHO[j][i]     - half_nu*(F1[j+1][i]-F1[j][i]);
	  
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
    free(F1[j]); free(F2[j]); free(F3[j]); free(F4[j]);
    free(CV->rhoIx[j]); free(CV->uIx[j]); free(CV->vIx[j]); free(CV->pIx[j]);
  }
  for(j = 0; j < m; ++j)
  {
    free(G1[j]); free(G2[j]); free(G3[j]); free(G4[j]);
    free(CV->rhoIy[j]); free(CV->uIy[j]); free(CV->vIy[j]); free(CV->pIy[j]);
    free(CV->s_rho[j]); free(CV->s_u[j]); free(CV->s_v[j]); free(CV->s_p[j]);
    free(CV->t_rho[j]); free(CV->t_u[j]); free(CV->t_v[j]); free(CV->t_p[j]);
  }
    free(F1);   free(F2); free(F3); free(F4);
    free(G1);   free(G2); free(G3); free(G4);
    free(CV->rhoIx); free(CV->uIx); free(CV->vIx); free(CV->pIx);
    free(CV->rhoIy); free(CV->uIy); free(CV->vIy); free(CV->pIy); 
    free(CV->s_rho); free(CV->s_u); free(CV->s_v); free(CV->s_p);
    free(CV->t_rho); free(CV->t_u); free(CV->t_v); free(CV->t_p);
}
