#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN32
#define ISNAN(a) _isnan((a))
#elif __linux__
#define ISNAN(a) isnan((a))
#endif

#include "../include/finite_difference_solver.h"
#include "../include/Riemann_solver.h"

/* This function use Godunov scheme to solve 1-D Euler
 * equations of motion on Lagrange coordinate.
 *
 *[config] is the array of configuration data, the detail
 *         could be seen in the comments of the main function.
 *[m]      is the number of the grids.
 */

void Godunov_solver_source
(double * config, int m, 
 double * RHO[], double * U[], double * P[], double * E[], double * X[], double * MASS,
 double * RHOL, double * UL, double * PL, double * RHOR, double * UR, double * PR, double * cpu_time)
{
  int j = 0, k = 1;  /* j is a frequently used index for
		      * spatial variables. n is a frequ-
		      * ently used index for the time
		      * step.
		      */

  clock_t tic, toc;
  double sum = 0.0;

  int const N = (int)(config[4]);  // the number of time steps
  double const eps = config[3];    // the largest value could be
                                   // seen as zero
  double const h = config[2];      // the length of the initial spatial grids
  double tau = config[1];          // the length of the time step
  double const gamma = config[0];      // the constant of the perfect gas
  double const t_all = config[5];      // the total time
  double const CFL   = config[6];      // the CFL number
  
  double const zeta = (gamma-1.0)/(gamma+1.0);

  double u_L, p_L, rho_L;
  double u_R, p_R, rho_R;
  double c_L, c_R;               // the speed of sound
  double c_star_L,c_star_R;
  double u_star, p_star;
  double *u_mid = malloc((m + 1) * sizeof(double)); // the Riemann solutions
  double *p_mid = malloc((m + 1) * sizeof(double));
  int CRW[2];
  double h_L, h_R;
  double h_S_max, time_n = 0.0;

//------------THE MAIN LOOP-------------
  for(k = 1; k <= N; ++k)
  {
      h_S_max = 1.0/0.0;
    tic = clock();

    for(j = 0; j <= m; ++j)
    { /*
       *  j-1          j          j+1
       * j-1/2  j-1  j+1/2   j   j+3/2  j+1
       *   o-----X-----o-----X-----o-----X--...
       */

      if(j)
      { rho_L = RHO[k-1][j-1];
	u_L   =   U[k-1][j-1];
	p_L   =   P[k-1][j-1];
	h_L   = X[k-1][j] - X[k-1][j-1];
      }
      else
      { rho_L = RHOL[k-1];
	u_L   = UL[k-1];
	p_L   = PL[k-1];
	h_L   = h;
      }

      if(j < m)
      { rho_R = RHO[k-1][j];
	u_R   =   U[k-1][j];
	p_R   =   P[k-1][j];
	h_R   = X[k-1][j+1] - X[k-1][j];
      }
      else
      { rho_R = RHOR[k-1];
	u_R = UR[k-1];
	p_R = PR[k-1];
	h_R = h;
      } /* initialize the values */

      c_L = sqrt(gamma * p_L / rho_L);
      c_R = sqrt(gamma * p_R / rho_R);

      h_S_max = fmin(h_S_max,h_L/c_L);
      h_S_max = fmin(h_S_max,h_R/c_R);

      if((p_L < eps) || (p_R < eps) || (rho_L < eps) || (rho_R < eps)||ISNAN(p_L)||ISNAN(p_R)||ISNAN(u_L)||ISNAN(u_R)||ISNAN(rho_L)||ISNAN(rho_R))
	  printf("error on %d\t%d(t_n,x)\n", k, j);

//===============Solve Riemann Problem==============

	Riemann_solver_exact(&u_star, &p_star, gamma, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, eps, 500);

	u_mid[j] = u_star;
	p_mid[j] = p_star;
    }

//==================================================
	
    tau = CFL*h_S_max;
    //printf("tau=%g\n",tau);
    if ((time_n + tau) >= t_all)
        tau = t_all - time_n + eps;
    
    for(j = 0; j <= m; ++j)
	  X[k][j] = X[k-1][j] + tau*u_mid[j];//motion along the contact discontinuity

    
//===============THE CORE ITERATION=================(On Lagrange Coordinate)
    for(j = 0; j < m; ++j)
    { /*
       *  j-1          j          j+1
       * j-1/2  j-1  j+1/2   j   j+3/2  j+1
       *   o-----X-----o-----X-----o-----X--...
       */
		RHO[k][j] = 1 / ( 1/RHO[k-1][j] + tau/MASS[j]*(u_mid[j+1] - u_mid[j]));
		U[k][j] = U[k-1][j] - tau/MASS[j]*(p_mid[j+1] - p_mid[j]);
		E[k][j] = E[k-1][j] - tau/MASS[j]*(p_mid[j+1]*u_mid[j+1] - p_mid[j]*u_mid[j]);
		P[k][j] = (E[k][j] - 0.5*U[k][j]*U[k][j]) * (gamma - 1.0) * RHO[k][j];
    } //forward Euler


//==================================================

    toc = clock();
    cpu_time[k] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    sum += cpu_time[k];

    time_n += tau;
    if(time_n > t_all)
	{
	    printf("Time is up in time step %d.\n", k);
	    break;
	}
  }

  free(u_mid);
  free(p_mid);

  printf("The cost of CPU time for 1D-Godunov scheme for this problem is %g seconds.\n", sum);
//------------END OF THE MAIN LOOP-------------

}
