#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN32
#define ISNAN(a) _isnan((a))
#elif __linux__
#define ISNAN(a) isnan((a))
#endif

#include "../include/Riemann_solver.h"


/* This function use GRP scheme to solve 1-D Euler
 * equations of motion on Lagrange coordinate.
 *
 *[config] is the array of configuration data, the detail
 *         could be seen in the comments of the main function.
 *[m]      is the number of the grids.
 */

void GRP_solver_source
(double * config, int m, double * RHO[], double * U[], double * P[],
 double * E[], double * X[], double * MASS,
 double * RHOL, double * UL, double * PL,
 double * RHOR, double * UR, double * PR,
 double * SRHOL, double * SUL, double * SPL,
 double * SRHOR, double * SUR, double * SPR,
 double * cpu_time)
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
  double tau = config[1];    // the length of the time step
  double const gamma = config[0];      // the constant of the perfect gas
  double const t_all = config[5];      // the total time
  double const CFL   = config[6];      // the CFL number

  double s_rho[m], s_u[m], s_p[m];
  double s_L, s_R;
  double u_L, p_L, rho_L;
  double u_R, p_R, rho_R;
  double c_L, c_R;
  double t_u_L, t_p_L, t_rho_L;
  double t_u_R, t_p_R, t_rho_R;
  double dire[4], mid[4];  //RHO_L_t,U_t,P_t,RHO_R_t.
  double s_u_L, s_p_L, s_rho_L;
  double s_u_R, s_p_R, s_rho_R;
  double slope_temp;
  double X_mass[m];


  double U_next[m+1], P_next[m+1], RHO_next_L[m+1], RHO_next_R[m+1];
  
  double U_F[m+1], P_F[m+1];// the numerical flux at t_{n+1/2}

  double s_time = 0.0;//1.9; // the paramater in slope limiters

  double h_L, h_R;
  double h_S_max, time_n = 0.0;

//=============initialize the values===========

  for(j = 0; j < m; ++j)
  {
    if(j)
    {
      s_u_L = (U[0][j] - U[0][j-1]) / h;
      s_p_L = (P[0][j] - P[0][j-1]) / h;
      s_rho_L = (RHO[0][j] - RHO[0][j-1]) / h;
    }
    else
    {
      s_u_L = (U[0][j] - UL[0]) / h;
      s_p_L = (P[0][j] - PL[0]) / h;
      s_rho_L = (RHO[0][j] - RHOL[0]) / h;
    }
    if(j < m-1)
    {
      s_u_R = (U[0][j+1] - U[0][j]) / h;
      s_p_R = (P[0][j+1] - P[0][j]) / h;
      s_rho_R = (RHO[0][j+1] - RHO[0][j]) / h;
    }
    else
    {
      s_u_R = (UR[0] - U[0][j]) / h;
      s_p_R = (PR[0] - P[0][j]) / h;
      s_rho_R = (RHOR[0] - RHO[0][j]) / h;
    }

    if(s_u_L * s_u_R < 0.0)
      s_u[j] = 0.0;
    else if(fabs(s_u_R) < fabs(s_u_L))
      s_u[j] = s_u_R;
    else
      s_u[j] = s_u_L;

    if(s_p_L * s_p_R < 0.0)
      s_p[j] = 0.0;
    else if(fabs(s_p_R) < fabs(s_p_L))
      s_p[j] = s_p_R;
    else
      s_p[j] = s_p_L;

    if(s_rho_L * s_rho_R < 0.0)
      s_rho[j] = 0.0;
    else if(fabs(s_rho_R) < fabs(s_rho_L))
      s_rho[j] = s_rho_R;
    else
      s_rho[j] = s_rho_L;
  }

//===============================================

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
      { rho_L = RHO[k-1][j-1] + 0.5*(X[k-1][j]-X[k-1][j-1])*s_rho[j-1];
	u_L   =   U[k-1][j-1] + 0.5*(X[k-1][j]-X[k-1][j-1])*s_u[j-1];
	p_L   =   P[k-1][j-1] + 0.5*(X[k-1][j]-X[k-1][j-1])*s_p[j-1];
	h_L   = X[k-1][j] - X[k-1][j-1];
      }
      else
      { rho_L = RHOL[k-1];
	u_L = UL[k-1];
	p_L = PL[k-1];
	h_L   = h;
      }

      if(j < m)
      { rho_R = RHO[k-1][j] - 0.5*(X[k-1][j+1]-X[k-1][j])*s_rho[j];
	u_R   =   U[k-1][j] - 0.5*(X[k-1][j+1]-X[k-1][j])*s_u[j];
	p_R   =   P[k-1][j] - 0.5*(X[k-1][j+1]-X[k-1][j])*s_p[j];
 	h_R   = X[k-1][j+1] - X[k-1][j];
     }
      else
      { rho_R = RHOR[k-1];
	u_R = UR[k-1];
	p_R = PR[k-1];
 	h_R = h;
     }

      c_L = sqrt(gamma * p_L / rho_L);
      c_R = sqrt(gamma * p_R / rho_R);

      h_S_max = fmin(h_S_max,h_L/c_L);
      h_S_max = fmin(h_S_max,h_R/c_R);

      if(j)
      { t_u_L   = s_u[j-1]/rho_L;
	t_p_L   = s_p[j-1]/rho_L;
	t_rho_L = s_rho[j-1]/rho_L;
      }
      else
      { t_rho_L = SRHOL[k-1]/rho_L;
	t_u_L = SUL[k-1]/rho_L;
	t_p_L = SPL[k-1]/rho_L;
      }

      if(j < m)
      { t_u_R   = s_u[j]/rho_R;
	t_p_R   = s_p[j]/rho_R;
	t_rho_R = s_rho[j]/rho_R;
      }
      else
      { t_rho_R = SRHOR[k-1]/rho_R;
	t_u_R = SUR[k-1]/rho_R;
	t_p_R = SPR[k-1]/rho_R;
      }//calculate the material derivative

      if((p_L < eps) || (p_R < eps) || (rho_L < eps) || (rho_R < eps)||ISNAN(p_L)||ISNAN(p_R)||ISNAN(u_L)||ISNAN(u_R)||ISNAN(rho_L)||ISNAN(rho_R))
	  printf("error on (%d,%d) (t_n,x)\n", k, j);

//===============GRP scheme======================

      linear_GRP_solver_LAG(dire, mid, rho_L, rho_R, t_rho_L, t_rho_R, u_L, u_R, t_u_L, t_u_R, p_L, p_R, t_p_L, t_p_R, gamma, 1e-9);
    }
    
    tau = CFL*h_S_max;
    //printf("tau=%g\n",tau);
    if ((time_n + tau) >= t_all)
        tau = t_all - time_n + eps;
    
    for(j = 0; j <= m; ++j)
	{
	    U_F[j] = mid[1] + 0.5*tau*dire[1];
	    P_F[j] = mid[2] + 0.5*tau*dire[2];         

	    X[k][j] = X[k-1][j] + tau*U_F[j];

	    RHO_next_L[j] = mid[0] + tau*dire[0];
	    RHO_next_R[j] = mid[3] + tau*dire[3];
	    U_next[j] = mid[1] + tau*dire[1];
	    P_next[j] = mid[2] + tau*dire[2];

//===============================================
 
	    if(j)
		X_mass[j-1] = 0.5*(X[k][j-1]+X[k][j]);
	}

//===============THE CORE ITERATION=================
    for(j = 0; j < m; ++j)
    { /*
       *  j-1          j          j+1
       * j-1/2  j-1  j+1/2   j   j+3/2  j+1
       *   o-----X-----o-----X-----o-----X--...
       */
        RHO[k][j] = 1 / ( 1/RHO[k-1][j] + tau/MASS[j]*(U_F[j+1] - U_F[j]) );
	U[k][j] = U[k-1][j] - tau/MASS[j]*(P_F[j+1] - P_F[j]);
	E[k][j] = E[k-1][j] - tau/MASS[j]*(P_F[j+1]*U_F[j+1] - P_F[j]*U_F[j]);
	P[k][j] = (E[k][j] - 0.5*U[k][j]*U[k][j]) * (gamma - 1.0) * RHO[k][j];
	/* forward Euler */

//-----------------------determind the slope-----------------------------

        s_u[j] = (  U_next[j+1] -   U_next[j])/(X[k][j+1]-X[k][j]);
        s_p[j] = (  P_next[j+1] -   P_next[j])/(X[k][j+1]-X[k][j]);
        s_rho[j] = (RHO_next_L[j+1] - RHO_next_R[j])/(X[k][j+1]-X[k][j]);
   }

    for(j = 0; j < m; ++j)
    {
      if(j)
      {
	s_u_L = s_time*(U[k][j] - U[k][j-1]) / (X_mass[j]-X_mass[j-1]);
	s_p_L = s_time*(P[k][j] - P[k][j-1]) / (X_mass[j]-X_mass[j-1]);
	s_rho_L = s_time*(RHO[k][j] - RHO[k][j-1]) / (X_mass[j]-X_mass[j-1]);
      }
      else
      {
	s_u_L = s_time*(U[k][j] - UL[k]) / (X_mass[j]-X[k][j])/2;
	s_p_L = s_time*(P[k][j] - PL[k]) / (X_mass[j]-X[k][j])/2;
	s_rho_L = s_time*(RHO[k][j] - RHOL[k]) / (X_mass[j]-X[k][j])/2;
      }
      if(j < m-1)
      {
	s_u_R = s_time*(U[k][j+1] - U[k][j]) / (X_mass[j+1]-X_mass[j]);
	s_p_R = s_time*(P[k][j+1] - P[k][j]) / (X_mass[j+1]-X_mass[j]);
	s_rho_R = s_time*(RHO[k][j+1] - RHO[k][j]) / (X_mass[j+1]-X_mass[j]);
      }
      else
      {
	s_u_R = s_time*(UR[k] - U[k][j]) / (X[k][j+1]-X_mass[j])/2;
	s_p_R = s_time*(PR[k] - P[k][j]) / (X[k][j+1]-X_mass[j])/2;
	s_rho_R = s_time*(RHOR[k] - RHO[k][j]) / (X[k][j+1]-X_mass[j])/2;
      }

      if(s_u_L * s_u_R < 0.0)
	s_u[j] = 0.0;
      else if((s_u_L * s_u[j]) < 0.0)
        s_u[j] = 0.0;
      else
      {
	slope_temp = ( (fabs(s_u_L) < fabs(s_u_R)) ? s_u_L : s_u_R );
	if( fabs(slope_temp) < fabs(s_u[j]) )
	    s_u[j] = slope_temp;
      }

      if(s_p_L * s_p_R < 0.0)
	s_p[j] = 0.0;
      else if((s_p_L * s_p[j]) < 0.0)
        s_p[j] = 0.0;
      else
      {
	slope_temp = ((fabs(s_p_L) < fabs(s_p_R)) ? s_p_L : s_p_R);
	if( fabs(slope_temp) < fabs(s_p[j]) )
	    s_p[j] = slope_temp;
      }

      if(s_rho_L * s_rho_R < 0.0)
	s_rho[j] = 0.0;
      else if((s_rho_L * s_rho[j]) < 0.0)
        s_rho[j] = 0.0;
      else
      {
	slope_temp = ((fabs(s_rho_L) < fabs(s_rho_R)) ? s_rho_L : s_rho_R);
	if( fabs(slope_temp) < fabs(s_rho[j]) )
	    s_rho[j] = slope_temp;
      }
    }


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

  printf("The cost of CPU time for 1D-GRP scheme for this problem is %g seconds.\n", sum);
//------------END OF THE MAIN LOOP-------------

}
