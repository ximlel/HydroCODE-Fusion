/**
 * @file  bound_cond_slope_limiter.c
 * @brief This is a function to set boundary conditions and use the slope limiter in one dimension.
 */
#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>

#include "../include/var_struc.h"
#include "../include/inter_process.h"


/**
 * @brief This function apply the minmod limiter to the slope in one dimension.
 * @param[in] NO_h:       Whether there are moving grid point coordinates.
 *                  - true: There are moving spatial grid point coordinates *X.
 *                  - false: There is fixed spatial grid length.
 * @param[in] m:          Number of the grids.
 * @param[in] nt:         Current plot time step for computing updates of conservative variables.
 * @param[in] CV:         Structure of cell variable data.
 * @param[in,out] bfv_L:  Fluid variables at left boundary.
 * @param[in,out] bfv_R:  Fluid variables at right boundary.
 * @param[in] find_bound: Whether the boundary conditions have been found.
 * @param[in] Slope:      Are there slopes? (true: 2nd-order / false: 1st-order)
 * @param[in] t_c:        Time of current time step.
 * @param[in] ...:        Variable parameter if NO_h is true.
 *            - \b double \c *X: Array of moving spatial grid point coordinates.
 * @return find_bound:    Whether the boundary conditions have been found.
 */
_Bool bound_cond_slope_limiter(const _Bool NO_h, const int m, const int nt, struct cell_var_stru CV,
			       struct b_f_var * bfv_L, struct b_f_var * bfv_R, _Bool find_bound, const _Bool Slope, const double t_c, ...)
{
    va_list ap;
    va_start(ap, t_c);
    int const bound = (int)(config[17]);// the boundary condition in x-direction
    double const h  = config[10];       // the length of the initial x-spatial grids
    double * X;
    if (NO_h)
	X  = va_arg(ap, double *);

    switch (bound)
	{
	case -1: // initial boudary conditions
	    if(find_bound)
		break;
	    else
		printf("Initial boudary conditions in x direction at time %g .\n", t_c);
	    bfv_L->U   =   CV.U[0][0]; bfv_R->U   =   CV.U[0][m-1];
	    bfv_L->P   =   CV.P[0][0]; bfv_R->P   =   CV.P[0][m-1];
	    bfv_L->RHO = CV.RHO[0][0]; bfv_R->RHO = CV.RHO[0][m-1];
	    break;
	case -2: // reflective boundary conditions
	    if(!find_bound)
		printf("Reflective boudary conditions in x direction.\n");
	    bfv_L->U   = - CV.U[nt][0]; bfv_R->U   = - CV.U[nt][m-1];
	    bfv_L->P   =   CV.P[nt][0]; bfv_R->P   =   CV.P[nt][m-1];
	    bfv_L->RHO = CV.RHO[nt][0]; bfv_R->RHO = CV.RHO[nt][m-1];
	    break;
	case -4: // free boundary conditions
	    if(!find_bound)
		printf("Free boudary conditions in x direction.\n");
	    bfv_L->U   =   CV.U[nt][0]; bfv_R->U   =   CV.U[nt][m-1];
	    bfv_L->P   =   CV.P[nt][0]; bfv_R->P   =   CV.P[nt][m-1];
	    bfv_L->RHO = CV.RHO[nt][0]; bfv_R->RHO = CV.RHO[nt][m-1];
	    break;
	case -5: // periodic boundary conditions
	    if(!find_bound)
		printf("Periodic boudary conditions in x direction.\n");
	    bfv_L->U   =   CV.U[nt][m-1]; bfv_R->U   =   CV.U[nt][0];
	    bfv_L->P   =   CV.P[nt][m-1]; bfv_R->P   =   CV.P[nt][0];
	    bfv_L->RHO = CV.RHO[nt][m-1]; bfv_R->RHO = CV.RHO[nt][0];
	    break;
	case -24: // reflective + free boundary conditions
	    if(!find_bound)
		printf("Reflective + Free boudary conditions in x direction.\n");
	    bfv_L->U   = - CV.U[nt][0]; bfv_R->U   =   CV.U[nt][m-1];
	    bfv_L->P   =   CV.P[nt][0]; bfv_R->P   =   CV.P[nt][m-1];
	    bfv_L->RHO = CV.RHO[nt][0]; bfv_R->RHO = CV.RHO[nt][m-1];
	    break;
	default:
	    printf("No suitable boundary coditions in x direction!\n");
	    return false;
	}

    if (NO_h)
	{
	    switch (bound)
		{
		case -1: // initial boudary conditions
		    bfv_L->H  = h; bfv_R->H = h;
		    break;
		case -5: // periodic boundary conditions
		    bfv_L->H = X[m] - X[m-1];
		    bfv_R->H = X[1] - X[0];
		    break;
		case -2: case -4: case -24:
		    bfv_L->H = X[1] - X[0];
		    bfv_R->H = X[m] - X[m-1];
		    break;
		}
	}
//=================Initialize slopes=====================
      // Reconstruct slopes
    if (Slope)
	{
    if (NO_h)
	{
	    minmod_limiter(NO_h, m, find_bound, CV.d_u,   CV.U[nt],   bfv_L->U,   bfv_R->U,   bfv_L->H, bfv_R->H, X);
	    minmod_limiter(NO_h, m, find_bound, CV.d_p,   CV.P[nt],   bfv_L->P,   bfv_R->P,   bfv_L->H, bfv_R->H, X);
	    minmod_limiter(NO_h, m, find_bound, CV.d_rho, CV.RHO[nt], bfv_L->RHO, bfv_R->RHO, bfv_L->H, bfv_R->H, X);
	}
	else
	{
	    minmod_limiter(NO_h, m, find_bound, CV.d_u,   CV.U[nt],   bfv_L->U,   bfv_R->U,   h);
	    minmod_limiter(NO_h, m, find_bound, CV.d_p,   CV.P[nt],   bfv_L->P,   bfv_R->P,   h);
	    minmod_limiter(NO_h, m, find_bound, CV.d_rho, CV.RHO[nt], bfv_L->RHO, bfv_R->RHO, h);
	}

	    switch(bound)
		{
		case -2: // reflective boundary conditions
		    bfv_L->SU = CV.d_u[0]; bfv_R->SU = CV.d_u[m-1];
		    break;
		case -5: // periodic boundary conditions
		    bfv_L->SU   =   CV.d_u[m-1]; bfv_R->SU   =   CV.d_u[0];
		    bfv_L->SP   =   CV.d_p[m-1]; bfv_R->SP   =   CV.d_p[0];
		    bfv_L->SRHO = CV.d_rho[m-1]; bfv_R->SRHO = CV.d_rho[0];
		    break;
		case -24: // reflective + free boundary conditions
		    bfv_L->SU = CV.d_u[0];
		    break;
		}
	}
    va_end(ap);
    return true;
}
