/**
 * @file  bound_cond_slope_limiter_x.c
 * @brief This is a function to set boundary conditions and use the slope limiter in x-direction of two dimension.
 */
#include <stdio.h>
#include <stdbool.h>

#include "../include/var_struc.h"
#include "../include/inter_process.h"


/**
 * @brief This function apply the minmod limiter to the slope in the x-direction of two dimension.
 * @param[in] m:          Number of the x-grids: n_x.
 * @param[in] n:          Number of the y-grids: n_y.
 * @param[in] nt:         Current plot time step for computing updates of conservative variables.
 * @param[in] CV:         Structure of cell variable data.
 * @param[in,out] bfv_L:  Fluid variables at left boundary.
 * @param[in,out] bfv_R:  Fluid variables at right boundary.
 * @param[in,out] bfv_D:  Fluid variables at downside boundary.
 * @param[in,out] bfv_U:  Fluid variables at upper boundary.
 * @param[in] find_bound_x: Whether the boundary conditions in x-direction have been found.
 * @param[in] Slope:      Are there slopes? (true: 2nd-order / false: 1st-order)
 * @param[in] t_c:        Time of current time step.
 * @return find_bound_x:  Whether the boundary conditions in x-direction have been found.
 */
_Bool bound_cond_slope_limiter_x(const int m, const int n, const int nt, struct cell_var_stru * CV, struct b_f_var * bfv_L, struct b_f_var * bfv_R,
				 struct b_f_var * bfv_D, struct b_f_var * bfv_U, _Bool find_bound_x, const _Bool Slope, const double t_c)
{
    int const bound_x = (int)(config[17]);// the boundary condition in x-direction
    int const bound_y = (int)(config[18]);// the boundary condition in y-direction
    double const h_x  = config[10];       // the length of the initial x-spatial grids
    int i, j;
    for(i = 0; i < n; ++i)
	switch (bound_x)
	    {
	    case -1: // initial boudary conditions
		if(find_bound_x)
		    break;
		else if(!i)
		    printf("Initial boudary conditions in x direction at time %g .\n", t_c);
		bfv_L[i].U   =   CV->U[0][i]; bfv_R[i].U   =   CV->U[m-1][i];
		bfv_L[i].V   =   CV->V[0][i]; bfv_R[i].V   =   CV->V[m-1][i];
		bfv_L[i].P   =   CV->P[0][i]; bfv_R[i].P   =   CV->P[m-1][i];
		bfv_L[i].RHO = CV->RHO[0][i]; bfv_R[i].RHO = CV->RHO[m-1][i];
		break;
	    case -2: // reflective boundary conditions
		if(!find_bound_x && !i)
		    printf("Reflective boudary conditions in x direction.\n");
		bfv_L[i].U   = - CV[nt].U[0][i]; bfv_R[i].U   = - CV[nt].U[m-1][i];
		bfv_L[i].V   =   CV[nt].V[0][i]; bfv_R[i].V   =   CV[nt].V[m-1][i];
		bfv_L[i].P   =   CV[nt].P[0][i]; bfv_R[i].P   =   CV[nt].P[m-1][i];
		bfv_L[i].RHO = CV[nt].RHO[0][i]; bfv_R[i].RHO = CV[nt].RHO[m-1][i];
		break;
	    case -4: // free boundary conditions
		if(!find_bound_x && !i)
		    printf("Free boudary conditions in x direction.\n");
		bfv_L[i].U   =   CV[nt].U[0][i]; bfv_R[i].U   =   CV[nt].U[m-1][i];
		bfv_L[i].V   =   CV[nt].V[0][i]; bfv_R[i].V   =   CV[nt].V[m-1][i];
		bfv_L[i].P   =   CV[nt].P[0][i]; bfv_R[i].P   =   CV[nt].P[m-1][i];
		bfv_L[i].RHO = CV[nt].RHO[0][i]; bfv_R[i].RHO = CV[nt].RHO[m-1][i];
		break;
	    case -7: // periodic boundary conditions
		if(!find_bound_x && !i)
		    printf("Periodic boudary conditions in x direction.\n");
		bfv_L[i].U   =   CV[nt].U[m-1][i]; bfv_R[i].U   =   CV[nt].U[0][i];
		bfv_L[i].V   =   CV[nt].V[m-1][i]; bfv_R[i].V   =   CV[nt].V[0][i];
		bfv_L[i].P   =   CV[nt].P[m-1][i]; bfv_R[i].P   =   CV[nt].P[0][i];
		bfv_L[i].RHO = CV[nt].RHO[m-1][i]; bfv_R[i].RHO = CV[nt].RHO[0][i];
		break;
	    case -24: // reflective + free boundary conditions
		if(!find_bound_x && !i)
		    printf("Reflective + Free boudary conditions in x direction.\n");
		bfv_L[i].U   = - CV[nt].U[0][i]; bfv_R[i].U   =   CV[nt].U[m-1][i];
		bfv_L[i].V   =   CV[nt].V[0][i]; bfv_R[i].V   =   CV[nt].V[m-1][i];
		bfv_L[i].P   =   CV[nt].P[0][i]; bfv_R[i].P   =   CV[nt].P[m-1][i];
		bfv_L[i].RHO = CV[nt].RHO[0][i]; bfv_R[i].RHO = CV[nt].RHO[m-1][i];
		break;
	    default:
		printf("No suitable boundary coditions in x direction!\n");
		return false;
	    }
    if (Slope)
	{
#pragma acc parallel loop
#pragma omp parallel for  schedule(dynamic, 8)
	    for(i = 0; i < n; ++i)
		{
		    minmod_limiter_2D_x(false, m, i, find_bound_x, CV->s_u,   CV[nt].U,   bfv_L[i].U,   bfv_R[i].U,   h_x);
		    minmod_limiter_2D_x(false, m, i, find_bound_x, CV->s_v,   CV[nt].V,   bfv_L[i].V,   bfv_R[i].V,   h_x);
		    minmod_limiter_2D_x(false, m, i, find_bound_x, CV->s_p,   CV[nt].P,   bfv_L[i].P,   bfv_R[i].P,   h_x);
		    minmod_limiter_2D_x(false, m, i, find_bound_x, CV->s_rho, CV[nt].RHO, bfv_L[i].RHO, bfv_R[i].RHO, h_x);
		} // End of parallel region

	    for(i = 0; i < n; ++i)
		switch(bound_x)
		    {
		    case -2: // reflective boundary conditions
			bfv_L[i].SU   =   CV->s_u[0][i];   bfv_R[i].SU   =   CV->s_u[m-1][i];
			break;
		    case -7: // periodic boundary conditions
			bfv_L[i].SU   =   CV->s_u[m-1][i]; bfv_R[i].SU   =   CV->s_u[0][i];
			bfv_L[i].SV   =   CV->s_v[m-1][i]; bfv_R[i].SV   =   CV->s_v[0][i];
			bfv_L[i].SP   =   CV->s_p[m-1][i]; bfv_R[i].SP   =   CV->s_p[0][i];
			bfv_L[i].SRHO = CV->s_rho[m-1][i]; bfv_R[i].SRHO = CV->s_rho[0][i];
			break;
		    case -24: // reflective + free boundary conditions
			bfv_L[i].SU   =   CV->s_u[0][i];
			break;
		    }

	    for(j = 0; j < m; ++j)
		switch(bound_y)
		    {
		    case -2: case -4: case -24: // reflective OR free boundary conditions in y-direction
			bfv_D[j].SU   =   CV->s_u[j][0];   bfv_U[j].SU   =   CV->s_u[j][n-1];
			bfv_D[j].SV   =   CV->s_v[j][0];   bfv_U[j].SV   =   CV->s_v[j][n-1];
			bfv_D[j].SP   =   CV->s_p[j][0];   bfv_U[j].SP   =   CV->s_p[j][n-1];
			bfv_D[j].SRHO = CV->s_rho[j][0];   bfv_U[j].SRHO = CV->s_rho[j][n-1];
			break;
		    case -7: // periodic boundary conditions in y-direction
			bfv_D[j].SU   =   CV->s_u[j][n-1]; bfv_U[j].SU   =   CV->s_u[j][0];
			bfv_D[j].SV   =   CV->s_v[j][n-1]; bfv_U[j].SV   =   CV->s_v[j][0];
			bfv_D[j].SP   =   CV->s_p[j][n-1]; bfv_U[j].SP   =   CV->s_p[j][0];
			bfv_D[j].SRHO = CV->s_rho[j][n-1]; bfv_U[j].SRHO = CV->s_rho[j][0];
			break;
		    }
	}
    return true;
}
