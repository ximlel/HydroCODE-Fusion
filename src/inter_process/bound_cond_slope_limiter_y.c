/**
 * @file  bound_cond_slope_limiter_y.c
 * @brief This is a function to set boundary conditions and use the slope limiter in y-direction of two dimension.
 */
#include <stdio.h>
#include <stdbool.h>

#include "../include/var_struc.h"
#include "../include/inter_process.h"


/**
 * @brief This function apply the minmod limiter to the slope in the y-direction of two dimension.
 * @param[in] m:          Number of the x-grids: n_x.
 * @param[in] n:          Number of the y-grids: n_y.
 * @param[in] nt:         Current plot time step for computing updates of conservative variables.
 * @param[in] CV:         Structure of cell variable data.
 * @param[in,out] bfv_L:  Fluid variables at left boundary.
 * @param[in,out] bfv_R:  Fluid variables at right boundary.
 * @param[in,out] bfv_D:  Fluid variables at downside boundary.
 * @param[in,out] bfv_U:  Fluid variables at upper boundary.
 * @param[in] find_bound_y: Whether the boundary conditions in y-direction have been found.
 * @param[in] Slope:      Are there slopes? (true: 2nd-order / false: 1st-order)
 * @param[in] t_c:        Time of current time step.
 * @return find_bound_y:  Whether the boundary conditions in y-direction have been found.
 */
_Bool bound_cond_slope_limiter_y(const int m, const int n, const int nt, struct cell_var_stru * CV, struct b_f_var * bfv_L, struct b_f_var * bfv_R,
				 struct b_f_var * bfv_D, struct b_f_var * bfv_U, _Bool find_bound_y, const _Bool Slope, const double t_c)
{
    int const bound_x = (int)(config[17]);// the boundary condition in x-direction
    int const bound_y = (int)(config[18]);// the boundary condition in y-direction
    double const h_y  = config[11];       // the length of the initial y-spatial grids
    int i, j;
    for(j = 0; j < m; ++j)
	switch (bound_y)
	    {
	    case -1: // initial boudary conditions
		if(find_bound_y)
		    break;
		else if (!j)
		    printf("Initial boudary conditions in y direction at time %g .\n", t_c);
		bfv_D[j].U   =   CV->U[j][0]; bfv_U[j].U   =   CV->U[j][n-1];
		bfv_D[j].V   =   CV->V[j][0]; bfv_U[j].V   =   CV->V[j][n-1];
		bfv_D[j].P   =   CV->P[j][0]; bfv_U[j].P   =   CV->P[j][n-1];
		bfv_D[j].RHO = CV->RHO[j][0]; bfv_U[j].RHO = CV->RHO[j][n-1];
		break;
	    case -2: // reflective boundary conditions
		if(!find_bound_y && !j)
		    printf("Reflective boudary conditions in y direction.\n");
		bfv_D[j].U   =   CV[nt].U[j][0]; bfv_U[j].U   =   CV[nt].U[j][n-1];
		bfv_D[j].V   = - CV[nt].V[j][0]; bfv_U[j].V   = - CV[nt].V[j][n-1];
		bfv_D[j].P   =   CV[nt].P[j][0]; bfv_U[j].P   =   CV[nt].P[j][n-1];
		bfv_D[j].RHO = CV[nt].RHO[j][0]; bfv_U[j].RHO = CV[nt].RHO[j][n-1];
		break;
	    case -4: // free boundary conditions
		if(!find_bound_y && !j)
		    printf("Free boudary conditions in y direction.\n");
		bfv_D[j].U   =   CV[nt].U[j][0]; bfv_U[j].U   =   CV[nt].U[j][n-1];
		bfv_D[j].V   =   CV[nt].V[j][0]; bfv_U[j].V   =   CV[nt].V[j][n-1];
		bfv_D[j].P   =   CV[nt].P[j][0]; bfv_U[j].P   =   CV[nt].P[j][n-1];
		bfv_D[j].RHO = CV[nt].RHO[j][0]; bfv_U[j].RHO = CV[nt].RHO[j][n-1];
		break;
	    case -7: // periodic boundary conditions
		if(!find_bound_y && !j)
		    printf("Periodic boudary conditions in y direction.\n");
		bfv_D[j].U   =   CV[nt].U[j][n-1]; bfv_U[j].U   =   CV[nt].U[j][0];
		bfv_D[j].V   =   CV[nt].V[j][n-1]; bfv_U[j].V   =   CV[nt].V[j][0];
		bfv_D[j].P   =   CV[nt].P[j][n-1]; bfv_U[j].P   =   CV[nt].P[j][0];
		bfv_D[j].RHO = CV[nt].RHO[j][n-1]; bfv_U[j].RHO = CV[nt].RHO[j][0];
		break;
	    case -24: // reflective + free boundary conditions
		if(!find_bound_y && !j)
		    printf("Reflective + Free boudary conditions in y direction.\n");
		bfv_D[j].U   =   CV[nt].U[j][0]; bfv_U[j].U   =   CV[nt].U[j][n-1];
		bfv_D[j].V   = - CV[nt].V[j][0]; bfv_U[j].V   =   CV[nt].V[j][n-1];
		bfv_D[j].P   =   CV[nt].P[j][0]; bfv_U[j].P   =   CV[nt].P[j][n-1];
		bfv_D[j].RHO = CV[nt].RHO[j][0]; bfv_U[j].RHO = CV[nt].RHO[j][n-1];
		break;
	    default:
		printf("No suitable boundary coditions in y direction!\n");
		return false;
	    }
    if (Slope)
	{
#pragma omp parallel for  schedule(dynamic, 8)
	    for(j = 0; j < m; ++j)
		{
		    minmod_limiter(false, n, find_bound_y, CV->t_u[j],   CV[nt].U[j],   bfv_D[j].U,   bfv_U[j].U,   h_y);
		    minmod_limiter(false, n, find_bound_y, CV->t_v[j],   CV[nt].V[j],   bfv_D[j].V,   bfv_U[j].V,   h_y);
		    minmod_limiter(false, n, find_bound_y, CV->t_p[j],   CV[nt].P[j],   bfv_D[j].P,   bfv_U[j].P,   h_y);
		    minmod_limiter(false, n, find_bound_y, CV->t_rho[j], CV[nt].RHO[j], bfv_D[j].RHO, bfv_U[j].RHO, h_y);
		} // End of parallel region

	    for(j = 0; j < m; ++j)
		switch(bound_y)
		    {
		    case -2: // reflective boundary conditions
			bfv_D[j].TV   =   CV->t_v[j][0];   bfv_U[j].TV   =   CV->t_v[j][n-1];
			break;
		    case -7: // periodic boundary conditions
			bfv_D[j].TU   =   CV->t_u[j][n-1]; bfv_U[j].TU   =   CV->t_u[j][0];
			bfv_D[j].TV   =   CV->t_v[j][n-1]; bfv_U[j].TV   =   CV->t_v[j][0];
			bfv_D[j].TP   =   CV->t_p[j][n-1]; bfv_U[j].TP   =   CV->t_p[j][0];
			bfv_D[j].TRHO = CV->t_rho[j][n-1]; bfv_U[j].TRHO = CV->t_rho[j][0];
			break;
		    case -24: // reflective + free boundary conditions
			bfv_D[j].TV   =   CV->t_v[j][0];
			break;
		    }

	    for(i = 0; i < n; ++i)
		switch(bound_x)
		    {
		    case -2: case -4: case -24: // reflective OR free boundary conditions in x-direction
			bfv_L[i].TU   =   CV->t_u[0][i];   bfv_R[i].TU   =   CV->t_u[m-1][i];
			bfv_L[i].TV   =   CV->t_v[0][i];   bfv_R[i].TV   =   CV->t_v[m-1][i];
			bfv_L[i].TP   =   CV->t_p[0][i];   bfv_R[i].TP   =   CV->t_p[m-1][i];
			bfv_L[i].TRHO = CV->t_rho[0][i];   bfv_R[i].TRHO = CV->t_rho[m-1][i];
			break;
		    case -7: // periodic boundary conditions in x-direction
			bfv_L[i].TU   =   CV->t_u[m-1][i]; bfv_R[i].TU   =   CV->t_u[0][i];
			bfv_L[i].TV   =   CV->t_v[m-1][i]; bfv_R[i].TV   =   CV->t_v[0][i];
			bfv_L[i].TP   =   CV->t_p[m-1][i]; bfv_R[i].TP   =   CV->t_p[0][i];
			bfv_L[i].TRHO = CV->t_rho[m-1][i]; bfv_R[i].TRHO = CV->t_rho[0][i];
			break;
		    }
	}
    return true;
}
