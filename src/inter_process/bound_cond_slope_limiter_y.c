#include <stdio.h>
#include <stdbool.h>

#include "../include/var_struc.h"
#include "../include/inter_process.h"


_Bool bound_cond_slope_limiter_y(const int m, const int n, const int nt, struct cell_var_stru * CV,
				 struct b_f_var * bfv_D, struct b_f_var * bfv_U, _Bool find_bound_y, const _Bool Slope, const double t_c)
{
    int const bound_y = (int)(config[18]);// the boundary condition in y-direction
    double const h_y  = config[11];       // the length of the initial y-spatial grids
    int j;
    for(j = 0; j < m; ++j)
	switch (bound_y)
	    {
	    case -1: // initial boudary conditions
		if(find_bound_y)
		    break;
		else if (j == 1)
		    printf("Initial boudary conditions in y direction at time %g .\n", t_c);		  
		bfv_D[j].U   =   CV->U[j][0]; bfv_U[j].U   =   CV->U[j][n-1];
		bfv_D[j].V   =   CV->V[j][0]; bfv_U[j].V   =   CV->V[j][n-1];
		bfv_D[j].P   =   CV->P[j][0]; bfv_U[j].P   =   CV->P[j][n-1];
		bfv_D[j].RHO = CV->RHO[j][0]; bfv_U[j].RHO = CV->RHO[j][n-1];
		break;
	    case -2: // reflective boundary conditions
		if(!find_bound_y && j == 1)
		    printf("Reflective boudary conditions in y direction.\n");
		bfv_D[j].U   =   CV[nt].U[j][0]; bfv_U[j].U   =   CV[nt].U[j][n-1];
		bfv_D[j].V   = - CV[nt].V[j][0]; bfv_U[j].V   = - CV[nt].V[j][n-1];
		bfv_D[j].P   =   CV[nt].P[j][0]; bfv_U[j].P   =   CV[nt].P[j][n-1];
		bfv_D[j].RHO = CV[nt].RHO[j][0]; bfv_U[j].RHO = CV[nt].RHO[j][n-1];
		break;
	    case -4: // free boundary conditions
		if(!find_bound_y && j == 1)
		    printf("Free boudary conditions in y direction.\n");
		bfv_D[j].U   =   CV[nt].U[j][0]; bfv_U[j].U   =   CV[nt].U[j][n-1];
		bfv_D[j].V   =   CV[nt].V[j][0]; bfv_U[j].V   =   CV[nt].V[j][n-1];
		bfv_D[j].P   =   CV[nt].P[j][0]; bfv_U[j].P   =   CV[nt].P[j][n-1];
		bfv_D[j].RHO = CV[nt].RHO[j][0]; bfv_U[j].RHO = CV[nt].RHO[j][n-1];
		break;
	    case -5: // periodic boundary conditions
		if(!find_bound_y && j == 1)
		    printf("Periodic boudary conditions in y direction.\n");
		bfv_D[j].U   =   CV[nt].U[j][n-1]; bfv_U[j].U   =   CV[nt].U[j][0];
		bfv_D[j].V   =   CV[nt].V[j][n-1]; bfv_U[j].V   =   CV[nt].V[j][0];
		bfv_D[j].P   =   CV[nt].P[j][n-1]; bfv_U[j].P   =   CV[nt].P[j][0];
		bfv_D[j].RHO = CV[nt].RHO[j][n-1]; bfv_U[j].RHO = CV[nt].RHO[j][0];
		break;
	    case -24: // reflective + free boundary conditions
		if(!find_bound_y && j == 1)
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
	    for(j = 0; j < m; ++j)
		{
		    minmod_limiter(false, n, find_bound_y, CV->t_u[j],   CV[nt].U[j],   bfv_D[j].U,   bfv_U[j].U,   h_y);
		    minmod_limiter(false, n, find_bound_y, CV->t_v[j],   CV[nt].V[j],   bfv_D[j].V,   bfv_U[j].V,   h_y);
		    minmod_limiter(false, n, find_bound_y, CV->t_p[j],   CV[nt].P[j],   bfv_D[j].P,   bfv_U[j].P,   h_y);
		    minmod_limiter(false, n, find_bound_y, CV->t_rho[j], CV[nt].RHO[j], bfv_D[j].RHO, bfv_U[j].RHO, h_y);
		}

	    for(j = 0; j < m; ++j)
		switch(bound_y)
		    {
		    case -2: // reflective boundary conditions
			bfv_D[j].TV   =   CV->t_v[j][0]; bfv_U[j].TV   =   CV->t_v[j][n-1];
		    case -5: // periodic boundary conditions
			bfv_D[j].SU   =   CV->s_u[j][n-1]; bfv_U[j].SU   =   CV->s_u[j][0];
			bfv_D[j].SV   =   CV->s_v[j][n-1]; bfv_U[j].SV   =   CV->s_v[j][0];
			bfv_D[j].SP   =   CV->s_p[j][n-1]; bfv_U[j].SP   =   CV->s_p[j][0];
			bfv_D[j].SRHO = CV->s_rho[j][n-1]; bfv_U[j].SRHO = CV->s_rho[j][0];
			bfv_D[j].TU   =   CV->t_u[j][n-1]; bfv_U[j].TU   =   CV->t_u[j][0];
			bfv_D[j].TV   =   CV->t_v[j][n-1]; bfv_U[j].TV   =   CV->t_v[j][0];
			bfv_D[j].TP   =   CV->t_p[j][n-1]; bfv_U[j].TP   =   CV->t_p[j][0];
			bfv_D[j].TRHO = CV->t_rho[j][n-1]; bfv_U[j].TRHO = CV->t_rho[j][0];
			break;
		    case -24: // reflective + free boundary conditions
			bfv_D[j].TV   =   CV->t_v[j][0];
			break;
		    }
		    
	    for(j = 0; j < m; ++j)
		switch(bound_y)
		    {
		    case -2: case -4: case -24:
			bfv_D[j].SU   =   CV->s_u[j][0]; bfv_U[j].SU   =   CV->s_u[j][n-1];
			bfv_D[j].SV   =   CV->s_v[j][0]; bfv_U[j].SV   =   CV->s_v[j][n-1];
			bfv_D[j].SP   =   CV->s_p[j][0]; bfv_U[j].SP   =   CV->s_p[j][n-1];
			bfv_D[j].SRHO = CV->s_rho[j][0]; bfv_U[j].SRHO = CV->s_rho[j][n-1];
			break;
		    }
	}
    return true;
}
