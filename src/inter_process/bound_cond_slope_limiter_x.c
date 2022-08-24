#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "../include/var_struc.h"
#include "../include/inter_process.h"


_Bool bound_cond_slope_limiter_x(const int m, const int n, const int nt, struct cell_var_stru * CV,
				 struct b_f_var * bfv_L, struct b_f_var * bfv_R, _Bool find_bound_x)
{
    int const bound_x = (int)(config[17]);// the boundary condition in x-direction
    double const h_x  = config[10];       // the length of the initial x-spatial grids
    int i;
    for(i = 1; i < n; ++i)
	switch (bound_x)
	    {
	    case -1: // initial boudary conditions
		if(find_bound_x)
		    break;
		else
		    printf("Initial boudary conditions in x direction.\n");		  
		bfv_L[i].U   =   CV->U[0][i]; bfv_R[i].U   =   CV->U[m-1][i];
		bfv_L[i].V   =   CV->V[0][i]; bfv_R[i].V   =   CV->V[m-1][i];
		bfv_L[i].P   =   CV->P[0][i]; bfv_R[i].P   =   CV->P[m-1][i];
		bfv_L[i].RHO = CV->RHO[0][i]; bfv_R[i].RHO = CV->RHO[m-1][i];
		break;
	    case -2: // reflective boundary conditions
		if(!find_bound_x)
		    printf("Reflective boudary conditions in x direction.\n");
		bfv_L[i].U   = - CV[nt].U[0][i]; bfv_R[i].U   = - CV[nt].U[m-1][i];
		bfv_L[i].V   =   CV[nt].V[0][i]; bfv_R[i].V   =   CV[nt].V[m-1][i];
		bfv_L[i].P   =   CV[nt].P[0][i]; bfv_R[i].P   =   CV[nt].P[m-1][i];
		bfv_L[i].RHO = CV[nt].RHO[0][i]; bfv_R[i].RHO = CV[nt].RHO[m-1][i];
		break;
	    case -4: // free boundary conditions
		if(!find_bound_x)
		    printf("Free boudary conditions in x direction.\n");
		bfv_L[i].U   =   CV[nt].U[0][i]; bfv_R[i].U   =   CV[nt].U[m-1][i];
		bfv_L[i].V   =   CV[nt].V[0][i]; bfv_R[i].V   =   CV[nt].V[m-1][i];
		bfv_L[i].P   =   CV[nt].P[0][i]; bfv_R[i].P   =   CV[nt].P[m-1][i];
		bfv_L[i].RHO = CV[nt].RHO[0][i]; bfv_R[i].RHO = CV[nt].RHO[m-1][i];
		break;
	    case -5: // periodic boundary conditions
		if(!find_bound_x)
		    printf("Periodic boudary conditions in x direction.\n");
		bfv_L[i].U   =   CV[nt].U[m-1][i]; bfv_R[i].U   =   CV[nt].U[0][i];
		bfv_L[i].V   =   CV[nt].V[m-1][i]; bfv_R[i].V   =   CV[nt].V[0][i];
		bfv_L[i].P   =   CV[nt].P[m-1][i]; bfv_R[i].P   =   CV[nt].P[0][i];
		bfv_L[i].RHO = CV[nt].RHO[m-1][i]; bfv_R[i].RHO = CV[nt].RHO[0][i];
		break;
	    case -24: // reflective + free boundary conditions
		if(!find_bound_x)
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

    for(i = 0; i < n; ++i)
	{
	    minmod_limiter_2D_x(false, m, find_bound_x, i, CV->s_u,   CV[nt].U,   bfv_L[i].U,   bfv_R[i].U,   h_x);
	    minmod_limiter_2D_x(false, m, find_bound_x, i, CV->s_v,   CV[nt].V,   bfv_L[i].V,   bfv_R[i].V,   h_x);
	    minmod_limiter_2D_x(false, m, find_bound_x, i, CV->s_p,   CV[nt].P,   bfv_L[i].P,   bfv_R[i].P,   h_x);
	    minmod_limiter_2D_x(false, m, find_bound_x, i, CV->s_rho, CV[nt].RHO, bfv_L[i].RHO, bfv_R[i].RHO, h_x);
	}

    for(i = 1; i < n; ++i)
	switch(bound_x)
	    {
	    case -2: // reflective boundary conditions
		bfv_L[i].SU = - CV->s_u[0][i]; bfv_R[i].SU = - CV->s_u[m-1][i];
		break;
	    case -5: // periodic boundary conditions
		bfv_L[i].SU   =   CV->s_u[m-1][i]; bfv_R[i].SU   =   CV->s_u[0][i];
		bfv_L[i].SV   =   CV->s_v[m-1][i]; bfv_R[i].SV   =   CV->s_v[0][i];
		bfv_L[i].SP   =   CV->s_p[m-1][i]; bfv_R[i].SP   =   CV->s_p[0][i];
		bfv_L[i].SRHO = CV->s_rho[m-1][i]; bfv_R[i].SRHO = CV->s_rho[0][i];
		bfv_L[i].TU   =   CV->t_u[m-1][i]; bfv_R[i].TU   =   CV->t_u[0][i];
		bfv_L[i].TV   =   CV->t_v[m-1][i]; bfv_R[i].TV   =   CV->t_v[0][i];
		bfv_L[i].TP   =   CV->t_p[m-1][i]; bfv_R[i].TP   =   CV->t_p[0][i];
		bfv_L[i].TRHO = CV->t_rho[m-1][i]; bfv_R[i].TRHO = CV->t_rho[0][i];
		break;
	    case -24: // reflective + free boundary conditions
		bfv_L[i].SU = - CV->s_u[0][i];
		break;
	    }
    return true;
}
