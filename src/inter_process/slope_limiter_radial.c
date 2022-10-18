/**
 * @file  slope_limiter_radial.c
 * @brief This is a function of the minmod slope limiter in radially symmetric case.
 */
#include <stdio.h>
#include <stdlib.h>

#include "../include/var_struc.h"
#include "../include/tools.h"


/**
 * @brief This function apply the minmod limiter to the slope in radially symmetric case.
 * @param[in] Ncell:      Number of the r-grids.
 * @param[in] i_f_var_get: Whether the cell interfacial variables have been obtained.
 *                        - true: interfacial variables at t_{n+1} are available, 
 *                                and then trivariate minmod3() function is used.
 *                        - false: bivariate minmod2() function is used.
 * @param[in,out] s[]:    Spatial derivatives of the fluid variable are stored here.
 * @param[in] U[]: Array to store fluid variable values.
 * @param[in] rmv: Structure of radially symmetric meshing variable data.
 */
void minmod_limiter_radial(const int Ncell, const _Bool i_f_var_get, double s[],
			  const double U[], struct radial_mesh_var *rmv)
{
    double const Alpha       =      config[41]; // the paramater in slope limiters.
    int    const LIMITER_VIP = (int)config[42];
    double s_L, s_R; // spatial derivatives in coordinate x (slopes)

    //minmod limiter update
    for(int j = 1; j <= Ncell; ++j) // Reconstruct slopes
	{
	    if (abs(LIMITER_VIP)==1)
		{
		    s_L = (U[j]   - U[j-1]) / rmv->dRc[j];
		    s_R = (U[j+1] - U[j])   / rmv->dRc[j+1];
		}
	    else if (abs(LIMITER_VIP)==2)
		{
		    s_L = (U[j]   - U[j-1]) / 2.0 / rmv->DdrR[j];
		    s_R = (U[j+1] - U[j])   / 2.0 / rmv->DdrL[j];
		}
	    else
		{
		    fprintf(stderr, "ERROE! No suitable LIMITER_VIP Parameter.\n");
		    exit(2);
		}
	    if (i_f_var_get)
		s[j] = minmod3(Alpha*s_L, Alpha*s_R, s[j]);
	    else
		s[j] = minmod2(s_L, s_R);
	}
    s[0] = minmod2(s[0], s[1]);
    s[Ncell+1] = minmod2(s[Ncell], s[Ncell+1]);
}
