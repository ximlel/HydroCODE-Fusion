#include <stdio.h>
#include <stdlib.h>

#include "../include/var_struc.h"
#include "../include/tools.h"


void minmod_limiter_spher(const int Ncell, const _Bool find_bound, double s[],
			  const double U[], struct spher_mesh_var *smv)
{
    double const Alpha       =      config[41]; // the paramater in slope limiters.
    int    const LIMITER_VIP = (int)config[42];
    double s_L, s_R; // spatial derivatives in coordinate x (slopes)

    //minmod limiter update
    for(int j = 1; j <= Ncell; ++j) // Reconstruct slopes
	{
	    if (abs(LIMITER_VIP)==1)
		{
		    s_L = (U[j]   - U[j-1]) / smv->dRc[j];
		    s_R = (U[j+1] - U[j])   / smv->dRc[j+1];
		}
	    else if (abs(LIMITER_VIP)==2)
		{
		    s_L = (U[j]   - U[j-1]) / 2.0 / smv->DdrR[j];
		    s_R = (U[j+1] - U[j])   / 2.0 / smv->DdrL[j];
		}
	    else
		{
		    printf("LIMITER_VIP Parameter Error!\n");
		    exit(2);
		}
	    if (find_bound)
		s[j] = minmod3(Alpha*s_L, Alpha*s_R, s[j]);
	    else
		s[j] = minmod2(s_L, s_R);
	}
    s[0] = minmod2(s[0], s[1]);
    s[Ncell+1] = minmod2(s[Ncell], s[Ncell+1]);
}
