#include <stdio.h>
#include <math.h>

#include "../include/var_struc.h"


int ifvar_check(struct i_f_var *ifv_L, struct i_f_var *ifv_R, const int dim)
{
    double const eps = config[4];
    if(ifv_L->P < eps || ifv_R->P < eps || ifv_L->RHO < eps || ifv_R->RHO < eps)
	{
	    printf("<0.0 error - Reconstruction");
	    return 1;
	}
    if(dim == 1)
	{
	    if(!isfinite(ifv_L->d_p)|| !isfinite(ifv_R->d_p)|| !isfinite(ifv_L->d_u)|| !isfinite(ifv_R->d_u)|| !isfinite(ifv_L->d_rho)|| !isfinite(ifv_R->d_rho))
		{
		    printf("NAN or INFinite error - Slope"); 
		    return 2;
		}
	}
    else if (dim == 2)
	{
	    if(!isfinite(ifv_L->d_p)|| !isfinite(ifv_R->d_p)|| !isfinite(ifv_L->d_u)|| !isfinite(ifv_R->d_u)|| !isfinite(ifv_L->d_v)|| !isfinite(ifv_R->d_v)|| !isfinite(ifv_L->d_rho)|| !isfinite(ifv_R->d_rho))
		{
		    printf("NAN or INFinite error - d_Slope_x"); 
		    return 2;
		}
	    if(!isfinite(ifv_L->t_p)|| !isfinite(ifv_R->t_p)|| !isfinite(ifv_L->t_u)|| !isfinite(ifv_R->t_u)|| !isfinite(ifv_L->t_v)|| !isfinite(ifv_R->t_v)|| !isfinite(ifv_L->t_rho)|| !isfinite(ifv_R->t_rho))
		{
		    printf("NAN or INFinite error - t_Slope_x"); 
		    return 2;
		}
	}
    return 0;
}

int star_dire_check(double *mid, double *dire, const int dim)
{
    double const eps = config[4];
    if (dim == 1)
	{
	    if(mid[2] < eps || mid[0] < eps)
		{
		    printf("<0.0 error - STAR");
		    return 1;
		}
	    if(!isfinite(mid[1])|| !isfinite(mid[2])|| !isfinite(mid[0]))
		{
		    printf("NAN or INFinite error - STAR"); 
		    return 2;
		}
	    if(!isfinite(dire[1])|| !isfinite(dire[2])|| !isfinite(dire[0]))
		{
		    printf("NAN or INFinite error - DIRE"); 
		    return 3;
		}
	}
    else if(dim == 2)
	{
	    if(mid[3] < eps || mid[0] < eps)
		{
		    printf("<0.0 error - STAR");
		    return 1;
		}
	    if(!isfinite(mid[1])|| !isfinite(mid[2])|| !isfinite(mid[0])|| !isfinite(mid[3]))
		{
		    printf("NAN or INFinite error - STAR"); 
		    return 2;
		}
	    if(!isfinite(dire[1])|| !isfinite(dire[2])|| !isfinite(dire[0])|| !isfinite(dire[3]))
		{
		    printf("NAN or INFinite error - DIRE"); 
		    return 3;
		}
	}

    return 0;
}
