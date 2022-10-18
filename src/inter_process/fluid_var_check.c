/**
 * @file  fluid_var_check.c
 * @brief There are some functions to check whether fluid variables are within the value range.
 */
#include <stdio.h>
#include <math.h>

#include "../include/var_struc.h"


/**
 * @brief This function checks whether interfacial fluid variables are within the value range.
 * @param[in] ifv_L: Structure pointer of interfacial left state.
 * @param[in] ifv_R: Structure pointer of interfacial right state.
 * @param[in] dim:   Spatial dimension.
 * @return    miscalculation indicator.
 *   @retval  0: Successful calculation.
 *   @retval  1: < 0.0 error.
 *   @retval  2: NAN or INFinite error of Slope.
 */
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


/**
 * @brief This function checks whether fluid variables of mid[] and dire[] are within the value range.
 * @param[in] mid:  Intermediate Riemann solutions at t-axis OR in star region.
 * @param[in] dire: Temporal derivative of fluid variables.
 * @param[in] dim:  Spatial dimension.
 * @return    miscalculation indicator.
 *   @retval  0: Successful calculation.
 *   @retval  1: < 0.0 error of mid[].
 *   @retval  2: NAN or INFinite error of mid[].
 *   @retval  3: NAN or INFinite error of dire[].
 */
int star_dire_check(double *mid, double *dire, const int dim)
{
    double const eps = config[4];
    int    const el  = (int)config[8];
    double * star = NULL;
    if (dim == 1)
	{
	    switch(el)
		{
		case 1:
			star = mid;
		    if(star[2] < eps || star[0] < eps || star[3] < eps)
			{
			    printf("<0.0 error - STAR");
			    return 1;
			}
		    if(!isfinite(star[1])|| !isfinite(star[2])|| !isfinite(star[0])|| !isfinite(star[3]))
			{
			    printf("NAN or INFinite error - STAR"); 
			    return 2;
			}
		    if(!isfinite(dire[1])|| !isfinite(dire[2])|| !isfinite(dire[0])|| !isfinite(dire[3]))
			{
			    printf("NAN or INFinite error - DIRE"); 
			    return 3;
			}
		    break;
		default:
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
		    break;
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
