/**
 * @file  slope_limiter.c
 * @brief This is a function of the minmod slope limiter in one dimension.
 */
#include <stdio.h>
#include <stdarg.h>

#include "../include/var_struc.h"
#include "../include/tools.h"


/**
 * @brief This function apply the minmod limiter to the slope in one dimension.
 * @param[in] NO_h:       Whether there are moving grid point coordinates.
 *                  - true: There are moving spatial grid point coordinates *X.
 *                  - false: There is fixed spatial grid length.
 * @param[in] m:         Number of the x-grids: n_x.
 * @param[in] i_f_var_get: Whether the cell interfacial variables have been obtained.
 *                        - true: interfacial variables at t_{n+1} are available, 
 *                                and then trivariate minmod3() function is used.
 *                        - false: bivariate minmod2() function is used.
 * @param[in,out] s[]:    Spatial derivatives of the fluid variable are stored here.
 * @param[in] U[]: Array to store fluid variable values.
 * @param[in] UL:  Fluid variable value at left boundary.
 * @param[in] UR:  Fluid variable value at right boundary.
 * @param[in] HL:  Spatial grid length at left boundary OR fixed spatial grid length.
 * @param[in] ...: Variable parameter if NO_h is true.
 *            - \b double \c HR: Spatial grid length at right boundary.
 *            - \b double \c *X: Array of moving spatial grid point coordinates.
 */
void minmod_limiter(const _Bool NO_h, const int m, const _Bool i_f_var_get, double s[],
		    const double U[], const double UL, const double UR, const double HL, ...)
{
    va_list ap;
    va_start(ap, HL);
    double const alpha = config[41]; // the paramater in slope limiters.
    double s_L, s_R; // spatial derivatives in coordinate x (slopes) 
    double h = HL, HR, * X;
    if (NO_h)
	{
	    HR = va_arg(ap, double);
	    X  = va_arg(ap, double *);
	}
#ifdef _OPENACC
#pragma acc parallel loop private(s_L, s_R, h)
#endif
    for(int j = 0; j < m; ++j) // Reconstruct slopes
	{ /*
	   *  j-1          j          j+1
	   * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	   *   o-----X-----o-----X-----o-----X--...
	   */
	    if(j)
		{
		    if (NO_h)
			h = 0.5 * (X[j+1] - X[j-1]);
		    s_L = (U[j] - U[j-1]) / h;
		}
	    else
		{
		    if (NO_h)
			h = 0.5 * (X[j+1] - X[j] + HL);
		    s_L = (U[j] - UL) / h;
		}
	    if(j < m-1)
		{
		    if (NO_h)
			h = 0.5 * (X[j+2] - X[j]);
		    s_R = (U[j+1] - U[j]) / h;
		}
	    else
		{
		    if (NO_h)
			h = 0.5 * (X[j+1] - X[j] + HR);
		    s_R = (UR - U[j]) / h;
		}
	    if (i_f_var_get)
		s[j] = minmod3(alpha*s_L, alpha*s_R, s[j]);
	    else
		s[j] = minmod2(s_L, s_R);
	} // End of parallel region
    va_end(ap);
}
