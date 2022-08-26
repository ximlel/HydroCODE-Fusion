#include <stdio.h>
#include <stdarg.h>

#include "../include/var_struc.h"
#include "../include/tools.h"


void minmod_limiter_2D_x(const _Bool NO_h, const int m, const int i, const _Bool find_bound_x, double ** s,
			 double * U[], const double UL, const double UR, const double HL, ...)
{
    va_list ap;
    va_start(ap, HL);
    int j;
    double const alpha = config[41]; // the paramater in slope limiters.
    double s_L, s_R; // spatial derivatives in coordinate x (slopes) 
    double h = HL, HR, * X;
    if (NO_h)
	{
	    HR = va_arg(ap, double);
	    X  = va_arg(ap, double *);
	}
	
    for(j = 0; j < m; ++j) // Reconstruct slopes
	{ /*
	   *  j-1          j          j+1
	   * j-1/2  j-1  j+1/2   j   j+3/2  j+1
	   *   o-----X-----o-----X-----o-----X--...
	   */
	    if(j)
		{
		    if (NO_h)
			h = 0.5 * (X[j+1] - X[j-1]);
		    s_L = (U[j][i] - U[j-1][i]) / h;
		}
	    else
		{
		    if (NO_h)
			h = 0.5 * (X[j+1] - X[j] + HL);
		    s_L = (U[j][i] - UL) / h;
		}
	    if(j < m-1)
		{
		    if (NO_h)
			h = 0.5 * (X[j+2] - X[j]);
		    s_R = (U[j+1][i] - U[j][i]) / h;
		}
	    else
		{
		    if (NO_h)
			h = 0.5 * (X[j+1] - X[j] + HR);
		    s_R = (UR - U[j][i]) / h;
		}
	    if (find_bound_x)
		s[j][i] = minmod3(alpha*s_L, alpha*s_R, s[j][i]);
	    else
		s[j][i] = minmod2(s_L, s_R);
	}
    va_end(ap);
}
