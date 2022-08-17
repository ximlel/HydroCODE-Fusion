/**
 * @file tools.h
 * @brief This file is the header file of several independent tool functions.
 * @details This header file declares functions in files 'sys_pro.c' and 'math_algo.c',
 */

#ifndef TOOLS_H
#define TOOLS_H

void DispPro(double pro, int step);
int CreateDir(const char* pPath);

int rinv(double a[], int n);
void Gauss_elimination(int n, double (*a)[n+1], double *x);

/**
 * @brief Minmod limiter of two variables.
 */
inline double minmod2(double s_L, double s_R)
{
    if(s_L * s_R < 0.0)
	return 0.0;
    else if(fabs(s_R) < fabs(s_L))
	return s_R;
    else
	return s_L;
}

/**
 * @brief Minmod limiter of three variables.
 */
inline double minmod3(double s_L, double s_R, double s_m)
{
    if(s_L * s_m < 0.0 || s_R * s_m < 0.0)
	return 0.0;
    else if(fabs(s_m) < fabs(s_L) && fabs(s_m) < fabs(s_R))
	return s_m;
    else if(fabs(s_R) < fabs(s_L))
	return s_R;
    else
	return s_L;
}

#endif
