/**
 * @file tools.h
 * @brief This file is the header file of several independent tool functions.
 * @details This header file declares functions in the folder 'tools',
 */

#ifndef TOOLS_H
#define TOOLS_H

// sys_pro.c
void DispPro(const double pro, const int step);

int CreateDir(const char* pPath);

void init_mem (double * p[], const int n, int ** cell_pt);
void init_mem_int(int * p[], const int n, int ** cell_pt);


// math_algo.c
void mat_mul(const double A[], const double B[], double C[], const int m, const int p, const int n);
void mat_add(const double A[], const double B[], double C[], const int m,              const int n);
void mat_sub(const double A[], const double B[], double C[], const int m,              const int n);

int rinv(double a[], const int n);


/**
 * @brief Minmod limiter function of two variables.
 */
inline double minmod2(const double s_L, const double s_R)
{
    if(s_L * s_R <= 0.0)
	return 0.0;
    else if(s_R >  0.0 && s_R < s_L)
	return s_R;
    else if(s_R <= 0.0 && s_R > s_L)
	return s_R;
    else // fabs(s_R) > fabs(s_L)
	return s_L;
}

/**
 * @brief Minmod limiter function of three variables.
 */
inline double minmod3(const double s_L, const double s_R, const double s_m)
{
    if(s_L * s_m <= 0.0 || s_R * s_m <= 0.0)
	return 0.0;
    else if(s_m >  0.0 && s_m < s_L && s_m < s_R)
	return s_m;
    else if(s_m <= 0.0 && s_m > s_L && s_m > s_R)
	return s_m;
    else if(s_R >  0.0 && s_R < s_L)
	return s_R;
    else if(s_R <= 0.0 && s_R > s_L)
	return s_R;
    else
	return s_L;
}

#endif
