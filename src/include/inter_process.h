/**
 * @file inter_process.h
 * @brief This file is the header file of intermediate processes of finite volume scheme.
 * @details This header file declares functions in the folder 'inter_process'.
 */

#ifndef INTERPROCESS_H
#define INTERPROCESS_H

#include "../include/var_struc.h"

void minmod_limiter(_Bool NO_h, int m, _Bool find_bound_x, double * s, double * U, double UL, double UR, ...);


_Bool bound_cond_slope_limiter_x(const int m, const int n, const int nt, S_CVS * CV, S_BFV * bfv_L, S_BFV * bfv_R, _Bool find_bound_x);
_Bool bound_cond_slope_limiter_y(const int m, const int n, const int nt, S_CVS * CV, S_BFV * bfv_U, S_BFV * bfv_D, _Bool find_bound_y);

#endif
