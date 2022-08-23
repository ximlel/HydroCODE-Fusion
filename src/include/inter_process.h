    /**
 * @file inter_process.h
 * @brief This file is the header file of intermediate processes of finite volume scheme.
 * @details This header file declares functions in the folder 'inter_process'.
 */

#ifndef INTERPROCESS_H
#define INTERPROCESS_H

#include "../include/var_struc.h"

void minmod_limiter(_Bool NO_h, int m, _Bool find_bound_x, double * s, double * U, double UL, double UR, double HL, ...);
void minmod_limiter_2D_x(_Bool NO_h, int m, int i, _Bool find_bound_x, double ** s, double ** U, double UL, double UR, double HL, ...);

_Bool bound_cond_slope_limiter_x(const int m, const int n, const int nt, struct cell_var_stru * CV, struct b_f_var  * bfv_L, struct b_f_var  * bfv_R, _Bool find_bound_x);
_Bool bound_cond_slope_limiter_y(const int m, const int n, const int nt, struct cell_var_stru * CV, struct b_f_var  * bfv_D, struct b_f_var  * bfv_U, _Bool find_bound_y);

#endif
