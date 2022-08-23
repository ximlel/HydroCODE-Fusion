/**
 * @file flux_calc.h
 * @brief This file is the header file of intermediate processes of finite volume scheme.
 * @details This header file declares functions in the folder 'flux_calc'.
 */

#ifndef FLUXCALC_H
#define FLUXCALC_H

#include "../include/var_struc.h"

void flux_generator_x(int m, int n, int nt, double tau, struct cell_var_stru * CV,
		      struct b_f_var * bfv_L, struct b_f_var * bfv_R);
void flux_generator_y(int m, int n, int nt, double tau, struct cell_var_stru * CV,
		      struct b_f_var * bfv_D, struct b_f_var * bfv_U);

void GRP_2D_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, const double tau);

#endif
