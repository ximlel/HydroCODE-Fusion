    /**
 * @file inter_process.h
 * @brief This file is the header file of intermediate processes of finite volume scheme.
 * @details This header file declares functions in the folder 'inter_process'.
 */

#ifndef INTERPROCESS_H
#define INTERPROCESS_H

#include "../include/var_struc.h"


///////////////////////////////////
// fluid_var_check.c
///////////////////////////////////
int ifvar_check(struct i_f_var *ifv_L, struct i_f_var *ifv_R, const int dim);
int star_dire_check(double *mid, double *dire, const int dim);


///////////////////////////////////
// slope_limiter.c
///////////////////////////////////
void minmod_limiter(const _Bool NO_h, const int m, const _Bool i_f_var_get, double s[],
		    const double U[], const double UL, const double UR, const double HL, ...);
///////////////////////////////////
// slope_limiter_2D_x.c
///////////////////////////////////
void minmod_limiter_2D_x(const _Bool NO_h, const int m, const int i, const _Bool i_f_var_x_get, double ** s,
			 double ** U, const double UL, const double UR, const double HL, ...);
///////////////////////////////////
// slope_limiter_radial.c
///////////////////////////////////
void minmod_limiter_radial(const int Ncell, const _Bool i_f_var_get, double s[],
			  const double U[], struct radial_mesh_var *rmv);
///////////////////////////////////
// slope_VIP_limiter_radial.c
///////////////////////////////////
void VIP_limiter_radial(const int Ncell, const _Bool i_f_var_get, double DmU[], double TmV[],
		       const double UU[], struct radial_mesh_var *rmv);


/* Set boundary conditions & Use the slope limiter */
///////////////////////////////////
// bound_cond_slope_limiter.c
///////////////////////////////////
_Bool bound_cond_slope_limiter(const _Bool NO_h, const int m, const int nt, struct cell_var_stru * CV,
			       struct b_f_var * bfv_L, struct b_f_var * bfv_R, _Bool find_bound, const _Bool Slope, const double t_c, ...);
///////////////////////////////////
// bound_cond_slope_limiter_x.c
///////////////////////////////////
_Bool bound_cond_slope_limiter_x(const int m, const int n, const int nt, struct cell_var_stru * CV, struct b_f_var * bfv_L, struct b_f_var * bfv_R,
				 struct b_f_var * bfv_D, struct b_f_var * bfv_U, _Bool find_bound_x, const _Bool Slope, const double t_c);
///////////////////////////////////
// bound_cond_slope_limiter_y.c
///////////////////////////////////
_Bool bound_cond_slope_limiter_y(const int m, const int n, const int nt, struct cell_var_stru * CV, struct b_f_var * bfv_L, struct b_f_var * bfv_R,
				 struct b_f_var * bfv_D, struct b_f_var * bfv_U, _Bool find_bound_y, const _Bool Slope, const double t_c);

#endif
