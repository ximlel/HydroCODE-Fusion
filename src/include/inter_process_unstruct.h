/**
 * @file inter_process_unstruct.h
 * @brief This file is the header file of intermediate processes of finite volume scheme on unstructured grid.
 * @details This header file declares functions in the folder 'inter_process_unstruct'.
 */

#ifndef INTERPROCESS_UNSTRUCT_H
#define INTERPROCESS_UNSTRUCT_H

#include "../include/var_struc.h"


/////////////////////////
// cons_qty_calc.c
/////////////////////////
void cons_qty_init(const struct cell_var * cv, const struct flu_var * FV);
int cons2prim(struct i_f_var * ifv);
int cons_qty_update(const struct cell_var * cv, const struct mesh_var * mv,
					const struct flu_var *  FV, const double tau);
/////////////////////////
// cons_qty_update_P_ave.c
/////////////////////////
int cons_qty_update_corr_ave_P(struct cell_var * cv, const struct mesh_var * mv,
							   const struct flu_var * FV, double tau, const int RK);

/////////////////////////
// cell_init_free.c
/////////////////////////
void cell_mem_init_free(struct cell_var * cv, const struct mesh_var * mv, struct flu_var * FV, const int i_or_f);
void vol_comp(const struct cell_var * cv, const struct mesh_var * mv);
void cell_pt_clockwise(const struct mesh_var * mv);
void cell_rel(const struct cell_var * cv, const struct mesh_var * mv);
void cell_centroid(const struct cell_var * cv, const struct mesh_var * mv);

/////////////////////////
// slope_limiter_unstruct.c
/////////////////////////
void slope_limiter_prim(const struct cell_var * cv,const struct mesh_var * mv, const struct flu_var * FV);

/////////////////////////
// copy_func.c
/////////////////////////
void cons_qty_copy_cv2ifv(struct i_f_var * ifv, const struct cell_var * cv, const int c);
void cons_qty_copy_ifv2cv(const struct i_f_var * ifv, struct cell_var * cv, const int c);
void prim_var_copy_ifv2FV(const struct i_f_var * ifv, const struct flu_var * FV,const int c);
void flux_copy_ifv2cv(const struct i_f_var * ifv, const struct cell_var *cv, const int k, const int j);
void flux_add_ifv2cv(const struct i_f_var * ifv, const struct cell_var * cv, const int k, const int j);

/////////////////////////
// assist_func.c
/////////////////////////
int fluid_var_update(struct flu_var *FV, struct cell_var *cv);
int interface_var_init(const struct cell_var * cv, const struct mesh_var * mv,
					   struct i_f_var * ifv, struct i_f_var * ifv_R,
					   const int k, const int j, const int i, const double gauss);
double tau_calc(const struct cell_var * cv, const struct mesh_var * mv);

#endif
