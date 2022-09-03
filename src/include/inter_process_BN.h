    /**
 * @file inter_process_BN.h
 * @brief This file is the header file of intermediate processes of finite volume scheme
 *        for Baer-Nunziato two-phase flow model.
 * @details This header file declares functions in the folder 'inter_process_BN'.
 */

#ifndef INTERPROCESS_BN_H
#define INTERPROCESS_BN_H

#include "../include/var_struc.h"


void FV_2_C_init(struct center_var C, struct flu_var FV);
void NewtonRapshon(double * x_star, double * err, double fun, double dfun, double eps);
void NewtonRapshon_matrix(double * x_star, double * err, double * fun, double * dfun, double eps);
void RI2U_cal(struct U_var * U, const struct RI_var * RI, double z_s, const double rho_g_start);
void U2RI_cal(const struct U_var * U, struct RI_var * RI);
void primitive_comp(double * U, struct U_var * U_L, struct U_var * U_R, double z_sL, double z_sR, double z_sL_out, double z_sR_out, double area_L, double area_R);


void BN_C2U(struct center_var C, double *U, int i, int j, int x_or_y);
void BN_ULR2prim(struct U_var U_L, struct U_var U_R, struct center_var C, int i, int j, int x_or_y);
void BN_ULR2cons(struct U_var U_L, struct U_var U_R, struct center_var C, int i, int j, int x_or_y);
void RI_LR_ave(struct RI_var *RI, struct RI_var RI_L,struct RI_var RI_R);
void BN_RI2Cx(struct RI_var RI, struct center_var C, int i, int j);
void BN_RI2Cy(struct RI_var RI, struct center_var C, int i, int j);
void GRP_var_init(struct GRP_LR_var *G, struct slope_var SV, struct U_var U, double d, int i, int j, int pm_xy);
void GRP_RI_var_init(struct GRP_LR_var *G, struct slope_var SV, struct center_var C, double d, int i, int j, int pm_xy);
void G_LR_RI2U(struct GRP_LR_var *G, double z_s, int x_or_y);
void boundary_cond_x(struct center_var C, int cond, int l);
void boundary_cond_y(struct center_var C, int cond, int l);
void boundary_cond_slope_x(struct slope_var SV, int cond, int l);
void boundary_cond_slope_y(struct slope_var SV, int cond, int l);

#endif
