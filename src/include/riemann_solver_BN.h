/**
 * @file riemann_solver_BN.h
 * @brief This file is the header file of several Riemann solvers and GRP solvers
 *        for Baer-Nunziato two-phase flow model.
 * @details This header file declares functions in the folder 'Riemann_solver_BN'.
 */

#ifndef RIEMANNSOLVER_BN_H
#define RIEMANNSOLVER_BN_H

#include "../include/var_struc.h"

int Roe_GRP_solver_BN
(double *Dt_U_all, double *U_all,
 const double   z_g_L, const double   z_g_R, const double   d_z_g_L, const double   d_z_g_R, const double   t_z_g_L, const double   t_z_g_R,
 const double rho_g_L, const double rho_g_R, const double d_rho_g_L, const double d_rho_g_R, const double t_rho_g_L, const double t_rho_g_R,
 const double   u_g_L, const double   u_g_R, const double   d_u_g_L, const double   d_u_g_R, const double   t_u_g_L, const double   t_u_g_R,
 const double   v_g_L, const double   v_g_R, const double   d_v_g_L, const double   d_v_g_R, const double   t_v_g_L, const double   t_v_g_R,
 const double   p_g_L, const double   p_g_R, const double   d_p_g_L, const double   d_p_g_R, const double   t_p_g_L, const double   t_p_g_R,
 const double rho_l_L, const double rho_l_R, const double d_rho_l_L, const double d_rho_l_R, const double t_rho_l_L, const double t_rho_l_R,
 const double   u_l_L, const double   u_l_R, const double   d_u_l_L, const double   d_u_l_R, const double   t_u_l_L, const double   t_u_l_R,
 const double   v_l_L, const double   v_l_R, const double   d_v_l_L, const double   d_v_l_R, const double   t_v_l_L, const double   t_v_l_R,
 const double   p_l_L, const double   p_l_R, const double   d_p_l_L, const double   d_p_l_R, const double   t_p_l_L, const double   t_p_l_R,
 const double gamma_g, const double gamma_l, const double  eps);

int Roe_GRP_solver_BN_1D
(double *Dt_U_all, double *U_all,
 const double   z_g_L, const double   z_g_R, const double   d_z_g_L, const double   d_z_g_R, const double   t_z_g_L, const double   t_z_g_R,
 const double rho_g_L, const double rho_g_R, const double d_rho_g_L, const double d_rho_g_R, const double t_rho_g_L, const double t_rho_g_R,
 const double   u_g_L, const double   u_g_R, const double   d_u_g_L, const double   d_u_g_R, const double   t_u_g_L, const double   t_u_g_R,
 const double   v_g_L, const double   v_g_R, const double   d_v_g_L, const double   d_v_g_R, const double   t_v_g_L, const double   t_v_g_R,
 const double   p_g_L, const double   p_g_R, const double   d_p_g_L, const double   d_p_g_R, const double   t_p_g_L, const double   t_p_g_R,
 const double rho_l_L, const double rho_l_R, const double d_rho_l_L, const double d_rho_l_R, const double t_rho_l_L, const double t_rho_l_R,
 const double   u_l_L, const double   u_l_R, const double   d_u_l_L, const double   d_u_l_R, const double   t_u_l_L, const double   t_u_l_R,
 const double   v_l_L, const double   v_l_R, const double   d_v_l_L, const double   d_v_l_R, const double   t_v_l_L, const double   t_v_l_R,
 const double   p_l_L, const double   p_l_R, const double   d_p_l_L, const double   d_p_l_R, const double   t_p_l_L, const double   t_p_l_R,
 const double gamma_g, const double gamma_l, const double  eps);

int linear_GRP_RI_solver_BN
(struct RI_var *RI, const double D_z_s, const double z_s, const double *mid_g, const double *mid_s, 
 const struct GRP_LR_var GL, const struct GRP_LR_var GR,
 const double gamma_s, const double gamma_g, const double eps, const double tau, const int x_or_y);


#endif
