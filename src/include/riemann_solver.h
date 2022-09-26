/**
 * @file riemann_solver.h
 * @brief This file is the header file of several Riemann solvers and GRP solvers.
 * @details This header file declares functions in the folder 'Riemann_solver'.
 */

#ifndef RIEMANNSOLVER_H
#define RIEMANNSOLVER_H

#include "../include/var_struc.h"


/**
 * @brief Which solver is chosen as the exact Riemann solver for single-component flow.
 */
#ifndef Riemann_solver_exact_single
#define Riemann_solver_exact_single Riemann_solver_exact_Ben
#endif

/* exact Riemann solver (two-component flow) */
//////////////////////////////////////
// riemann_solver_exact_Ben.c
//////////////////////////////////////
double Riemann_solver_exact(double * U_star, double * P_star, const double gammaL, const double gammaR,
			    const double u_L, const double u_R, const double p_L, const double p_R, 
			    const double c_L, const double c_R, _Bool * CRW,
			    const double eps, const double tol, int N);
//////////////////////////////////////
// riemann_solver_starPU.c
//////////////////////////////////////
double Riemann_solver_starPU(double * U_star, double * P_star, const double GammaL, const double GammaR,
			   const double UL, const double UR, const double PL, const double PR,
			   const double CL, const double CR, _Bool * CRW,
			   const double eps, const double TOLPRE, const int NRITER);
/* exact Riemann solver (single-component flow) */
//////////////////////////////////////
// riemann_solver_exact_Ben.c
//////////////////////////////////////
double Riemann_solver_exact_Ben(double * U_star, double * P_star, const double gamma,
				const double u_L, const double u_R, const double p_L, const double p_R,
				const double c_L, const double c_R, _Bool * CRW,
				const double eps, const double tol, const int N);
//////////////////////////////////////
// riemann_solver_exact_Toro.c
//////////////////////////////////////
double Riemann_solver_exact_Toro(double * U_star, double * P_star, const double gamma,
				 const double U_l, const double U_r, const double P_l, const double P_r,
				 const double c_l, const double c_r, _Bool * CRW,
				 const double eps, const double tol, const int N);

//////////////////////////////////////
// linear_grp_solver_LAG.c
//////////////////////////////////////
/* 1-D GRP solver (Lagrangian, two-component flow) */
void linear_GRP_solver_LAG (double *D, double *U, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R, const double eps, const double  atc);
//////////////////////////////////////
// linear_grp_solver_Edir.c
//////////////////////////////////////
/* 1-D GRP solver (Eulerian, single-component flow) */
void linear_GRP_solver_Edir(double *D, double *U, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R, const double eps, const double  atc);

/* 2-D GRP solver (ALE, two-component flow) */
//////////////////////////////////////
// linear_grp_solver_Edir_Q1D.c
//////////////////////////////////////
void linear_GRP_solver_Edir_Q1D(double *wave_speed, double *D, double *U, double *U_star, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R, const double  eps, const double  atc);
//////////////////////////////////////
// linear_grp_solver_Edir_G2D.c
//////////////////////////////////////
void linear_GRP_solver_Edir_G2D(double *wave_speed, double *D, double *U, double *U_star, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R, const double  eps, const double  atc);

//////////////////////////////////////
// linear_grp_solver_radial_LAG.c
//////////////////////////////////////
void AcousticRLagTangent(double *dire, double *U_star, const struct i_f_var * ifv_L, const struct i_f_var * ifv_R,
			 double r, double M, const double eps);
void GRPsolverRLag(double *wave_speed, double *dire, double *U_star, const struct i_f_var * ifv_L, const struct i_f_var * ifv_R,
		   double r, double M, const double eps, const double  atc);


/* HLL solver (single-component flow) */
//////////////////////////////////////
// hll_2D_solver.c
//////////////////////////////////////
void HLL_2D_solver(double *F, double *lambda_max, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R);

/* Roe solver (single-component flow) */
//////////////////////////////////////
// roe_solver.c
//////////////////////////////////////
void Roe_solver   (double *F, double *lambda_max, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R, const double delta);
//////////////////////////////////////
// roe_2D_solver.c
//////////////////////////////////////
void Roe_2D_solver(double *F, double *lambda_max, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R, const double delta);
//////////////////////////////////////
// roe_HLL_solver.c
//////////////////////////////////////
void Roe_HLL_solver(double *V_mk, double *F, double *lambda_max, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R, const double delta);

#endif
