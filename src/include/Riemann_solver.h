/**
 * @file Riemann_solver.h
 * @brief This file is the header file of several Riemann solvers and GRP solvers.
 * @details This header file declares functions in the folder 'Riemann_solver'.
 */

#ifndef RIEMANNSOLVER_H
#define RIEMANNSOLVER_H

double Riemann_solver_exact(double * U_star, double * P_star, const double gammaL, const double gammaR,
			    const double u_L, const double u_R, const double p_L, const double p_R, 
			    const double c_L, const double c_R, _Bool * CRW,
			    const double eps, const double tol, int N);
double Riemann_solver_exact_Ben(double * U_star, double * P_star, const double gamma,
				const double u_L, const double u_R, const double p_L, const double p_R,
				const double c_L, const double c_R, _Bool * CRW,
				const double eps, const double tol, const int N);
double Riemann_solver_exact_Toro(double * U_star, double * P_star, const double gamma,
				 const double U_l, const double U_r, const double P_l, const double P_r,
				 const double c_l, const double c_r, _Bool * CRW,
				 const double eps, const double tol, const int N);

void linear_GRP_solver_LAG(double * D, double * U, const struct i_f_var ifv_L, const struct i_f_var ifv_R, const double eps, const double  atc);
void linear_GRP_solver_LAG(double * D, double * U, const struct i_f_var ifv_L, const struct i_f_var ifv_R, const double eps, const double  atc);
void linear_GRP_solver_Edir(double * D, double * U, const struct i_f_var ifv_L, const struct i_f_var ifv_R, const double eps);

void linear_GRP_solver_Edir_Q1D(double *wave_speed, double *D, double *U, double *U_star, const struct i_f_var ifv_L, const struct i_f_var ifv_R, const double  eps, const double  atc);
void linear_GRP_solver_Edir_G2D(double *wave_speed, double *D, double *U, double *U_star, const struct i_f_var ifv_L, const struct i_f_var ifv_R, const double  eps, const double  atc);

#ifndef Riemann_solver_exact_single
#define Riemann_solver_exact_single Riemann_solver_exact_Ben
#endif

#endif
