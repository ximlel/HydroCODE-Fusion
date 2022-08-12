/**
 * @file Riemann_solver.h
 * @brief This file is the header file of several Riemann solvers and GRP solvers.
 * @details This header file declares functions in files 'Riemann_solver_exact_Ben.c',
 *          'Riemann_solver_exact_Toro.c',  'linear_GRP_solver_LAG.c' and 'linear_GRP_solver_Edir.c'.
 */

#ifndef RIEMANNSOLVER_H
#define RIEMANNSOLVER_H

double Riemann_solver_exact_Ben(double * U_star, double * P_star, const double gamma,
			    const double u_L, const double u_R, const double p_L, const double p_R,
			    const double c_L, const double c_R, _Bool * CRW,
			    const double eps, const double tol, const int N);

double Riemann_solver_exact_Toro(double * U_star, double * P_star, const double gamma,
				 const double U_l, const double U_r, const double P_l, const double P_r,
				 const double c_l, const double c_r, _Bool * CRW,
				 const double eps, const double tol, const int N);

void linear_GRP_solver_LAG
(double * dire, double * mid,
 const double rho_L, const double rho_R, const double s_rho_L, const double s_rho_R,
 const double   u_L, const double   u_R, const double   s_u_L, const double   s_u_R,
 const double   p_L, const double   p_R, const double   s_p_L, const double   s_p_R,
 const double gamma, const double eps, const double  atc);

void linear_GRP_solver_Edir
(double * direvative, double * mid,
 const double rho_L, const double rho_R, const double s_rho_L, const double s_rho_R,
 const double   u_L, const double   u_R, const double   s_u_L, const double   s_u_R,
 const double   p_L, const double   p_R, const double   s_p_L, const double   s_p_R,
 const double gamma, const double eps);

#define Riemann_solver_exact Riemann_solver_exact_Ben

#endif
