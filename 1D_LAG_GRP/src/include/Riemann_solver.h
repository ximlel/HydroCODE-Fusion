/**
 * @file Riemann_solver.h
 * @brief This file is the header file of several Riemann solvers.
 * @details This header file declares functions in files 'Riemann_solver_exact.c'
 *          and 'Riemann_solver_exact_Toro.c'.
 */

double Riemann_solver_exact(double * U_star, double * P_star, const double gamma,
			    const double u_L, const double u_R, const double p_L, const double p_R,
			    const double c_L, const double c_R, int * CRW,
			    const double eps, const double tol, const int N);

double Riemann_solver_exact_Toro(double * U_star, double * P_star, const double gamma,
				 const double U_l, const double U_r, const double P_l, const double P_r,
				 const double c_l, const double c_r, int * CRW,
				 const double eps, const double tol, const int N);
