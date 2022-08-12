/**
 * @file finite_difference_solver.h
 * @brief This file is the header file of Lagrangian/Eulerian hydrocode in finite volume framework.
 * @details This header file declares functions in files 'Godunov_solver_source.c',
 *          and 'GRP_solver_source.c'.
 */

#ifndef FINITEVOLUME_H
#define FINITEVOLUME_H

void Godunov_solver_LAG_source
(double * config, const int m, double * RHO[], double * U[], double * P[],
 double * E[], double * X[], double * cpu_time);

void GRP_solver_LAG_source
(double * config, const int m, double * RHO[], double * U[], double * P[],
 double * E[], double * X[], double * cpu_time);

void Godunov_solver_EUL_source
(double * config, const int m, double * RHO[], double * U[], double * P[],
 double * E[], double * X[], double * cpu_time);

void GRP_solver_EUL_source
(double * config, const int m, double * RHO[], double * U[], double * P[],
 double * E[], double * X[], double * cpu_time);

#endif
