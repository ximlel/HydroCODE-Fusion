/**
 * @file finite_volume.h
 * @brief This file is the header file of Lagrangian/Eulerian hydrocode in finite volume framework.
 * @details This header file declares functions in the folder 'finite_volume'.
 */

#ifndef FINITEVOLUME_H
#define FINITEVOLUME_H

#include "../include/var_struc.h"

/* 1-D Godunov/GRP scheme (Lagrangian, single-component flow) */
//////////////////////////////////////
// godunov_solver_LAG_source.c
//////////////////////////////////////
void Godunov_solver_LAG_source(const int m, struct cell_var_stru CV, double * X[], double * cpu_time, int * N_plot, double time_plot[]);
//////////////////////////////////////
// grp_solver_LAG_source.c
//////////////////////////////////////
void     GRP_solver_LAG_source(const int m, struct cell_var_stru CV, double * X[], double * cpu_time, int * N_plot, double time_plot[]);

/* radially symmertric Godunov/GRP scheme (Lagrangian, two-component flow, radial structured grid) */
//////////////////////////////////////
// grp_solver_radial_LAG_source.c
//////////////////////////////////////
void GRP_solver_radial_LAG_source(struct cell_var_stru CV, struct radial_mesh_var * rmv, double * R[], const int M,
				  double * cpu_time, const char * problem, int N_T, int * N_plot , double time_plot[]);

/* 1-D Godunov/GRP scheme (Eulerian, single-component flow) */
//////////////////////////////////////
// godunov_solver_EUL_source.c
//////////////////////////////////////
void Godunov_solver_EUL_source(const int m, struct cell_var_stru CV, double * cpu_time, int * N_plot, double time_plot[]);
//////////////////////////////////////
// grp_solver_EUL_source.c
//////////////////////////////////////
void     GRP_solver_EUL_source(const int m, struct cell_var_stru CV, double * cpu_time, int * N_plot, double time_plot[]);

/* 2-D Godunov/GRP scheme (Eulerian, single-component flow, structured grid) */
//////////////////////////////////////
// grp_solver_2D_EUL_source.c
//////////////////////////////////////
void GRP_solver_2D_EUL_source(const int m, const int n, struct cell_var_stru * CV, double ** X, double **Y, 
			      double * cpu_time, const char * problem, int N_T, int * N_plot, double time_plot[]);
//////////////////////////////////////
// grp_solver_2D_split_EUL_source.c
//////////////////////////////////////
void GRP_solver_2D_split_EUL_source(const int m, const int n, struct cell_var_stru * CV, double ** X, double **Y, 
                                    double * cpu_time, const char * problem, int N_T, int * N_plot, double time_plot[]);

/* 2-D Godunov/GRP scheme (Eulerian, two-component flow, unstructured grid) */
//////////////////////////////////////
// finite_volume_scheme_unstruct.c
//////////////////////////////////////
void finite_volume_scheme_unstruct(struct flu_var * FV, const struct mesh_var * mv, const char * scheme, 
				   const char * problem, int * N_plot, const double time_plot[]);

/* 2-D Godunov/GRP scheme (Eulerian, Baer-Nunziato two-phase flow, structured grid) */
void finite_volume_scheme_GRP2D(struct flu_var * FV, const struct mesh_var * mv, const char * phase, const char * problem);

#endif
