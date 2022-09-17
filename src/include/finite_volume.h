/**
 * @file finite_volume.h
 * @brief This file is the header file of Lagrangian/Eulerian hydrocode in finite volume framework.
 * @details This header file declares functions in the folder 'finite_volume'.
 */

#ifndef FINITEVOLUME_H
#define FINITEVOLUME_H

#include "../include/var_struc.h"

// 1-D Godunov/GRP scheme (Lagrangian, single-component flow)
void Godunov_solver_LAG_source(const int m, struct cell_var_stru CV, double * X[], double * cpu_time, const int N_plot, double time_plot[]);
void     GRP_solver_LAG_source(const int m, struct cell_var_stru CV, double * X[], double * cpu_time, const int N_plot, double time_plot[]);

// 1-D Godunov/GRP scheme (Eulerian, single-component flow)
void Godunov_solver_EUL_source(const int m, struct cell_var_stru CV, double * cpu_time, const int N_plot, double time_plot[]);
void     GRP_solver_EUL_source(const int m, struct cell_var_stru CV, double * cpu_time, const int N_plot, double time_plot[]);

// 2-D Godunov/GRP scheme (Eulerian, single-component flow, structured grid)
void GRP_solver_2D_EUL_source      (const int m, const int n, struct cell_var_stru * CV, double * cpu_time, const int N_plot, double time_plot[]);
void GRP_solver_2D_split_EUL_source(const int m, const int n, struct cell_var_stru * CV, double * cpu_time, const int N_plot, double time_plot[]);

// 2-D Godunov/GRP scheme (Eulerian, two-component flow, unstructured grid)
void finite_volume_scheme_unstruct(struct flu_var * FV, const struct mesh_var * mv, const char * scheme, const char * problem,
				   const int N_plot , const double time_plot[]);

// 2-D Godunov/GRP scheme (Eulerian, Baer-Nunziato two-phase flow, structured grid)
void finite_volume_scheme_GRP2D(struct flu_var * FV, const struct mesh_var * mv, const char * phase, const char * problem);

void grp_solver_spher_LAG_source(struct flu_var *FV, struct spher_mesh_var *smv, const int M, double * cpu_time, const int N_plot , double time_plot[]);

#endif
