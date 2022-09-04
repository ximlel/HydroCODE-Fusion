/**
 * @file meshing.h
 * @brief This file is the header file of unstructured grid processing program.
 * @details This header file declares functions in the folder 'meshing'.
 */

#ifndef MESHING_H
#define MESHING_H


void Sod_mesh(struct mesh_var * mv);
void Shock_Bubble_mesh(struct mesh_var * mv);
void Shear_mesh(struct mesh_var * mv);
void free_mesh(struct mesh_var * mv);
void Shell_mesh(struct mesh_var * mv);
void RMI_mesh(struct mesh_var * mv);
void RMI_S_mesh(struct mesh_var * mv);
void R2D_mesh(struct mesh_var * mv);
void Vortex_mesh(struct mesh_var * mv);
void cylinder_mesh(struct mesh_var * mv);
void odd_even_mesh(struct mesh_var * mv);
void odd_even_periodic_mesh(struct mesh_var * mv);
void odd_even_inflow_mesh(struct mesh_var * mv);
void rand_disturb_inflow_mesh(struct mesh_var * mv);
void oblique_periodic_mesh(struct mesh_var * mv);
void Saltzman_mesh_Lag(struct mesh_var * mv);


int msh_read(FILE * fp, struct mesh_var * mv);


struct mesh_var mesh_load(const char *example, const char *mesh_name);


void period_cell_modi(struct mesh_var * mv);
void period_ghost(struct cell_var * cv, const struct mesh_var * mv, struct flu_var * FV, double t);

#endif
