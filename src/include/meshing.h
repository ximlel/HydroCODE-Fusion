/**
 * @file meshing.h
 * @brief This file is the header file of unstructured grid processing program.
 * @details This header file declares functions in the folder 'meshing'.
 */

#ifndef MESHING_H
#define MESHING_H

// quad_mesh.c
int quad_mesh(struct mesh_var * mv, const char * mesh_name);

void cylinder_mesh(struct mesh_var * mv);
void odd_even_mesh(struct mesh_var * mv);
void odd_even_periodic_mesh(struct mesh_var * mv);
void odd_even_inflow_mesh(struct mesh_var * mv);
void rand_disturb_inflow_mesh(struct mesh_var * mv);
void Saltzman_Lag_mesh(struct mesh_var * mv);


// msh_load.c
int msh_read(FILE * fp, struct mesh_var * mv);


// mesh_int_free.c
struct mesh_var mesh_init(const char *example, const char *mesh_name);
void mesh_mem_free(struct mesh_var * mv);


// ghost_cell.c
void period_cell_modify(struct mesh_var * mv);
void period_ghost(struct cell_var * cv, const struct mesh_var * mv, struct flu_var * FV, const double t);


// spher_mesh.c
struct spher_mesh_var spher_mesh_init(const char *example);
void spher_mesh_update  (struct spher_mesh_var *smv);
void spher_mesh_mem_free(struct spher_mesh_var *smv);

#endif
