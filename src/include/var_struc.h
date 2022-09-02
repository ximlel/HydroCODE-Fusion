/**
 * @file var_struc.h
 * @brief This file is the header file of some globally common variables and structural bodies.
 */

#ifndef VARSTRUC_H
#define VARSTRUC_H

#ifdef DOXYGEN_PREDEFINED
/**
 * @def MULTIFLUID_BASICS
 * @brief Switch whether to compute multi-fluids.
 */
#define MULTIFLUID_BASICS
/**
 * @def MULTIPHASE_BASICS
 * @brief Switch whether to compute multi-phase flow. (essential macro: MULTIFLUID_BASICS)
 */
#define MULTIPHASE_BASICS
#endif

//! If the system does not set, the default largest value can be seen as zero is EPS.
#ifndef EPS
#define EPS 1e-9
#endif

//! Define the number of configuration parameters.
#ifndef N_CONF
#define N_CONF 400
#endif

extern double config[]; //!< Initial configuration data array.

//! pointer structure of FLUid VARiables.
typedef struct flu_var {
	double * RHO,   * U,   * V,   * P;
#ifdef MULTIFLUID_BASICS
	double * Z_a;
#ifdef MULTIPHASE_BASICS
	double * RHO_b, * U_b, * V_b, * P_b;
#else
	double * PHI, * gamma;
#endif
#endif
} Fluid_Variable;

//! pointer structure of VARiables on STRUctural computational grid CELLs.
typedef struct cell_var_stru {
	double ** RHO, ** U, ** V, ** P, ** E;   //!< density, velocity components in direction x and y, pressure, specific total energy.
	double  * d_rho,  * d_u,          * d_p; //!< spatial derivatives in one dimension.
	double ** s_rho, ** s_u, ** s_v, ** s_p; //!< spatial derivatives in coordinate x (slopes).
	double ** t_rho, ** t_u, ** t_v, ** t_p; //!< spatial derivatives in coordinate y (slopes).
	double ** rhoIx, ** uIx, ** vIx, ** pIx; //!< interfacial variable values in coordinate x at t_{n+1}.
	double ** rhoIy, ** uIy, ** vIy, ** pIy; //!< interfacial variable values in coordinate y at t_{n+1}.
	double ** F_rho, ** F_e, ** F_u, ** F_v; //!< numerical fluxes at (x_{j-1/2}, t_{n}).
	double ** G_rho, ** G_e, ** G_u, ** G_v; //!< numerical fluxes at (y_{j-1/2}, t_{n}).
} Cell_Variable_Structured;

//! pointer structure of VARiables on unstructural computational grid CELLs.
struct cell_var {
	int **cell_cell;
	double **n_x, **n_y;
	double **F_rho, **F_e, **F_u, **F_v;
	double  *U_rho,  *U_e,  *U_u,  *U_v;
	double  **F_p_x,  **F_p_y, **RHO_p, **U_p, **V_p, **P_p;
	double **dt_U_p, **dt_V_p, **dt_F_p_x, **dt_F_p_y;
	double *X_c, *Y_c;
	double *vol, *c, *dist_p;
	double *gradx_rho, *grady_rho;
	double *gradx_e,   *grady_e;
	double *gradx_u,   *grady_u;
	double *gradx_v,   *grady_v;
	double **RHO_star,    **P_star,    **U_qt_star,    **V_qt_star,    **gamma_star;
	double **RHO_minus_c, **P_minus_c, **U_qt_minus_c, **V_qt_minus_c, **gamma_minus_c;
	double **RHO_add_c,   **P_add_c,   **U_qt_add_c,   **V_qt_add_c,   **gamma_add_c;
	double **u_star, **u_minus_c, **u_add_c;
#ifdef MULTIFLUID_BASICS
	double **Z_a_p;
	double *gradx_z_a, *grady_z_a;
#ifdef MULTIPHASE_BASICS

#else
	double **F_phi, **F_gamma, **F_e_a;
	double  *U_phi,  *U_gamma,  *U_e_a;
	double **PHI_p, **gamma_p;
	double *gradx_phi,   *grady_phi;
	double *gradx_gamma, *grady_gamma;
#endif
#endif
} Cell_Variable;

//! Interfacial Fluid VARiables.
typedef struct i_f_var {
	double n_x, n_y;
	double RHO,     P,     U,     V;     //!< variable values at t_{n}.
	double RHO_int, P_int, U_int, V_int; //!< interfacial variables at t_{n+1}.
	double F_rho, F_e, F_u, F_v;         //!< interfacial fluxes at t_{n+1/2}.
	double d_rho, d_p, d_u, d_v;         //!< normal spatial derivatives.
	double t_rho, t_p, t_u, t_v;         //!< tangential spatial derivatives OR spatial derivatives in Lagrangian coordinate ξ
	double lambda_u, lambda_v;           //!< grid moving velocity components in direction x and y
	double gamma;                        //!< specific heat ratio
#ifdef MULTIFLUID_BASICS
	double PHI, d_phi, t_phi; //!< Mass fraction of fluid a.
	double Z_a, d_z_a, t_z_a; //!< Volume fraction of fluid a.
#endif
} Interface_Fluid_Variable;

//! Fluid VARiables at Boundary.
typedef struct b_f_var {
	double  RHO,  P,  U,  V,  H;  //!< H is the grid cell width.
	double SRHO, SP, SU, SV;      //!< spatial derivatives in coordinate x (slopes).
	double TRHO, TP, TU, TV;      //!< spatial derivatives in coordinate y (slopes).
} Boundary_Fluid_Variable;

//mesh
typedef struct mesh_var {
	int num_pt, num_ghost, *cell_type, **cell_pt;
	int num_border[10], *border_pt, *border_cond, *peri_cell, *normal_v;
	double *X, *Y;
	void (*bc)(struct cell_var * cv, struct mesh_var mv, struct flu_var * FV, double t);
} Mesh_Variable;

#endif
