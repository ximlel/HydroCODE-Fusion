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
 * @def LAGRANGIAN_MAIRE
 * @brief Switch whether to use Maire's cell-centered schemes for Lagrangian hydrodynamics.
 */
#define LAGRANGIAN_MAIRE
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

//! pointer structure of FLUid VARiables array.
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
	double **    E;                      //!< specific total energy.
	double **  RHO, **  U, **  V, **  P; //!< density, velocity components in direction x and y, pressure.
	double * d_rho, * d_u,        * d_p; //!< spatial derivatives in one dimension.
	double **s_rho, **s_u, **s_v, **s_p; //!< spatial derivatives in coordinate x (slopes).
	double **t_rho, **t_u, **t_v, **t_p; //!< spatial derivatives in coordinate y (slopes).
	double **rhoIx, **uIx, **vIx, **pIx; //!< interfacial variable values in coordinate x at t_{n+1}.
	double **rhoIy, **uIy, **vIy, **pIy; //!< interfacial variable values in coordinate y at t_{n+1}.
	double **F_rho, **F_e, **F_u, **F_v; //!< numerical fluxes at (x_{j-1/2}, t_{n}).
	double **G_rho, **G_e, **G_u, **G_v; //!< numerical fluxes at (y_{j-1/2}, t_{n}).
} Cell_Variable_Structured;

//! pointer structure of VARiables on unstructured computational grid CELLs.
typedef struct cell_var {
	/**
	 * @brief Topological relationships between grid cells.
	 *  @arg  cell_cell[i][0], number of cells adjacent to the i-th cell.
	 *  @arg  cell_cell[i][j](j>0), serial number of j-th adjacent cell to the i-th cell.
	 */
	int    **cell_cell;  
	double * vol;        //!< area(volume) of each grid cell.
	double **n_x, **n_y; //!< x- and y-coordinates of the interfacial unit normal vector.
	double * X_c, * Y_c; //!< x- and y-coordinates of the center point of grid cells.
	double **   F_rho, **   F_e, **   F_u, **   F_v; //!< interfacial fluxes.
	double *    U_rho, *    U_e, *    U_u, *    U_v; //!< conservative variables.
	double **   RHO_p, **   U_p, **   V_p, **   P_p;
	double *gradx_rho, *gradx_e, *gradx_u, *gradx_v; //!< spatial derivatives in coordinate x (gradients).
	double *grady_rho, *grady_e, *grady_u, *grady_v; //!< spatial derivatives in coordinate y (gradients).
#ifdef MULTIFLUID_BASICS
	double **F_e_a,   *U_e_a,   **Z_a_p,   *gradx_z_a,   *grady_z_a;   //!< Total energy OR volume fraction of fluid a.
	double **F_phi,   *U_phi,   **PHI_p,   *gradx_phi,   *grady_phi;   //!< Mass fraction of fluid a.
	double **F_gamma, *U_gamma, **gamma_p, *gradx_gamma, *grady_gamma; //!< Specific heat ratio.
	double **P_star;
	double **U_qt_star,  **V_qt_star;
	double **U_qt_add_c, **V_qt_add_c;
#endif
#ifdef MULTIPHASE_BASICS_abandoned
	double **RHO_star,    **gamma_star,    **P_star;
	double **RHO_minus_c, **gamma_minus_c, **P_minus_c, **U_qt_minus_c, **V_qt_minus_c;
	double **RHO_add_c,   **gamma_add_c,   **P_add_c;
	double **u_star, **u_minus_c, **u_add_c;
#endif
#ifdef LAGRANGIAN_MAIRE
	double **F_p_x,  **F_p_y;
	double **dt_U_p, **dt_V_p, **dt_F_p_x, **dt_F_p_y;
	double  *c, *dist_p;
#endif
} Cell_Variable;

//! Interfacial Fluid VARiables.
typedef struct i_f_var {
	double     n_x,     n_y;             //!< x- and y-coordinates of the interfacial unit normal vector.
	double delta_x, delta_y;             //!< distance from the interfacial center to the grid cell center in direction x and y.
	double length;                       //!< length of the interface.
	double lambda_u, lambda_v;           //!< grid moving velocity components in direction x and y.
	double gamma;                        //!< specific heat ratio.
	double RHO,     P,     U,     V;     //!< primitive variable values at t_{n}.
	double RHO_int, P_int, U_int, V_int; //!< interfacial primitive variables at t_{n+1}.
	double F_rho, F_e, F_u, F_v;         //!< interfacial fluxes at t_{n+1/2}.
	double U_rho, U_e, U_u, U_v;         //!< interfacial conservative variables at t_{n}.
	double d_rho, d_e, d_u, d_v, d_p;    //!< normal spatial derivatives.
	double t_rho, t_e, t_u, t_v, t_p;    //!< tangential spatial derivatives OR spatial derivatives in Lagrangian coordinate ξ.
#ifdef MULTIFLUID_BASICS
	double Z_a, d_z_a,   t_z_a,   F_e_a,   U_e_a;   //!< Volume fraction OR total energy of fluid a.
	double PHI, d_phi,   t_phi,   F_phi,   U_phi;   //!< Mass fraction of fluid a.
	double      d_gamma, t_gamma, F_gamma, U_gamma; //!< Specific heat ratio.
	double P_star;
	double U_qt_star,  V_qt_star;
	double U_qt_add_c, V_qt_add_c;
#endif
#ifdef MULTIPHASE_BASICS_abandoned
	double RHO_star,    gamma_star;
	double RHO_minus_c, gamma_minus_c, P_minus_c, U_qt_minus_c, V_qt_minus_c;
	double RHO_add_c,   gamma_add_c,   P_add_c;
	double u_star, u_minus_c, u_add_c;
#endif
} Interface_Fluid_Variable;

//! Fluid VARiables at Boundary in one direction.
typedef struct b_f_var {
	double    H;              //!< cell width of the ghost grid at boundary.
	double  RHO,  P,  U,  V;  //!< primitive variables on the ghost grid cell at boundary.
	double SRHO, SP, SU, SV;  //!< spatial derivatives in coordinate x (slopes).
	double TRHO, TP, TU, TV;  //!< spatial derivatives in coordinate y (slopes).
} Boundary_Fluid_Variable;

//! MESHing VARiables.
typedef struct mesh_var {
	int num_pt;      //!< Total number of grid nodes.
	int num_ghost;   //!< Total number of ghost grid cells.
	int **cell_pt;   //!< Serial number of each grid node on a grid cell in the clockwise direction.
	int * cell_type; //!< @todo Grid cell type read from mesh file '*.msh' 
	/**
	 * @brief num_border[0] is the total number of connected boundaries on the entire computational domain.
	 * @details num_border[i](i = 1,2,…,9) is the total number of grid cell interfaces at the i-th connected boundary.
	 */
	int num_border[10]; 
	int *border_pt;   //!< Serial number of each grid node on the connected boundary in the clockwise direction.
	/**
	 * @brief Boudary condition indicator of i-th grid cell interface at the boundary.
	 *  @arg  border_cond[i] < 0, refer to 'config[17]' in'doc/config.csv'.
	 *  @arg  border_cond[i] > 0, serial number of the ghost grid cell at the periodic boundary is i.
	 */
	int *border_cond;
	int *period_cell;   //!< Serial number of ghost grid cells at the periodic boundary.
	double *normal_v; //!< @todo Normal velocity on grid cell interfaces at the boundary.
	double *X, *Y;    //!< x- and y-coordinates of the grid nodes with fixed serial number.
	//! Pointer to the boundary condition function.
	void (*bc)(struct cell_var * cv, const struct mesh_var * mv, struct flu_var * FV, double t);
} Mesh_Variable;


//! SPHERically symmetric MESHing VARiables.
typedef struct spher_mesh_var {
    double * Rb;
    double * Lb;
    double * RR;
    double * DdrL;
	double * DdrR;
	double * Ddr;
	double * dRc;
    double * vol;
} Spherical_Mesh_Variable;

#endif
