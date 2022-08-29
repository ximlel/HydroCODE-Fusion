#ifndef VARSTRUC_H
#define VARSTRUC_H

//! If the system does not set, the default largest value can be seen as zero is EPS.
#ifndef EPS
#define EPS 1e-9
#endif

//! Define the number of configuration parameters.
#ifndef N_CONF
#define N_CONF 400
#endif

extern double config[]; //!< Initial configuration data array.

//! pointer structural body of FLUid VARiables.
typedef struct flu_var {
	double * RHO, * U, * V, * P;
} Fluid_Variable;

//! pointer structural body of VARiables on STRUctural computational grid CELLs.
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

#endif
