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

//! Pointer structural body of fluid variables.
typedef struct flu_var {
	double * RHO, * U, * V, * P;
} Fluid_Variable;

//! Pointer structural body of variables on structural computational grid cells.
typedef struct cell_var_stru {
	double ** RHO, ** U, ** V, ** P, ** E;
	double ** s_rho, ** s_u, ** s_v, ** s_p; // spatial derivatives in coordinate x (slopes).
	double ** t_rho, ** t_u, ** t_v, ** t_p; // spatial derivatives in coordinate y (slopes).
	double ** rhoIx, ** uIx, ** vIx, ** pIx; // interfacial variable values in coordinate x.
	double ** rhoIy, ** uIy, ** vIy, ** pIy; // interfacial variable values in coordinate y.
	double ** F_rho, ** F_e, ** F_u, ** F_v; // numerical fluxes at (x_{j-1/2}, t_{n}).
	double ** G_rho, ** G_e, ** G_u, ** G_v; // numerical fluxes at (y_{j-1/2}, t_{n}).
} Cell_Variable_Structured;

//! Interfacial fluid variables.
typedef struct i_f_var {
	double n_x, n_y;
	double RHO,     P,     U,     V;
	double RHO_int, P_int, U_int, V_int;
	double F_rho, F_e, F_u, F_v; // interfacial flux
	double d_rho, d_p, d_u, d_v; // normal spatial derivatives
	double t_rho, t_p, t_u, t_v; // tangential spatial derivatives
} Interface_Fluid_Variable;

//! Fluid variables at boundary.
typedef struct b_f_var {
	double  RHO,  P,  U,  V;
	double SRHO, SP, SU, SV;
	double TRHO, TP, TU, TV;
} Boundary_Fluid_Variable;

#endif
