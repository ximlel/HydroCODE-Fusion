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
	double *RHO, *U, *V, *P;
} S_FV;

//! Pointer structural body of variables on structural computational grid cells.
typedef struct cell_var_stru {
	double **RHO, **U, **V, **P, **E;
} S_CVS;

#endif
