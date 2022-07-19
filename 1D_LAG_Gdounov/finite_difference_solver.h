#ifndef U_MIN_BURGERS
#define U_MIN_BURGERS 0.0
#endif /* U_MIN_BURGERS */


void Godunov_solver_source
(double * config, int m, double * RHO[], double * U[], double * P[],
 double * E[], double * X[],
 double * MASS, double * RHOL, double * UL, double * PL,
 double * RHOR, double * UR, double * PR, double * cpu_time);
