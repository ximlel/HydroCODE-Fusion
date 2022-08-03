void Godunov_solver_source
(double * config, int m, double * RHO[], double * U[], double * P[],
 double * E[], double * X[],
 double * MASS, double * RHOL, double * UL, double * PL,
 double * RHOR, double * UR, double * PR, double * cpu_time);

void linear_GRP_solver_LAG
(double * direvative, double * mid,
 double rho_L, double rho_R, double s_rho_L, double s_rho_R,
 double   u_L, double   u_R, double   s_u_L, double   s_u_R,
 double   p_L, double   p_R, double   s_p_L, double   s_p_R,
 double gamma, double eps);

void GRP_solver_source
(double * config, int m, double * RHO[], double * U[], double * P[],
 double * E[], double * X[], double * MASS,
 double * RHOL, double * UL, double * PL,
 double * RHOR, double * UR, double * PR,
 double * SRHOL, double * SUL, double * SPL,
 double * SRHOR, double * SUR, double * SPR,
 double * cpu_time);
