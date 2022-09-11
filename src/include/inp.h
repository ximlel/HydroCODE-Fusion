#ifndef _INP_H
#define _INP_H


void wrin2s(FILE *out,double **R,double **Z,double **D,double **U,double **P,double **Gamma,double t);//for mesh
void Write(FILE *out,double *XX,int N);// write data in file

void GuessP(double *PStar,double DL,double DR,double UL,double UR,double PL,double PR,double CL,double CR,double GammaL, double GammaR);
//provide a guess value for pressure PM in the star region.AIRS approximate Riemann solvers.
void PreFun(double *F,double *FD,double P,double DK,double PK,double CK,double Gamma);
// evaluate the pressure functions FL and FR in exact Riemann solver
void StarPU(double *POut,double *U,double *DML,double*DMR,double DL,double DR,double UL,double UR,double PL,double PR,double CL,double CR,double GammaL, double GammaR);
// compute the solution for pressure and velocity in the star region
void AcousticSLag(double *DtDL,double *DtDR,double *DtU,double *DtP,double DUL,double DUR,double DPL,double DPR,
					double D,double U,double C_star,double r);
//GRP solver for acoustic case. Lagrangian version(moving mesh) for cylindrical case
void AcousticSLagTangent(double *DtP,double *DtU,double DU,double DP,double D,double U,double C_star,double r);
//GRP solver for tangential case. Lagrangian version(moving mesh) for cylindrical case
void GRPsolverSLag(double *DtDL,double *DtDR,double *DtU,double *DtP,double UM,double PM,
					 double DL,double DR,double UL,double UR,double PL,double PR,
					 double DDL,double DDR,double DUL,double DUR,double DPL,double DPR,
					 double TDSL,double TDSR,double DpsiL,double DphiR,double r,double GammaL,double GammaR);
//GRP solver for non acoustic case. Lagrangian version(moving mesh) cylindrical case

#endif
