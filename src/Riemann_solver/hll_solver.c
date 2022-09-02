#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define max(x,y)  ( x>y?x:y )
#define min(x,y)  (x>y?y:x)



void HLL_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L,double n_x, double n_y, double P_R, double RHO_R, double U_R,double V_R,double *lambda_max)
{	
	double H_L, H_R;
	H_L = gamma/(gamma-1.0)*P_L/RHO_L + 0.5*(U_L*U_L+V_L*V_L);
	H_R = gamma/(gamma-1.0)*P_R/RHO_R + 0.5*(U_R*U_R+V_R*V_R);

	double RHO_S, U_S, V_S, H_S, C_S;
	RHO_S = sqrt(RHO_L*RHO_R);
	U_S = (U_L*sqrt(RHO_L)+U_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	V_S = (V_L*sqrt(RHO_L)+V_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	H_S = (H_L*sqrt(RHO_L)+H_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	C_S = sqrt((gamma-1.0)*(H_S-0.5*U_S*U_S-0.5*V_S*V_S));
	double C_L, C_R;
	C_L = sqrt(gamma*P_L/RHO_L);
	C_R = sqrt(gamma*P_R/RHO_R);
	
	double qn_S, qt_S;
	qn_S = U_S*n_x + V_S*n_y;
	qt_S = -U_S*n_y + V_S*n_x;
	double qn_R, qt_R;
	qn_R = U_R*n_x + V_R*n_y;
	qt_R = -U_R*n_y + V_R*n_x;
	double qn_L, qt_L;
	qn_L = U_L*n_x + V_L*n_y;
	qt_L = -U_L*n_y + V_L*n_x;

	double E_L, E_R;
	E_L = 1.0/(gamma-1.0)*P_L/RHO_L + 0.5*(U_L*U_L+V_L*V_L);
	E_R = 1.0/(gamma-1.0)*P_R/RHO_R + 0.5*(U_R*U_R+V_R*V_R);

	double S_L, S_R;
	S_L = min(qn_L-C_L, qn_S-C_S);
	S_R = max(qn_R+C_R, qn_S+C_S);
	
	S_L = min(0,S_L);
	S_R = max(0,S_R);
	
		
	F[0] = (S_R*RHO_L*U_L-S_L*RHO_R*U_R)*n_x + (S_R*RHO_L*V_L-S_L*RHO_R*V_R)*n_y;
	F[0] = F[0]/(S_R-S_L)+S_R*S_L/(S_R-S_L)*(RHO_R - RHO_L);
	F[1] = (S_R*RHO_L*U_L*U_L+S_R*P_L-S_L*RHO_R*U_R*U_R-S_L*P_R)*n_x + (S_R*RHO_L*U_L*V_L-S_L*RHO_R*U_R*V_R)*n_y;
	F[1] = F[1]/(S_R-S_L)+S_R*S_L/(S_R-S_L)*(RHO_R*U_R - RHO_L*U_L);
	F[2] = (S_R*RHO_L*U_L*V_L-S_L*RHO_R*U_R*V_R)*n_x + (S_R*RHO_L*V_L*V_L+S_R*P_L-S_L*RHO_R*V_R*V_R-S_L*P_R)*n_y;
	F[2] = F[2]/(S_R-S_L)+S_R*S_L/(S_R-S_L)*(RHO_R*V_R - RHO_L*V_L);
	F[3] = (S_R*RHO_L*U_L*H_L-S_L*RHO_R*U_R*H_R)*n_x + (S_R*RHO_L*V_L*H_L-S_L*RHO_R*V_R*H_R)*n_y;
	F[3] = F[3]/(S_R-S_L)+S_R*S_L/(S_R-S_L)*(RHO_R*E_R - RHO_L*E_L);
		   

	* lambda_max = fabs(qn_S)+C_S;	  
}


