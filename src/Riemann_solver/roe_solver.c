#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/var_struc.h"


void Roe_solver(double * F, double * lambda_max, const struct i_f_var ifv_L, const struct i_f_var ifv_R, const double delta)
{
	const double gamma = ifv_L.gamma;
	const double P_L   = ifv_L.P,   P_R = ifv_R.P;
	const double RHO_L = ifv_L.RHO, RHO_R = ifv_R.RHO;
	const double U_L   = ifv_L.U,   U_R = ifv_R.U;

	double C_L, C_R;
	C_L = sqrt(gamma*P_L/RHO_L);
	C_R = sqrt(gamma*P_R/RHO_R);
	double z = 0.5 *  (gamma-1.0) / gamma;
	double P_star, U_star_L, U_star_R, C_star_L, C_star_R, lambda_L_1, lambda_R_1, lambda_L_3, lambda_R_3;
	P_star = pow((C_L + C_R - (gamma-1.0)/2.0*(U_R - U_L))/(C_L/pow(P_L,z) + C_R/pow(P_R,z)), 1.0/z);
	C_star_L = C_L * pow(P_star/P_L, z);
	U_star_L = U_L + 2.0/(gamma - 1.0)*(C_L - C_star_L);
	C_star_R = C_R * pow(P_star/P_R, z);
	U_star_R = U_R + 2.0/(gamma - 1.0)*(C_star_R - C_R);
	lambda_L_1 = U_L - C_L;
	lambda_R_1 = U_star_L - C_star_L;
	lambda_L_3 = U_star_R + C_star_R;
	lambda_R_3 = U_R + C_R;
	

	double H_L, H_R;
	H_L = gamma/(gamma-1.0)*P_L/RHO_L + 0.5*(U_L*U_L);
	H_R = gamma/(gamma-1.0)*P_R/RHO_R + 0.5*(U_R*U_R);

	F[0] = 0.5*(RHO_L*U_L+RHO_R*U_R);
	F[1] = 0.5*(RHO_L*U_L*U_L+P_L+RHO_R*U_R*U_R+P_R);
	F[2] = 0.5*(RHO_L*U_L*H_L+RHO_R*U_R*H_R);

	double RHO_S, U_S, H_S, C_S;
	RHO_S = sqrt(RHO_L*RHO_R);
	U_S = (U_L*sqrt(RHO_L)+U_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	H_S = (H_L*sqrt(RHO_L)+H_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	C_S = sqrt((gamma-1.0)*(H_S-0.5*U_S*U_S));
	
	double R[3][3];
	double lambda[3], W[3];
	R[0][0] = 1.0;
	R[0][1] = 1.0;
	R[0][2] = 1.0;
	R[1][0] = U_S - C_S;
	R[1][1] = U_S;
	R[1][2] = U_S + C_S;
	R[2][0] = H_S - U_S*C_S;
	R[2][1] = 0.5*(U_S*U_S);
	R[2][2] = H_S + U_S*C_S;

	int i, j;
	
	W[0] = 0.5*((P_R-P_L)-RHO_S*C_S*(U_R-U_L))/(C_S*C_S);
	W[1] = (RHO_R-RHO_L)-(P_R-P_L)/(C_S*C_S);
	W[2] = 0.5*((P_R-P_L)+RHO_S*C_S*(U_R-U_L))/(C_S*C_S);
	
		
	lambda[0] = fabs(U_S - C_S);
	lambda[1] = fabs(U_S);
	lambda[2] = fabs(U_S + C_S);
	
	*lambda_max = fabs(U_S) + C_S;
	
	if(lambda_L_1<0&&lambda_R_1>0)
		{
			F[0] = RHO_L*U_L;
			F[1] = RHO_L*U_L*U_L+P_L;
			F[2] = RHO_L*U_L*H_L;
			lambda[0] = lambda_L_1*(lambda_R_1-(U_S-C_S))/(lambda_R_1-lambda_L_1);
			for(i = 0; i < 3; i++)
				{					
					F[i] += lambda[0]*W[0]*R[i][0];				
				}
		}
	else if(lambda_L_3<0&&lambda_R_3>0)
		{
			F[0] = RHO_R*U_R;
			F[1] = RHO_R*U_R*U_R+P_R;
			F[2] = RHO_R*U_R*H_R;
			lambda[2] = lambda_R_3*((U_S+C_S)-lambda_L_3)/(lambda_R_3-lambda_L_3);
			for(i = 0; i < 3; i++)
				{
					F[i] += -lambda[2]*W[2]*R[i][2];				
				}
		}
	else
		{			
			for(i = 0; i < 3; i++)
				{
					for(j = 0; j < 3 ; j++)
						{
							F[i] += -0.5*lambda[j]*W[j]*R[i][j];				
						}
				}
		}
}
