#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define max(x,y)  ( x>y?x:y )
#define min(x,y)  ( x<y?x:y )


void Roe_Goundov_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double P_R, double RHO_R, double U_R, double *lambda_max, double delta)
{
	//	double const  Q_user = 2.0;
	
	
	double C_L, C_R;
	C_L = sqrt(gamma*P_L/RHO_L);
	C_R = sqrt(gamma*P_R/RHO_R);
	double z = 0.5 *  (gamma-1.0) / gamma;
	
	/*
	double Q, P_pvrs, P_max, P_min, RHO_bar, C_bar;
	P_min = min(P_L,P_R);
	P_max = max(P_L,P_R);
	Q = P_max/P_min;
	RHO_bar = 0.5*(RHO_L+RHO_R);
	C_bar = 0.5*(C_L+C_R);
	P_pvrs = 0.5*(P_L+P_R)+0.5*(U_L-U_R)*RHO_bar*C_bar;

	double A_L,A_R,B_L,B_R;
	A_L = 2.0/(gamma+1.0)/RHO_L;
	A_R = 2.0/(gamma+1.0)/RHO_R;
	B_L = (gamma-1)/(gamma+1)*P_L;
	B_R = (gamma-1)/(gamma+1)*P_R;

	double P_star, U_star, U_star_L, U_star_R, RHO_star_L, RHO_star_R, C_star_L, C_star_R, P_0, g_L_0, g_R_0, lambda_L_1, lambda_R_1, lambda_L_3, lambda_R_3;

	if(Q<Q_user&&P_min<P_pvrs&&P_pvrs<P_max) //PVRS
		{
			P_star = max(0,P_pvrs);
			U_star = 0.5*(U_L+U_R)+0.5*(P_L-P_R)/(RHO_bar*C_bar);
			RHO_star_L = RHO_L + (U_L-U_star)*RHO_bar/C_bar;
			RHO_star_R = RHO_R + (U_star - U_R)*RHO_bar/C_bar;
			C_star_L = sqrt(gamma*P_star/RHO_star_L);
			C_star_R = sqrt(gamma*P_star/RHO_star_R);
			U_star_L = U_star;
			U_star_R = U_star;
		}
	else if(P_pvrs<P_min) //TRRS
		{	   	
			P_star = pow((C_L + C_R - (gamma-1.0)/2.0*(U_R - U_L))/(C_L/pow(P_L,z) + C_R/pow(P_R,z)), 1.0/z);
			C_star_L = C_L * pow(P_star/P_L, z);
			U_star_L = U_L + 2.0/(gamma - 1.0)*(C_L - C_star_L);
			C_star_R = C_R * pow(P_star/P_R, z);
			U_star_R = U_R + 2.0/(gamma - 1.0)*(C_star_R - C_R);
		}
	else //TSRS
		{
			P_0 = max(0,P_pvrs);
			g_L_0 = sqrt(A_L/(P_0+B_L));
			g_R_0 = sqrt(A_R/(P_0+B_R));
			P_star = (g_L_0*P_L+g_R_0*P_R-(U_R-U_L))/(g_L_0+g_R_0);
			U_star = 0.5*(U_R+U_L)+0.5*((P_star-P_R)*g_R_0-(P_star-P_L)*g_L_0);
			RHO_star_L = RHO_L*(P_star/P_L+(gamma-1.0)/(gamma+1.0))/((gamma-1.0)*P_star/(gamma+1.0)/P_L+1.0);
			RHO_star_R = RHO_R*(P_star/P_R+(gamma-1.0)/(gamma+1.0))/((gamma-1.0)*P_star/(gamma+1.0)/P_R+1.0);
			C_star_L = sqrt(gamma*P_star/RHO_star_L);
			C_star_R = sqrt(gamma*P_star/RHO_star_R);
			U_star_L = U_star;
			U_star_R = U_star;
		}
	
	lambda_L_1 = U_L - C_L;
	lambda_R_1 = U_star_L - C_star_L;
	lambda_L_3 = U_star_R + C_star_R;
	lambda_R_3 = U_R + C_R;
	*/

	double H_L, H_R;
	H_L = gamma/(gamma-1.0)*P_L/RHO_L + 0.5*(U_L*U_L);
	H_R = gamma/(gamma-1.0)*P_R/RHO_R + 0.5*(U_R*U_R);

	double E_L, E_R;
	E_L = 1.0/(gamma-1.0)*P_L/RHO_L + 0.5*(U_L*U_L);
	E_R = 1.0/(gamma-1.0)*P_R/RHO_R + 0.5*(U_R*U_R);

	double U[3];
	U[0] = RHO_L;
	U[1] = RHO_L*U_L;
	U[2] = RHO_L*E_L;

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
	
		
	lambda[0] = U_S - C_S;
	lambda[1] = U_S;
 	lambda[2] = U_S + C_S;
	

	*lambda_max = fabs(U_S) + C_S;
	
	/*
	if(lambda_L_1<0&&lambda_R_1>0)
		{
			for(i = 0; i < 3; i++)
				{					
					U[i] += (lambda_R_1-(U_S-C_S))/(lambda_R_1-lambda_L_1)*W[0]*R[i][0];				
				}
		}
	else if(lambda_L_3<0&&lambda_R_3>0)
		{
			U[0] = RHO_R;
			U[1] = RHO_R*U_R;
			U[2] = RHO_R*E_R;
			for(i = 0; i < 3; i++)
				{
					U[i] += -((U_S+C_S)-lambda_L_3)/(lambda_R_3-lambda_L_3)*W[2]*R[i][2];				
				}
		}
	else
		{			
	*/		for(j = 0; j < 3; j++)
				{
					if(lambda[j]<=0)
						{						   						
							for(i = 0; i < 3 ; i++)								
								{
									U[i] += W[j]*R[i][j];				
								}
						}
				}
	//		}
	
	F[0] = U[0];
	F[1] = U[1]/U[0];
	F[2] = (U[2]/U[0] - 0.5*F[1]*F[1])*(gamma-1.0)*F[0];
		
}
