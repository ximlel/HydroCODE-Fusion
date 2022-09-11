#ifndef _RIEMANN_H
#define _RIEMANN_H

double Riemann_solver_exact(double &P_star, double &U_star, double &rho_starL,double &rho_starR,double rho_L,double rho_R, double u_L, double u_R, double p_L, double p_R, double c_L, double c_R, double gamma_L, double gamma_R)
{
	if(P_star<=p_L)//Left rarefaction wave
		rho_starL=rho_L*pow(P_star/p_L,1./gamma_L);
	else//Left shock wave
		rho_starL=rho_L*(P_star/p_L+(gamma_L-1.)/(gamma_L+1.))/(P_star/p_L*(gamma_L-1.)/(gamma_L+1.)+1.);
	if(P_star<=p_R)//Right rarefaction wave
		rho_starR=rho_R*pow(P_star/p_R,1./gamma_R);
	else//Right shock wave
		rho_starR=rho_R*(P_star/p_R+(gamma_R-1.)/(gamma_R+1.))/(P_star/p_R*(gamma_R-1.)/(gamma_R+1.)+1.);
}
#endif
