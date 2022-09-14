#include <stdio.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/riemann_solver.h"


//GRP solver for tangential case. Lagrangian version(moving mesh) for cylindrical case.
void AcousticSLagTangent(double *dire, double *U_star, const struct i_f_var * ifv_L, const struct i_f_var * ifv_R,
			 double r, double M, const double eps)
{
	double GammaL = ifv_L->gamma, GammaR = ifv_R->gamma;
	double DL = ifv_L->RHO, DR = ifv_R->RHO;
	double UL = ifv_L->U,   UR = ifv_R->U;
	double PL = ifv_L->P,   PR = ifv_R->P;
	double DU = ifv_L->d_u;
	double DP = ifv_L->d_p;

	double CL,CR;
	CL = sqrt(GammaL*PL/DL);
	CR = sqrt(GammaR*PR/DR);

	_Bool CRW[2];
	double UM, PM, DM, DML, DMR;
	//P_star, U_star, rho_star, rho_starL, rho_starR
	Riemann_solver_starPU(&UM, &PM, GammaL, GammaR, UL, UR, PL, PR, CL, CR, CRW, eps, eps, 100);
	if(PM<=PL)//Left rarefaction wave
	    DML=DL*pow(PM/PL,1./GammaL);
	else//Left shock wave
	    DML=DL*(PM/PL+(GammaL-1.)/(GammaL+1.))/(PM/PL*(GammaL-1.)/(GammaL+1.)+1.);
	if(PM<=PR)//Right rarefaction wave
	    DMR=DR*pow(PM/PR,1./GammaR);
	else//Right shock wave
	    DMR=DR*(PM/PR+(GammaR-1.)/(GammaR+1.))/(PM/PR*(GammaR-1.)/(GammaR+1.)+1.);
	DM = 0.5*(DML+DMR);
	double C_star = 0.5*(sqrt(GammaL*PM/DML)+sqrt(GammaR*PM/DMR));

	U_star[1] = UM;
	U_star[2] = PM;
	dire[1] = -DP/DM;
	dire[2] = -DM*C_star*C_star*(DU+(M-1)*UM/r);
}


// Lagrangian version(moving mesh) cylindrical case.
void GRPsolverSLag(double *wave_speed, double *dire, double *U_star, const struct i_f_var * ifv_L, const struct i_f_var * ifv_R,
		   double r, double M, const double eps, const double  atc)
{
	const double GammaL = ifv_L->gamma, GammaR = ifv_R->gamma;
	const double DL = ifv_L->RHO,    DR = ifv_R->RHO;
	const double UL = ifv_L->U,      UR = ifv_R->U;
	const double PL = ifv_L->P,      PR = ifv_R->P;
	const double DDL = ifv_L->d_rho, DDR = ifv_R->d_rho;
	const double DUL = ifv_L->d_u,   DUR = ifv_R->d_u;
	const double DPL = ifv_L->d_p,   DPR = ifv_R->d_p;

	_Bool CRW[2];
	double dist;
	double CL, CR, C_starL, C_starR, C_star, theta;
	CL=sqrt(GammaL*PL/DL);
	CR=sqrt(GammaR*PR/DR);

	double UM, PM, DML, DMR, DM;
	//P_star, U_star, rho_starL, roh_starR, max wave speed

	double aL,bL,dL,aR,bR,dR;
	double phic, phi1, phi2, phi3;
	double sigmaL, sigmaR;
	
	const double musL=(GammaL-1.)/(GammaL+1.);
	const double musR=(GammaR-1.)/(GammaR+1.);
	const double TDSL=-CL*CL/(DL*(GammaL-1.))*DDL+1./(DL*(GammaL-1.))*DPL;
	const double TDSR=-CR*CR/(DR*(GammaR-1.))*DDR+1./(DR*(GammaR-1.))*DPR;
	const double DpsiL=DUL+GammaL/((GammaL-1.)*CL*DL)*DPL-CL/(DL*(GammaL-1.))*DDL;
	const double DphiR=DUR-GammaR/((GammaR-1.)*CR*DR)*DPR+CR/(DR*(GammaR-1.))*DDR;

	Riemann_solver_starPU(&UM, &PM, GammaL, GammaR, UL, UR, PL, PR, CL, CR, CRW, eps, eps, 100);
	if(PM<=PL)//left fan
	    {
		DML=DL*pow(PM/PL,1./GammaL);
		wave_speed[0]=UL-CL;
	    }
	else//left shock
	    {
		DML=DL*(PM/PL+(GammaL-1.)/(GammaL+1.))/(PM/PL*(GammaL-1.)/(GammaL+1.)+1.);
		wave_speed[0]=UL-CL*sqrt(PM/PL*(GammaL+1)/(2.*GammaL)+(GammaL-1.)/(2.*GammaL));
	    }
	C_starL=sqrt(GammaL*PM/DML);
	if(PM<=PR)//right fan
	    {
		DMR=DR*pow(PM/PR,1./GammaR);
		wave_speed[1]=UR+CR;
	    }
	else//right shock
	    {
		DMR=DR*(PM/PR+(GammaR-1.)/(GammaR+1.))/(PM/PR*(GammaR-1.)/(GammaR+1.)+1.);
		wave_speed[1]=UR+CR*sqrt(PM/PR*(GammaR+1)/(2.*GammaR)+(GammaR-1.)/(2.*GammaR));
	    }
	C_starR=sqrt(GammaR*PM/DMR);
	U_star[1] = UM;
	U_star[2] = PM;
	U_star[0] = DML;
	U_star[3] = DMR;

	dist = (DL-DR)*(DL-DR)+(UL-UR)*(UL-UR)+(PL-PR)*(PL-PR)+(CL-CR)*(CL-CR);
	if(dist < atc) //GRP solver for acoustic case.
	    {
		DM = 0.5*(DML + DMR);
		C_star = 0.5*(C_starL+C_starR);
		dire[2] =  0.5*C_star*(DPR-DPL)-DM*C_star*C_star*(0.5*(DUL+DUR)+(M-1)*UM/r);
		dire[1] = -0.5*(DPR+DPL)/DM+0.5*C_star*(DUR-DUL);
		dire[0] = dire[2]/C_star/C_star;
		dire[3] = dire[0];
		return;
	    }

	// GRP solver for non acoustic case.
	if(PM<=PL)//Left rarefaction wave
		{ // middle left state
			theta=C_starL/CL;
			aL=1.;
			bL=1./(DML*C_starL);
			if(fabs(GammaL-5./3.)<eps)
				phic=-2.*(3.*C_starL*log(theta)+(UL+2.*CL/(GammaL-1.))*(1.-theta));
			else if(fabs(GammaL-3.)<eps)
				phic=CL-C_starL+(UL+2.*CL/(GammaL-1.))*log(theta);
			else
				phic=(musL-1.)*C_starL/(musL*(4.*musL-1.))*(1.-pow(theta,(1.-4.*musL)/(2.*musL)))+(UL+2.*CL/(GammaL-1.))/(2.*musL-1.)*(1.-pow(theta,(1.-2*musL)/(2.*musL)));
			dL=((1.+musL)/(1.+2.*musL)*pow(theta,0.5/musL)+musL/(1.+2.*musL)*pow(theta,(1.+musL)/musL))*TDSL
				-pow(theta,0.5/musL)*CL*(DpsiL+(M-1.)/(2.*r)*UL)+(M-1.)/(2.*r)*C_starL*(phic-UM);
		}
	else//Left shock wave
		{ //middle Left state(behind shock)
			phi1=0.5*sqrt((1.-musL)/(DL*(PM+musL*PL)))*(PM+PL*(1.+2.*musL))/(PM+musL*PL);
			phi2=-0.5*sqrt((1.-musL)/(DL*(PM+musL*PL)))*(PM*(2.+musL)+musL*PL)/(PM+musL*PL);
			phi3=-0.5*(PM-PL)/DL*sqrt((1.-musL)/(DL*(PM+musL*PL)));
			sigmaL=(DML*UM-DL*UL)/(DML-DL);
			aL=1.-DML*(sigmaL-UM)*phi1;
			bL=phi1-(sigmaL-UM)/(DML*C_starL*C_starL);
			dL=(-(sigmaL-UL)*phi2-1./DL)*DPL+(sigmaL-UL+DL*phi3+CL*CL*DL*phi2)*DUL-(sigmaL-UL)*phi3*DDL+(DL*UL*CL*CL*phi2+DL*UL*phi3+(sigmaL-UM)*UM)*(M-1.)/r;
		}
	if(PM<=PR)//Right rarefaction wave
		{ //middle right state
			aR=1.;
			bR=-1./(DMR*C_starR);
			theta=C_starR/CR;
			if(fabs(GammaR-5./3.)<eps)
				phic=-2.*(3.*C_starR*log(theta)-(UR-2.*CR/(GammaR-1.))*(1.-theta));
			else if(fabs(GammaR-3.)<eps)
				phic=CR-C_starR-(UR-2.*CR/(GammaR-1.))*log(theta);
			else
				phic=(musR-1.)*C_starR/(musR*(4.*musR-1.))*(1.-pow(theta,(1.-4.*musR)/(2.*musR)))-(UR-2.*CR/(GammaR-1.))/(2.*musR-1.)*(1.-pow(theta,(1.-2*musR)/(2.*musR)));
			dR=((1.+musR)/(1.+2.*musR)*pow(theta,0.5/musR)+musR/(1.+2.*musR)*pow(theta,(1.+musR)/musR))*TDSR
				+pow(theta,0.5/musR)*CR*(DphiR+(M-1.)/(2.*r)*UR)+(M-1.)/(2.*r)*C_starR*(phic+UM);
		}
	else//Right shock wave
		{ //middle Right state (behind shock)
			phi1=0.5*sqrt((1.-musR)/(DR*(PM+musR*PR)))*(PM+PR*(1.+2.*musR))/(PM+musR*PR);
			phi2=-0.5*sqrt((1.-musR)/(DR*(PM+musR*PR)))*(PM*(2.+musR)+musR*PR)/(PM+musR*PR);
			phi3=-0.5*(PM-PR)/DR*sqrt((1.-musR)/(DR*(PM+musR*PR)));
			sigmaR=(DMR*UM-DR*UR)/(DMR-DR);
			aR=1.+DMR*(sigmaR-UM)*phi1;
			bR=-phi1-(sigmaR-UM)/(DMR*C_starR*C_starR);
			dR=((sigmaR-UR)*phi2-1./DR)*DPR+(sigmaR-UR-DR*phi3-CR*CR*DR*phi2)*DUR+(sigmaR-UR)*phi3*DDR-(DR*UR*CR*CR*phi2+DR*UR*phi3-(sigmaR-UM)*UM)*(M-1.)/r;
		}
	dire[1] = (dL*bR-dR*bL)/(aL*bR-aR*bL);
	dire[2] = (dL*aR-dR*aL)/(bL*aR-bR*aL);
	dire[0] = dire[2]/C_starL/C_starL;
	dire[3] = dire[2]/C_starR/C_starR;
}
