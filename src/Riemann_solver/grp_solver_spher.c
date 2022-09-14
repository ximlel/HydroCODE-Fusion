#include <stdio.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/riemann_solver.h"


void AcousticSLagTangent(double *dire, const struct i_f_var * ifv, double r, double M)
//GRP solver for tangential case. Lagrangian version(moving mesh) for cylindrical case
{
	double DU  = ifv->d_u;
	double DP  = ifv->d_p;
	double RHO = ifv->RHO;
	double U   = ifv->U;
	double C_star = sqrt(ifv->gamma*ifv->P/RHO);

	dire[2] = -RHO*C_star*C_star*(DU+(M-1)*U/r);
	dire[1] = -DP/RHO;
}

static void AcousticSLag(double *dire, const struct i_f_var * ifv_L, const struct i_f_var * ifv_R, double r, double M)
//GRP solver for acoustic case. Lagrangian version(moving mesh) for cylindrical case
{
	double DUL = ifv_L->d_u;
	double DUR = ifv_R->d_u;
	double DPL = ifv_L->d_p;
	double DPR = ifv_R->d_p;
	double RHO = 0.5*(ifv_L->RHO+ifv_R->RHO);
	double U   = 0.5*(ifv_L->U  +ifv_R->U);
	double C_star = 0.5*(sqrt(ifv_L->gamma*ifv_L->P/ifv_L->RHO)+sqrt(ifv_R->gamma*ifv_R->P/ifv_R->RHO));

	dire[2] =  0.5*C_star*(DPR-DPL)-RHO*C_star*C_star*(0.5*(DUL+DUR)+(M-1)*U/r);
	dire[1] = -0.5*(DPR+DPL)/RHO+0.5*C_star*(DUR-DUL);
	dire[0] = dire[2]/C_star/C_star;
	dire[3] = dire[0];
}

void GRPsolverSLag(double *wave_speed, double *dire, double *U_star, const struct i_f_var * ifv_L, const struct i_f_var * ifv_R, double r, double M)
//GRP solver for non acoustic case. Lagrangian version(moving mesh) cylindrical case
{
	double GammaL = ifv_L->gamma, GammaR = ifv_R->gamma;
	double DL = ifv_L->RHO,    DR = ifv_R->RHO;
	double UL = ifv_L->U,      UR = ifv_R->U;
	double PL = ifv_L->P,      PR = ifv_R->P;
	double DDL = ifv_L->d_rho, DDR = ifv_R->d_rho;
	double DUL = ifv_L->d_u,   DUR = ifv_R->d_u;
	double DPL = ifv_L->d_p,   DPR = ifv_R->d_p;

	StarPU(U_star, ifv_L, ifv_R);
	double UM = U_star[1], PM = U_star[2];
	//P_star, U_star, rho_starL, roh_starR, max wave speed

	const double eps = config[4];
	double aL,bL,dL,aR,bR,dR,CL,CR,C_starL,C_starR,D,U,P,theta,phic,phi1,phi2,phi3,sigmaL,sigmaR,musL,musR;
	musL=(GammaL-1.)/(GammaL+1.);
	musR=(GammaR-1.)/(GammaR+1.);
	CL=sqrt(GammaL*PL/DL);
	CR=sqrt(GammaR*PR/DR);
	double TDSL=-CL*CL/(DL*(GammaL-1.))*DDL+1./(DL*(GammaL-1.))*DPL;
	double TDSR=-CR*CR/(DR*(GammaR-1.))*DDR+1./(DR*(GammaR-1.))*DPR;
	double DpsiL=DUL+GammaL/((GammaL-1.)*CL*DL)*DPL-CL/(DL*(GammaL-1.))*DDL;
	double DphiR=DUR-GammaR/((GammaR-1.)*CR*DR)*DPR+CR/(DR*(GammaR-1.))*DDR;
	if(PM>PL) //left shock
		wave_speed[0]=UL-CL*sqrt(PM/PL*(GammaL+1)/(2.*GammaL)+(GammaL-1.)/(2.*GammaL));
	else //left fan
		wave_speed[0]=UL-CL;
	if(PM>PR) //right shock
		wave_speed[1]=UR+CR*sqrt(PM/PR*(GammaR+1)/(2.*GammaR)+(GammaR-1.)/(2.*GammaR));
	else //right fan
		wave_speed[1]=UR+CR;

	if((DL-DR)*(DL-DR)+(UL-UR)*(UL-UR)+(PL-PR)*(PL-PR)+(CL-CR)*(CL-CR)<eps)//Acoustic case
		{
			AcousticSLag(dire, ifv_L, ifv_R, r, M);
			return;
		}

	// non-Acoustic case
	if(PM<=PL)//Left rarefaction wave
		{	// middle left state
			D=DL*pow(PM/PL,1./GammaL);
			U=UM;
			P=PM;
			C_starL=sqrt(GammaL*P/D);
			theta=C_starL/CL;
			aL=1.;
			bL=1./(D*C_starL);
			if(fabs(GammaL-5./3.)<eps)
				phic=-2.*(3.*C_starL*log(theta)+(UL+2.*CL/(GammaL-1.))*(1.-theta));
			else if(fabs(GammaL-3.)<eps)
				phic=CL-C_starL+(UL+2.*CL/(GammaL-1.))*log(theta);
			else
				phic=(musL-1.)*C_starL/(musL*(4.*musL-1.))*(1.-pow(theta,(1.-4.*musL)/(2.*musL)))+(UL+2.*CL/(GammaL-1.))/(2.*musL-1.)*(1.-pow(theta,(1.-2*musL)/(2.*musL)));
			dL=((1.+musL)/(1.+2.*musL)*pow(theta,0.5/musL)+musL/(1.+2.*musL)*pow(theta,(1.+musL)/musL))*TDSL
				-pow(theta,0.5/musL)*CL*(DpsiL+(M-1.)/(2.*r)*UL)+(M-1.)/(2.*r)*C_starL*(phic-U);
		}
	else//Left shock wave
		{   //middle Left state(behind shock)
			D=DL*(PM/PL+(GammaL-1.)/(GammaL+1.))/(PM/PL*(GammaL-1.)/(GammaL+1.)+1.);
			U=UM;
			P=PM;
			C_starL=sqrt(GammaL*P/D);
			phi1=0.5*sqrt((1.-musL)/(DL*(PM+musL*PL)))*(PM+PL*(1.+2.*musL))/(PM+musL*PL);
			phi2=-0.5*sqrt((1.-musL)/(DL*(PM+musL*PL)))*(PM*(2.+musL)+musL*PL)/(PM+musL*PL);
			phi3=-0.5*(PM-PL)/DL*sqrt((1.-musL)/(DL*(PM+musL*PL)));
			sigmaL=(D*U-DL*UL)/(D-DL);
			aL=1.-D*(sigmaL-U)*phi1;
			bL=phi1-(sigmaL-U)/(D*C_starL*C_starL);
			dL=(-(sigmaL-UL)*phi2-1./DL)*DPL+(sigmaL-UL+DL*phi3+CL*CL*DL*phi2)*DUL-(sigmaL-UL)*phi3*DDL+(DL*UL*CL*CL*phi2+DL*UL*phi3+(sigmaL-U)*U)*(M-1.)/r;
		}

	if(PM<=PR)//Right rarefaction wave
		{	//middle right state
			D=DR*pow(PM/PR,1./GammaR);
			U=UM;
			P=PM;
			C_starR=sqrt(GammaR*P/D);
			aR=1.;
			bR=-1./(D*C_starR);
			theta=C_starR/CR;
			if(fabs(GammaR-5./3.)<eps)
				phic=-2.*(3.*C_starR*log(theta)-(UR-2.*CR/(GammaR-1.))*(1.-theta));
			else if(fabs(GammaR-3.)<eps)
				phic=CR-C_starR-(UR-2.*CR/(GammaR-1.))*log(theta);
			else
				phic=(musR-1.)*C_starR/(musR*(4.*musR-1.))*(1.-pow(theta,(1.-4.*musR)/(2.*musR)))-(UR-2.*CR/(GammaR-1.))/(2.*musR-1.)*(1.-pow(theta,(1.-2*musR)/(2.*musR)));
			dR=((1.+musR)/(1.+2.*musR)*pow(theta,0.5/musR)+musR/(1.+2.*musR)*pow(theta,(1.+musR)/musR))*TDSR
				+pow(theta,0.5/musR)*CR*(DphiR+(M-1.)/(2.*r)*UR)+(M-1.)/(2.*r)*C_starR*(phic+U);
		}
	else//Right shock wave
		{   //middle Right state (behind shock)
			D=DR*(PM/PR+(GammaR-1.)/(GammaR+1.))/(PM/PR*(GammaR-1.)/(GammaR+1.)+1.);
			U=UM;
			P=PM;
			C_starR=sqrt(GammaR*P/D);
			phi1=0.5*sqrt((1.-musR)/(DR*(PM+musR*PR)))*(PM+PR*(1.+2.*musR))/(PM+musR*PR);
			phi2=-0.5*sqrt((1.-musR)/(DR*(PM+musR*PR)))*(PM*(2.+musR)+musR*PR)/(PM+musR*PR);
			phi3=-0.5*(PM-PR)/DR*sqrt((1.-musR)/(DR*(PM+musR*PR)));
			sigmaR=(D*U-DR*UR)/(D-DR);
			aR=1.+D*(sigmaR-UM)*phi1;
			bR=-phi1-(sigmaR-UM)/(D*C_starR*C_starR);
			dR=((sigmaR-UR)*phi2-1./DR)*DPR+(sigmaR-UR-DR*phi3-CR*CR*DR*phi2)*DUR+(sigmaR-UR)*phi3*DDR-(DR*UR*CR*CR*phi2+DR*UR*phi3-(sigmaR-U)*U)*(M-1.)/r;
		}
	dire[1] = (dL*bR-dR*bL)/(aL*bR-aR*bL);
	dire[2] = (dL*aR-dR*aL)/(bL*aR-bR*aL);
	dire[0] = dire[2]/C_starL/C_starL;
	dire[3] = dire[2]/C_starR/C_starR;
}
