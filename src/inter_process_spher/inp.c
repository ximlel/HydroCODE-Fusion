#include <stdio.h>
#include <math.h>

#include "../include/var_struc.h"


#define Tcell_plot 100 // Output zoom
#define Ncell 3000 // Number of computing cells in r direction
#define M (2.)    // m=1 planar; m=2 cylindrical; m=3 spherical


void wrin2s(FILE *out,double **R,double **Z,double **D,double **U,double **P,double **Gamma,double t)//for mesh
{
	fprintf(out,"Solution For Euler Equation\n");
	fprintf(out,"variables=Z,R,D,U,P,Gamma\n");
	fprintf(out,"zone I=%d,J=%d,F=POINT,SOLUTIONTIME=%lf\n",Tcell_plot+1,Ncell+1,t);
	int it,jr;
	for(jr=0;jr<=Ncell;jr++)
			for(it=0;it<=Tcell_plot;it++)
				fprintf(out,"%lf %lf %lf %lf %lf %lf\n",Z[jr][it],R[jr][it],D[jr][it],U[jr][it],P[jr][it],Gamma[jr][it]);
}
void Write(FILE *out,double *XX,int N)// write data in file
{
	int i;
	for(i=0;i<N;i++)
		{
			fprintf(out,"%lf ",XX[i]);
		}
}

void GuessP(double *PStar,double DL,double DR,double UL,double UR,double PL,double PR,double CL,double CR,double GammaL, double GammaR)
//provide a guess value for pressure PM in the star region.AIRS approximate Riemann solvers.
{
	const double eps = config[4];
	double PM;
	double CUP,GEL,GER,PMAX,PMIN,PPV,QMAX;
	const double QUSER=2.0;

	CUP=0.25*(DL+DR)*(CL+CR);
	PPV=0.5*(PL+PR)+0.5*(UL-UR)*CUP;
	PPV=(eps>PPV?eps:PPV);
	PMIN=(PL<PR?PL:PR);
	PMAX=(PL>PR?PL:PR);
	QMAX=PMAX/PMIN;
	if(QMAX<=QUSER&&PMIN<=PPV&&PPV<=PMAX)//PVRS solution
		PM=PPV;
	/*
			 else if(PPV<PMIN)//Two rarefaction
			 {
			 //PQ=pow(PL/PR,(Gamma-1.)/(2.*Gamma));
			 //UM=(PQ*UL/CL+UR/CR+2./(Gamma-1.)*(PQ-1.))/(PQ/CL+1./CR);
			 //PTL=1.+(Gamma-1.)/2.*(UL-UM)/CL;
			 //PTR=1.+(Gamma-1.)/2.*(UM-UR)/CR;
			 //PM=0.5*(pow(PL*PTL,2.*Gamma/(Gamma-1.))+pow(PR*PTR,2.*Gamma/(Gamma-1.)));
			 PNU=CL+CR-0.5*(Gamma-1.)*(UR-UL);
			 PDE=CL/(pow(PL,0.5*(Gamma-1.)/Gamma))+CR/(pow(PR,0.5*(Gamma-1.)/Gamma));
			 PM=pow(PNU/PDE,2.*Gamma/(Gamma-1.));
			 }
	*/
	else if(PPV<PMIN)//No 2-GAMMA TRRS, We use PVRS too.
		PM=PPV;
	else//Two shock
		{
			GEL=sqrt((2./((GammaL+1.)*DL))/((GammaL-1.)/(GammaL+1.)*PL+PPV));
			GER=sqrt((2./((GammaR+1.)*DR))/((GammaR-1.)/(GammaR+1.)*PR+PPV));
			PM=(GEL*PL+GER*PR-(UR-UL))/(GEL+GER);
			PM=(eps>PM?eps:PM);
		}
	*PStar = PM;
}
void PreFun(double *F,double *FD,double P,double DK,double PK,double CK,double Gamma)
// evaluate the pressure functions FL and FR in exact Riemann solver
{
	double AK,BK,PRAT,QRT;
	if(P<=PK)//Rarefaction wave
		{
			PRAT=P/PK;
			*F=2./(Gamma-1.)*CK*(pow(PRAT,(Gamma-1.)/(2.*Gamma))-1.);
			*FD=1./(DK*CK)*pow(PRAT,-(Gamma+1.)/(2.*Gamma));
		}
	else//Shock wave
		{
			AK=2./((Gamma+1.)*DK);
			BK=PK*(Gamma-1.)/(Gamma+1.);
			QRT=sqrt(AK/(BK+P));
			*F=(P-PK)*QRT;
			*FD=(1.-0.5*(P-PK)/(BK+P))*QRT;
		}
}
void StarPU(double *POut,double *U,double *DML,double*DMR,double DL,double DR,double UL,double UR,double PL,double PR,double CL,double CR,double GammaL, double GammaR)
// compute the solution for pressure and velocity in the star region
{
	const double eps = config[4];
	double P;
	double change,FL,FR,FLD,FRD,POLD,PSTART,UDIFF;
	int i,NRITER=100;
	GuessP(&PSTART,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
	POLD=PSTART;
	UDIFF=UR-UL;
	for(i=1;i<=NRITER;i++)
		{
			PreFun(&FL,&FLD,POLD,DL,PL,CL,GammaL);
			PreFun(&FR,&FRD,POLD,DR,PR,CR,GammaR);
			P=POLD-(FL+FR+UDIFF)/(FLD+FRD);
			change=2.*fabs((P-POLD)/(P+POLD));
			if(change<=eps)//compute velocity in star region
				{
					*U=0.5*(UL+UR+FR-FL);
					break;
				}
			if(P<0)
				P=eps;
			POLD=P;
		}
	if(P<=PL)//Left rarefaction wave
		*DML=DL*pow(P/PL,1./GammaL);
	else//Left shock wave
		*DML=DL*(P/PL+(GammaL-1.)/(GammaL+1.))/(P/PL*(GammaL-1.)/(GammaL+1.)+1.);
	if(P<=PR)//Right rarefaction wave
		*DMR=DR*pow(P/PR,1./GammaR);
	else//Right shock wave
		*DMR=DR*(P/PR+(GammaR-1.)/(GammaR+1.))/(P/PR*(GammaR-1.)/(GammaR+1.)+1.);
	*POut = P;
}
void AcousticSLag(double *DtDL,double *DtDR,double *DtU,double *DtP,double DUL,double DUR,double DPL,double DPR,
					double D,double U,double C_star,double r)
//GRP solver for acoustic case. Lagrangian version(moving mesh) for cylindrical case
{
	*DtP= 0.5*C_star*(DPR-DPL)-D*C_star*C_star*(0.5*(DUL+DUR)+(M-1)*U/r);
	*DtU=-0.5*(DPR+DPL)/D+0.5*C_star*(DUR-DUL);
	*DtDL=*DtP/C_star/C_star;
	*DtDR=*DtDL;
}
void AcousticSLagTangent(double *DtP,double *DtU,double DU,double DP,double D,double U,double C_star,double r)
//GRP solver for tangential case. Lagrangian version(moving mesh) for cylindrical case
{
	*DtP= -D*C_star*C_star*(DU+(M-1)*U/r);
	*DtU= -DP/D;
}
void GRPsolverSLag(double *DtDL,double *DtDR,double *DtU,double *DtP,double UM,double PM,
					 double DL,double DR,double UL,double UR,double PL,double PR,
					 double DDL,double DDR,double DUL,double DUR,double DPL,double DPR,
					 double TDSL,double TDSR,double DpsiL,double DphiR,double r,double GammaL,double GammaR)
//GRP solver for non acoustic case. Lagrangian version(moving mesh) cylindrical case
{
	const double eps = config[4];
	double aL,bL,dL,aR,bR,dR,CL,CR,C_starL,C_starR,D,U,P,theta,phic,phi1,phi2,phi3,sigmaL,sigmaR,musL,musR;
	musL=(GammaL-1.)/(GammaL+1.);
	musR=(GammaR-1.)/(GammaR+1.);
	CL=sqrt(GammaL*PL/DL);
	CR=sqrt(GammaR*PR/DR);
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
	*DtU=(dL*bR-dR*bL)/(aL*bR-aR*bL);
	*DtP=(dL*aR-dR*aL)/(bL*aR-bR*aL);
	*DtDL=*DtP/C_starL/C_starL;
	*DtDR=*DtP/C_starR/C_starR;
}
