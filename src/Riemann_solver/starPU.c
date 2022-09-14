#include <stdio.h>
#include <math.h>

#include "../include/var_struc.h"


static void GuessP(double *PStar,double DL,double DR,double UL,double UR,double PL,double PR,double CL,double CR,double GammaL, double GammaR)
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


static void PreFun(double *F,double *FD,double P,double DK,double PK,double CK,double Gamma)
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


void StarPU(double *U_star, const struct i_f_var *ifv_L, const struct i_f_var *ifv_R)
// compute the solution for pressure and velocity in the star region
{
	double GammaL = ifv_L->gamma, GammaR = ifv_R->gamma;
	double DL = ifv_L->RHO, DR = ifv_R->RHO;
	double UL = ifv_L->U,   UR = ifv_R->U;
	double PL = ifv_L->P,   PR = ifv_R->P;
	double CL, CR;
	CL=sqrt(GammaL*PL/DL);
	CR=sqrt(GammaR*PR/DR);
	double P, U, DML, DMR;
	const double eps = config[4];
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
					U=0.5*(UL+UR+FR-FL);
					break;
				}
			if(P<0)
				P=eps;
			POLD=P;
		}
	if(P<=PL) //Left rarefaction wave
		DML=DL*pow(P/PL,1./GammaL);
	else //Left shock wave
		DML=DL*(P/PL+(GammaL-1.)/(GammaL+1.))/(P/PL*(GammaL-1.)/(GammaL+1.)+1.);
	if(P<=PR) //Right rarefaction wave
		DMR=DR*pow(P/PR,1./GammaR);
	else //Right shock wave
		DMR=DR*(P/PR+(GammaR-1.)/(GammaR+1.))/(P/PR*(GammaR-1.)/(GammaR+1.)+1.);

	U_star[0] = DML;
	U_star[3] = DMR;
	U_star[1] = U;
	U_star[2] = P;
} 
