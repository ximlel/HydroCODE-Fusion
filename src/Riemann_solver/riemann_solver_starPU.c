#include <stdio.h>
#include <math.h>
#include <stdbool.h>


static double GuessP(const double DL, const double DR, const double UL, const double UR, const double PL, const double PR,
		     const double CL, const double CR, const double GammaL, const double GammaR, const double eps)
//provide a guess value for pressure PM in the star region.AIRS approximate Riemann solvers.
{
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
	return PM;
}


static void PreFun(double *F, double *FD, const double P, const double DK, const double PK, const double CK, const double Gamma)
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

// NRITER=100
void Riemann_solver_starPU(double * U_star, double * P_star, const double GammaL, const double GammaR,
	    const double UL, const double UR, const double PL, const double PR,
	    const double CL, const double CR, _Bool * CRW,
	    const double eps, const double TOLPRE, const int NRITER)
// compute the solution for pressure and velocity in the star region
{
	double DL=GammaL*PL/CL/CL;
	double DR=GammaR*PR/CR/CR;
	double P, U;
	double change,FL,FR,FLD,FRD,POLD,PSTART,UDIFF;
	PSTART = GuessP(DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR,eps);
	POLD=PSTART;
	UDIFF=UR-UL;
	for(int i=1;i<=NRITER;i++)
		{
			PreFun(&FL,&FLD,POLD,DL,PL,CL,GammaL);
			PreFun(&FR,&FRD,POLD,DR,PR,CR,GammaR);
			P=POLD-(FL+FR+UDIFF)/(FLD+FRD);
			change=2.*fabs((P-POLD)/(P+POLD));
			if(change<=TOLPRE)//compute velocity in star region
				{
					U=0.5*(UL+UR+FR-FL);
					break;
				}
			if(P<0)
				P=TOLPRE;
			POLD=P;
		}

	if(P<=PL) //Left rarefaction wave
	    CRW[0]=true;
	else //Left shock wave
	    CRW[0]=false;
	if(P<=PR) //Right rarefaction wave
	    CRW[1]=true;
	else //Right shock wave
	    CRW[1]=false;

	* U_star = U;
	* P_star = P;
} 
