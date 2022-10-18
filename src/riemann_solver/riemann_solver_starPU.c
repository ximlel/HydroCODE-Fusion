/**
 * @file  riemann_solver_starPU.c
 * @brief This is an exact two-component Riemann solver in Toro's book.
 */
#include <stdio.h>
#include <math.h>
#include <stdbool.h>


/**
 * @brief Provide a guess value for pressure PM in the Star Region.
 * @details The choice is made according to adaptive Riemann solver using the
 *          PVRS (or TRRS or TSRS) ApproxImate Riemann Solvers (AIRS).
 * @param[in]  DL, UL, PL, CL: Initial Density/Velocity/Pressure/Sound_speed on left  state.
 * @param[in]  DR, UR, PR, CR: Initial Density/Velocity/Pressure/Sound_speed on right state.
 * @param[in]  GammaL, GammaR: Ratio of specific heats.
 * @param[in]  eps: The largest value can be seen as zero.
 * @return   \b PM: A guess value for pressure in the Star Region.
 */
static double GuessP(const double DL, const double DR, const double UL, const double UR, const double PL, const double PR,
		     const double CL, const double CR, const double GammaL, const double GammaR, const double eps)
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


/**
 * @brief Evaluate the pressure functions FL and FR in exact Riemann solver.
 * @param[out] F:    Pressure function FL or FR.
 * @param[out] FD:   First derivative of F with respect to pressure P.
 * @param[in]  P:    Pressure.
 * @param[in]  DK:   Density.
 * @param[in]  PK:   Pressure.
 * @param[in]  CK:   Sound speed.
 * @param[in] Gamma: Ratio of specific heats.
 */
static void PreFun(double *F, double *FD, const double P, const double DK, const double PK, const double CK, const double Gamma)
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


/**
 * @brief EXACT RIEMANN SOLVER FOR Two-Component γ-Law Gas
 * @details The purpose of this function is to compute the Riemann solution for pressure and velocity in the Star Region,
 *          for the time dependent one dimensional Euler equations for two-component γ-law gas.
 * @param[out] U_star, P_star: Velocity/Pressure in star region.
 * @param[in]  UL, PL, CL:     Initial Velocity/Pressure/Sound_speed on left  state.
 * @param[in]  UR, PR, CR:     Initial Velocity/Pressure/Sound_speed on right state.
 * @param[in]  GammaL, GammaR: Ratio of specific heats.
 * @param[out] CRW: Centred Rarefaction Wave (CRW) Indicator of left and right waves.
 *                  - true: CRW
 *                  - false: Shock wave
 * @param[in]  eps:    The largest value can be seen as zero.
 * @param[in]  TOLPRE: Condition value of 'gap' at the end of the iteration.
 * @param[in]  NRITER: Maximum iteration step (Recommended Value: 100).
 * @return  \b change: Relative pressure change after the last iteration.
 * @author E. F. Toro
 * @date February 1st 1999
 * @sa   Theory is found in Chapter 4 of Reference [1]. \n
 *       [1] E. F. Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics". 
 *           Springer-Verlag, Second Edition, 1999
 * @copyright This program is part of NUMERICA —— \n
 *            A Library of Source Codes for Teaching, Research and Applications, by E. F. Toro \n
 *            Published by NUMERITEK LTD
 */
double Riemann_solver_starPU(double * U_star, double * P_star, const double GammaL, const double GammaR,
	    const double UL, const double UR, const double PL, const double PR,
	    const double CL, const double CR, _Bool * CRW,
	    const double eps, const double TOLPRE, const int NRITER)
{
	if((2.0*CL/(GammaL-1.0)+2.0*CR/(GammaR-1.0)) <= (UR-UL))
		printf("Error: Vacuum is generated in Riemann solver!\n");

	double DL=GammaL*PL/CL/CL;
	double DR=GammaR*PR/CR/CR;
	double P=0.5*(PL+PR), U=0.5*(UL+UR);
	double change = 0.0,FL,FR,FLD,FRD,POLD,PSTART,UDIFF;
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

	return change;
} 
