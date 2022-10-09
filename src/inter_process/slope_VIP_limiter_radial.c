#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include_cpp/inter_process.hpp"


void VIP_limiter_radial(const int Ncell, const _Bool i_f_var_get, double DmU[], double TmV[],
		       const double UU[], struct radial_mesh_var *smv)
{
    double const dtheta = config[11]; //initial d_angle
    double const Alpha  = config[41]; // the paramater in slope limiters.
    /*
     * - LIMITER<0: add VIP limiter,
     * - LIMITER>0: only minmod limiter;
     * - abs(LIMITER)=1: original minmod limiter,
     * - abs(LIMITER)=2: VIP-like minmod limiter.
     */
    int const LIMITER_VIP = (int)config[42];
    double * Rb   = smv->Rb;
    double * RR   = smv->RR;
    double * DdrL = smv->DdrL;
    double * DdrR = smv->DdrR;
    double * dRc  = smv->dRc;

    double sU, sV;
    double VIP_lim, Vave[4][2], V0[2], Vp1[2], Vp2[2], Vp3[2];//VIP limiter

    //VIP limiter update
    for(int i = 1; i <= Ncell; i++)
	{
	    //sV=0.0;
	    //sV=VLmin[i]/(0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta));
	    sV=UU[i]/(0.5*(Rb[i]+Rb[i+1]));
	    sU=DmU[i];
	    Vave[0][0] = UU[i+1];
	    Vave[0][1] = 0.0;
	    Vave[1][0] = UU[i-1];
	    Vave[1][1] = 0.0;
	    Vave[2][0] = UU[i]*cos(dtheta);
	    Vave[2][1] = UU[i]*sin(dtheta);
	    Vave[3][0] = UU[i]*cos(dtheta);
	    Vave[3][1] =-UU[i]*sin(dtheta);
	    V0[0] = UU[i];
	    V0[1] = 0.0;
	    Vp1[0] = UU[i]+(0.5*(Rb[i]+Rb[i+1])-RR[i])*sU;
	    Vp1[1] = 0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta)*sV;
	    Vp2[0] = UU[i]+DdrL[i]*sU;
	    Vp2[1] = 0.0;
	    Vp3[0] = UU[i]-DdrR[i]*sU;
	    Vp3[1] = 0.0;
	    VIP_lim = fmin(1.0,     useVIPLimiter(4, Vave, V0, Vp1));
	    VIP_lim = fmin(VIP_lim, useVIPLimiter(4, Vave, V0, Vp2));
	    VIP_lim = fmin(VIP_lim, useVIPLimiter(4, Vave, V0, Vp3));
	    DmU[i]=VIP_lim*sU;
	    TmV[i]=VIP_lim*sV;
	    if (abs(LIMITER_VIP)==1)
		{
		    if(LIMITER_VIP>0)
			DmU[i]=minmod3(Alpha*(UU[i]-UU[i-1])/dRc[i],sU,Alpha*(UU[i+1]-UU[i])/dRc[i+1]);
		}
	    else if (abs(LIMITER_VIP)==2)
		{
		    if(LIMITER_VIP>0)
			DmU[i]=minmod3(Alpha*(UU[i]-UU[i-1])/2./DdrR[i],sU,Alpha*(UU[i+1]-UU[i])/2./DdrL[i]);
		}
	    if (LIMITER_VIP>0)
		TmV[i]=minmod2(Alpha*(UU[i]*sin(dtheta))/2./(0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta)),sV);
	}
    //sV=0.0;
    //sV=VLmin[1]/(0.5*Rb[1]*tan(0.5*dtheta));
    sV=UU[0]/Rb[1];
    sU=DmU[0];
    Vave[0][0] = UU[1];
    Vave[0][1] = 0.0;
    Vave[1][0] = UU[0]*cos(dtheta);
    Vave[1][1] = UU[0]*sin(dtheta);
    Vave[2][0] = UU[0]*cos(dtheta);
    Vave[2][1] =-UU[0]*sin(dtheta);
    V0[0] = UU[0];
    V0[1] = 0.0;
    Vp1[0] = UU[0]+(0.5*Rb[1]-RR[0])*sU;
    Vp1[1] = 0.5*Rb[1]*tan(0.5*dtheta)*sV;
    Vp2[0] = UU[0]+DdrL[0]*sU;
    Vp2[1] = 0.0;
    VIP_lim = fmin(1.0,     useVIPLimiter(3, Vave, V0, Vp1));
    VIP_lim = fmin(VIP_lim, useVIPLimiter(3, Vave, V0, Vp2));
    DmU[0]=VIP_lim*sU;
    TmV[0]=VIP_lim*sV;
    if (LIMITER_VIP>0)
	{
	    DmU[0]=minmod2(sU,DmU[1]);
	    TmV[0]=minmod2(sV,TmV[1]);
	}
    DmU[Ncell+1]=minmod2(DmU[Ncell],DmU[Ncell+1]); 
    TmV[Ncell+1]=TmV[Ncell];
}
