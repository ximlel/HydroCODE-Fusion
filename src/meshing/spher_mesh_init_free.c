#include <math.h>
#include <stdlib.h>

#include "../include/var_struc.h"

#define Epsilon (1.) // r_0=Epsilon*dr


struct spher_mesh_var spher_mesh_init(void)
{
    int    const Ncell  = (int)config[13]; // Number of computing cells in r direction
    double const dr     = config[10]; //initial d_raidus
    double const dtheta = config[11]; //initial d_angle

    struct spher_mesh_var mv = {NULL};
    double * Rb   = mv.Rb;
    double * Lb   = mv.Lb;
    double * RR   = mv.RR;
    double * DdrL = mv.DdrL;
    double * DdrR = mv.DdrR;
    double * Ddr  = mv.Ddr;
    double * dRc  = mv.dRc;
    double * vol  = mv.vol;

    int i;
	for(i=1;i<=Ncell;i++)//center cell is cell 0
		{
			Rb[i]=(Epsilon+i-1.)*dr*cos(0.5*dtheta);//outer cell boundary
			Lb[i]=2.*(Epsilon+i-1.)*dr*sin(0.5*dtheta);
			RR[i]=(Epsilon+i-(3.*Epsilon+3.*i-2.)/(6.*Epsilon+6.*i-3.))*dr*cos(0.5*dtheta);
			//centroid of triangular and trapezoid
		}
	Rb[0]=0.;
	Lb[0]=0.;
	RR[0]=(2./3.*Epsilon)*dr*cos(0.5*dtheta);

	for(i=1;i<Ncell;i++)
		{
			DdrL[i]=Rb[i+1]-RR[i];
			DdrR[i]=RR[i]-Rb[i];
			Ddr[i] =DdrL[i]+DdrR[i];
			dRc[i] =RR[i]-RR[i-1];
			vol[i] =RR[i]*0.5*(Lb[i]+Lb[i+1])*(Rb[i+1]-Rb[i]);//m=2.
		}
	DdrL[0]=Rb[1]-RR[0];
	Ddr[0] =Rb[1];
	dRc[0] =DdrL[0];
	vol[0] =RR[0]*0.5*Lb[1]*Rb[1];//m=2.
	DdrR[Ncell]=DdrR[Ncell-1];
	Ddr[Ncell] =Ddr[Ncell-1];
	dRc[Ncell] =DdrL[Ncell-1];//boundary condition
	/*
	  DdrR[Ncell]=RR[Ncell]-Rb[Ncell];
	  Ddr[Ncell] =Ddr[Ncell-1];
	  dRc[Ncell] =RR[Ncell]-RR[Ncell-1];//boundary condition
	*/
    return mv;
}
