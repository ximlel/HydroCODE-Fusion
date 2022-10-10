#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "../include_cii/mem.h"
#include "../include/var_struc.h"


struct radial_mesh_var radial_mesh_init(const char *example)
{
    double const dr     = config[10];       // initial d_raidus
    double const dtheta = config[11];       // initial d_angle
    double const r_0    = config[20];
    int    const Md     = (int)config[3]+2; // max vector dimension

    struct radial_mesh_var smv;
    smv.Rb   = (double*)ALLOC(Md*sizeof(double)); //radius and length of outer cell boundary
    smv.Lb   = (double*)ALLOC(Md*sizeof(double));
    smv.RR   = (double*)ALLOC(Md*sizeof(double)); //centroidal radius and variable in cells
    smv.DdrL = (double*)ALLOC(Md*sizeof(double)); //distance from boundary to center in a cell
    smv.DdrR = (double*)ALLOC(Md*sizeof(double));
    smv.Ddr  = (double*)ALLOC(Md*sizeof(double));
    smv.dRc  = (double*)ALLOC(Md*sizeof(double)); //(derivative)centers distance
    smv.vol  = (double*)ALLOC(Md*sizeof(double));

	smv.Rb[0]=0.0;
	smv.Lb[0]=0.0;
	smv.RR[0]=(2.0/3.0*r_0)*cos(0.5*dtheta);
	for(int i = 1; i < Md; i++)//center cell is cell 0
		{
			smv.Rb[i]=(r_0+(i-1)*dr)*cos(0.5*dtheta);//outer cell boundary
			smv.Lb[i]=(r_0+(i-1)*dr)*sin(0.5*dtheta)*2.0;
			//centroid of triangular and trapezoid
			smv.RR[i]=(r_0+i*dr-(3*r_0/dr+3*i-2)/(6*r_0/dr+6*i-3)*dr)*cos(0.5*dtheta);
		    //smv.RR[i]=smv.Rb[i+1]-(2.0*smv.Lb[i]+smv.Lb[i+1])/(3.0*(smv.Lb[i]+smv.Lb[i+1]))*(smv.Rb[i+1]-smv.Rb[i]);
		}

	printf("Mesh of %s has been constructed!\n", example);
	return smv;
}


void radial_mesh_update(struct radial_mesh_var *smv)
{
    int const Ncell = (int)config[3]; // Number of computing cells in r direction

    double *Rb   = smv->Rb;  //radius and length of outer cell boundary
    double *Lb   = smv->Lb;
    double *RR   = smv->RR;  //centroidal radius and variable in cells
    /* Variables to be updated  */
    double *DdrL = smv->DdrL;//distance from boundary to center in a cell
    double *DdrR = smv->DdrR;
    double *Ddr  = smv->Ddr;
    double *dRc  = smv->dRc; //(derivative)centers distance
    double *vol  = smv->vol;

	for(int i = 0; i <= Ncell; i++)
		{
			dRc[i]  = RR[i]-RR[i-1];
			DdrL[i] = Rb[i+1]-RR[i]; // Right side length of cell i
			DdrR[i] = RR[i]-Rb[i];   // Left  side length of cell i
			Ddr[i]  = DdrL[i]+DdrR[i];
			if(Ddr[i] < 0.0)
				{
					fprintf(stderr, "ERROR! deltar_r < 0 in cell %d.\n", i);
					exit(3);
				}
			vol[i] = RR[i]*0.5*(Lb[i]+Lb[i+1])*Ddr[i]; //m=2.
		}
	dRc[Ncell+1]  = RR[Ncell+1]-RR[Ncell]; // boundary condition
	DdrR[Ncell+1] = RR[Ncell+1]-Rb[Ncell+1];
	Ddr[Ncell+1]  = Ddr[Ncell];
}


void radial_mesh_mem_free(struct radial_mesh_var *smv)
{
    FREE(smv->Rb);
    FREE(smv->Lb);
    FREE(smv->RR);
    FREE(smv->DdrL);
    FREE(smv->DdrR);
    FREE(smv->Ddr);
    FREE(smv->dRc);
    FREE(smv->vol);
}
