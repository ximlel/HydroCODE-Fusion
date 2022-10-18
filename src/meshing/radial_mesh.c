/**
 * @file  radial_mesh.c
 * @brief There are some handler functions of radially symmetric mesh.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "../include_cii/mem.h"
#include "../include/var_struc.h"


/**
 * @brief This function initializes radially symmetric meshing variables.
 * @param[in]  example: Name of the test example.
 * @return  \b rmv:     Structure of radially symmetric meshing variable data.
 */
struct radial_mesh_var radial_mesh_init(const char *example)
{
    double const dr     = config[10];       // initial d_raidus
    double const dtheta = config[11];       // initial d_angle
    double const r_0    = config[20];
    int    const Md     = (int)config[3]+2; // max vector dimension

    struct radial_mesh_var rmv;
    rmv.Rb   = (double*)ALLOC(Md*sizeof(double)); //radius and length of outer cell boundary
    rmv.Lb   = (double*)ALLOC(Md*sizeof(double));
    rmv.RR   = (double*)ALLOC(Md*sizeof(double)); //centroidal radius and variable in cells
    rmv.DdrL = (double*)ALLOC(Md*sizeof(double)); //distance from boundary to center in a cell
    rmv.DdrR = (double*)ALLOC(Md*sizeof(double));
    rmv.Ddr  = (double*)ALLOC(Md*sizeof(double));
    rmv.dRc  = (double*)ALLOC(Md*sizeof(double)); //(derivative)centers distance
    rmv.vol  = (double*)ALLOC(Md*sizeof(double));

	rmv.Rb[0]=0.0;
	rmv.Lb[0]=0.0;
	rmv.RR[0]=(2.0/3.0*r_0)*cos(0.5*dtheta);
	for(int i = 1; i < Md; i++)//center cell is cell 0
		{
			rmv.Rb[i]=(r_0+(i-1)*dr)*cos(0.5*dtheta);//outer cell boundary
			rmv.Lb[i]=(r_0+(i-1)*dr)*sin(0.5*dtheta)*2.0;
			//centroid of triangular and trapezoid
			rmv.RR[i]=(r_0+i*dr-(3*r_0/dr+3*i-2)/(6*r_0/dr+6*i-3)*dr)*cos(0.5*dtheta);
		    //rmv.RR[i]=rmv.Rb[i+1]-(2.0*rmv.Lb[i]+rmv.Lb[i+1])/(3.0*(rmv.Lb[i]+rmv.Lb[i+1]))*(rmv.Rb[i+1]-rmv.Rb[i]);
		}

	printf("Mesh of %s has been constructed!\n", example);
	return rmv;
}


/**
 * @brief This function updates radially symmetric meshing variables after each time step update.
 * @param[in,out] rmv: Structure of radially symmetric meshing variable data.
 */
void radial_mesh_update(struct radial_mesh_var *rmv)
{
    int const Ncell = (int)config[3]; // Number of computing cells in r direction

    double *Rb   = rmv->Rb;  //radius and length of outer cell boundary
    double *Lb   = rmv->Lb;
    double *RR   = rmv->RR;  //centroidal radius and variable in cells
    /* Variables to be updated  */
    double *DdrL = rmv->DdrL;//distance from boundary to center in a cell
    double *DdrR = rmv->DdrR;
    double *Ddr  = rmv->Ddr;
    double *dRc  = rmv->dRc; //(derivative)centers distance
    double *vol  = rmv->vol;

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


/**
 * @brief This function free memory for storing radially symmetric meshing variables.
 * @param[out] rmv: Structure of radially symmetric meshing variable data.
 */
void radial_mesh_mem_free(struct radial_mesh_var *rmv)
{
    FREE(rmv->Rb);
    FREE(rmv->Lb);
    FREE(rmv->RR);
    FREE(rmv->DdrL);
    FREE(rmv->DdrR);
    FREE(rmv->Ddr);
    FREE(rmv->dRc);
    FREE(rmv->vol);
}
