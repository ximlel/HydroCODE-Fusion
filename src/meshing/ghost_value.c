#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"



#define CV_COPY(var)  cv->var[i] = cv->var[pc[i]]
#define FV_COPY(var)  FV->var[i] = FV->var[pc[i]]

void period_ghost(struct cell_var * cv, const struct mesh_var * mv, struct flu_var * FV, double t)
{
	const int order = (int)config[9];
	const int num_cell = (int)config[3];
	const int num_cell_ghost = mv->num_ghost + num_cell;
	const int *pc = mv->peri_cell;	
	
	for(int i = num_cell; i < num_cell_ghost; i++)
		{
			CV_COPY(U_rho);			
			CV_COPY(U_e);
			CV_COPY(U_u);
			CV_COPY(U_v);
			FV_COPY(RHO);
			FV_COPY(P);
			FV_COPY(U);
			FV_COPY(V);
			if (order > 1)
			    {
				CV_COPY(gradx_rho);
				CV_COPY(gradx_e);
				CV_COPY(gradx_u);
				CV_COPY(gradx_v);
				CV_COPY(grady_rho);
				CV_COPY(grady_e);
				CV_COPY(grady_u);
				CV_COPY(grady_v);
			    }
#ifdef MULTIFLUID_BASICS
			CV_COPY(U_e_a);
			CV_COPY(U_phi);
			CV_COPY(U_gamma);
			FV_COPY(PHI);
			FV_COPY(gamma);
			if (order > 1)
			    {
				CV_COPY(gradx_phi);
				CV_COPY(grady_phi);
			    }
			FV_COPY(Z_a);
			if (order > 1)
			    {
				CV_COPY(gradx_z_a);
				CV_COPY(grady_z_a);
			    }
#endif
		}
}



void period_cell_modi(struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];
	int *pc = mv->peri_cell;
	
	int *per_num = malloc(num_cell*sizeof(int));
	int per_n = 0;
	int i, j;
	for (i = 0; i < num_cell; i++)
		{
			if (pc[i] >= 0)
				per_n++;
			per_num[i] = per_n;
		}

	for (i = 0; i < num_cell; i++)
		if (pc[i] >= 0)
			{
				pc[i] -= per_num[pc[i]];
				if (pc[i] < 0)
					pc[i] = 0;
			}

	int *cc_tmp, pc_tmp; 
	for (i = 1; i < num_cell; i++)
		{
			j = i;
			for(j = i; pc[j-1] >= 0 && j >= 1; j--)
				{
					if(pc[j] < 0)
						{
							pc_tmp = pc[j-1];
							pc[j-1] = pc[j];
							pc[j] = pc_tmp;
							cc_tmp = mv->cell_pt[j];
							mv->cell_pt[j] = mv->cell_pt[j-1];
							mv->cell_pt[j-1] = cc_tmp;
						}
					else
						break;
				}	 				
		}
	
	mv->bc = period_ghost;
}
