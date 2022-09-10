/**
 * @file  ghost_cell.c
 * @brief This is a set of functions which manipulate the ghost cell and the data on it.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include_cii/mem.h"
#include "../include/var_struc.h"


/**
 * @brief Copy the grid variable data on the corresponding cell to the i-th ghost cell.
 */
#define CV_COPY(var)  cv->var[i] = cv->var[pc[i]]
/**
 * @brief Copy the fluid variable data on the corresponding cell to the i-th ghost cell.
 */
#define FV_COPY(var)  FV->var[i] = FV->var[pc[i]]


/**
 * @brief Copy the grid and fluid variable data in struct 'cv' and 'FV' on the corresponding cell to the i-th ghost cell.
 * @param[in] cv: Structure of grid variable data in computational grid cells.
 * @param[in] mv: Structure of meshing variable data.
 * @param[in] FV: Structure of fluid variable data array pointer.
 * @param[in] t:  Current computational time.
 */
void period_ghost(struct cell_var * cv, const struct mesh_var * mv, struct flu_var * FV, const double t)
{
	const int order = (int)config[9];
	const int num_cell = (int)config[3];
	const int num_cell_ghost = mv->num_ghost + num_cell;
	const int *pc = mv->period_cell;
	
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


/**
 * @brief Recount 'mv->cell_pt' according to the grid cell on the periodic boundary 'mv->period_cell'.
 * @param[in] mv: Structure of meshing variable data.
 */
void period_cell_modify(struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];
	int *pc = mv->period_cell;
	
	int *per_num = (int*)ALLOC(num_cell*sizeof(int));
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
	FREE(per_num);

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
