#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"
#include "../include/meshing.h"



static void cell_pt_clockwise(const struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];
	const double *X = mv->X, *Y = mv->Y;
	int **cp = mv->cell_pt;
	int p_p, p, p_n;

	int X_max, n_max, tmp;
	for(int k = 0; k < num_cell; k++)
	    {
		n_max = 1;
		p = cp[k][n_max];
		X_max = X[p];			
		for(int j = 2; j <= cp[k][0]; j++)
		    {
			n_max = X[cp[k][j]] > X_max ? j : n_max;
			p = cp[k][n_max];
			X_max = X[p];
		    }

		if(n_max == cp[k][0]) 
		    {
			p_p=cp[k][1];
			p_n=cp[k][n_max-1];
		    }
		else if(n_max == 1)
		    {
			p_p=cp[k][n_max+1];
			p_n=cp[k][cp[k][0]];
		    }
		else
		    {
			p_p=cp[k][n_max+1];
			p_n=cp[k][n_max-1];
		    }
			
		if ((X[p_p]-X[p])*(Y[p_n]-Y[p]) - (Y[p_p]-Y[p])*(X[p_n]-X[p]) < 0.0)
		    for(int j = 1; j < cp[k][0]/2; j++)
			{
			    tmp = cp[k][j];
			    cp[k][j] = cp[k][cp[k][0]+1-j];
			    cp[k][cp[k][0]+1-j] = tmp;
			}			
	    }
}


struct mesh_var mesh_load(const char *example, const char *mesh_name)
{	
	struct mesh_var mv = {0, 0, NULL, NULL, {1}, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

	char add_mkdir[FILENAME_MAX];
	example_io(example, add_mkdir, 1);
	char add[FILENAME_MAX];
	strcpy(add, add_mkdir);
	strcat(add, mesh_name);
	strcat(add, ".msh");

	FILE * fp;
	if ((fp = fopen(add, "r")) != NULL)
		{
		    if(msh_read(fp, &mv) == 0)
			{
			    fclose(fp);
			    exit(2);
			}
		    else
			{									
			    fclose(fp);
			    printf("Mesh file(%s.msh) has been read!\n", mesh_name);
			    return mv;
			}
		}
		
	if (strcmp(mesh_name,"Sod") == 0)
	    Sod_mesh(&mv);
	else if (strcmp(mesh_name,"Shear") == 0)
	    Shear_mesh(&mv);
	else if (strcmp(mesh_name,"Shock_Bubble") == 0)
	    Shock_Bubble_mesh(&mv);
	else if (strcmp(mesh_name,"free") == 0)
	    free_mesh(&mv);
	else if (strcmp(mesh_name,"Shell") == 0)
	    Shell_mesh(&mv);
	else if (strcmp(mesh_name,"RMI") == 0)
	    RMI_mesh(&mv);
	else if (strcmp(mesh_name,"RMI_S") == 0)
	    RMI_S_mesh(&mv);
	else if (strcmp(mesh_name,"R2D") == 0)
	    R2D_mesh(&mv);
	else if (strcmp(mesh_name,"Vortex") == 0)
	    Vortex_mesh(&mv);
	else if (strcmp(mesh_name,"cylinder") == 0)
	    cylinder_mesh(&mv);
	else if (strcmp(mesh_name,"odd_even") == 0)
	    odd_even_mesh(&mv);
	else if (strcmp(mesh_name,"odd_even_periodic") == 0)
	    odd_even_periodic_mesh(&mv);
	else if (strcmp(mesh_name,"odd_even_inflow") == 0)
	    odd_even_inflow_mesh(&mv);
	else if (strcmp(mesh_name,"rand_disturb_inflow") == 0)
	    rand_disturb_inflow_mesh(&mv);
	else if (strcmp(mesh_name,"oblique_periodic") == 0)
	    oblique_periodic_mesh(&mv);
	else if (strcmp(mesh_name,"Saltzman") == 0)
	    Saltzman_mesh_Lag(&mv);
	else
	    {
		fprintf(stderr, "No mesh setting!\n");
		exit(2);
	    }

	cell_pt_clockwise(&mv);
	
	return mv;
}

