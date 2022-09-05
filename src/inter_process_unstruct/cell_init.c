#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/inter_process_unstruct.h"


#define FV_RESET_MEM(v, n)						\
    do {								\
	FV->v = realloc(FV->v, (n) * sizeof(double));			\
	if(FV->v == NULL)						\
	    {								\
		fprintf(stderr, "Not enough memory in fluid var init!\n"); \
		goto return_NULL;					\
	    }								\
    } while (0)								\
		

#define CV_INIT_MEM(v, n)						\
    do {								\
	cv.v = calloc(n, sizeof(double));				\
	if(cv.v == NULL)						\
	    {								\
		fprintf(stderr, "Not enough memory in cell var init!\n"); \
		goto return_NULL;					\
	    }								\
    } while (0)								\

#define CP_INIT_MEM(v, n)						\
    do {								\
	cv.v = malloc((n) * sizeof(void *));				\
	if(cv.v == NULL)						\
	    {								\
		fprintf(stderr, "Not enough memory in cell var init!\n"); \
		goto return_NULL;					\
	    }								\
	init_mem(cv.v, n, mv->cell_pt);					\
    } while (0)								\

#define CP_INIT_MEM_INT(v, n)						\
    do {								\
	cv.v = malloc((n) * sizeof(void *));				\
	if(cv.v == NULL)						\
	    {								\
		fprintf(stderr, "Not enough memory in cell var init!\n"); \
		goto return_NULL;					\
	    }								\
	init_mem_int(cv.v, n, mv->cell_pt);				\
    } while (0)								\

struct cell_var cell_mem_init(const struct mesh_var * mv, struct flu_var * FV)
{
	const int order = (int)config[9];
	const int num_cell_ghost = mv->num_ghost + (int)config[3];
	const int num_cell = (int)config[3];
	
	struct cell_var cv;

	CP_INIT_MEM_INT(cell_cell, num_cell_ghost);
	CP_INIT_MEM(n_x, num_cell_ghost);
	CP_INIT_MEM(n_y, num_cell_ghost);
	CV_INIT_MEM(X_c, num_cell_ghost);
	CV_INIT_MEM(Y_c, num_cell_ghost);
	CV_INIT_MEM(vol, num_cell_ghost);

	CP_INIT_MEM(F_u,   num_cell);
	CP_INIT_MEM(F_v,   num_cell);
	CP_INIT_MEM(F_rho, num_cell);
	CP_INIT_MEM(F_e,   num_cell);
	CV_INIT_MEM(U_u,   num_cell_ghost);
	CV_INIT_MEM(U_v,   num_cell_ghost);
	CV_INIT_MEM(U_rho, num_cell_ghost);
	CV_INIT_MEM(U_e,   num_cell_ghost);
	FV_RESET_MEM(U,    num_cell_ghost);
	FV_RESET_MEM(V,    num_cell_ghost);
	FV_RESET_MEM(RHO,  num_cell_ghost);
	FV_RESET_MEM(P,    num_cell_ghost);

	CP_INIT_MEM(U_p,   num_cell_ghost);
	CP_INIT_MEM(V_p,   num_cell_ghost);
	CP_INIT_MEM(P_p,   num_cell);
	CP_INIT_MEM(RHO_p, num_cell);
	CP_INIT_MEM(F_p_x, num_cell);
	CP_INIT_MEM(F_p_y, num_cell);
	if (order > 1)
		{
			CV_INIT_MEM(gradx_rho, num_cell_ghost);
			CV_INIT_MEM(grady_rho, num_cell_ghost);
			CV_INIT_MEM(gradx_e,   num_cell_ghost);
			CV_INIT_MEM(grady_e,   num_cell_ghost);
			CV_INIT_MEM(gradx_u,   num_cell_ghost);			
			CV_INIT_MEM(grady_u,   num_cell_ghost);
			CV_INIT_MEM(gradx_v,   num_cell_ghost);
			CV_INIT_MEM(grady_v,   num_cell_ghost);

			CP_INIT_MEM(dt_U_p,    num_cell_ghost);
			CP_INIT_MEM(dt_V_p,    num_cell_ghost);
			CP_INIT_MEM(dt_F_p_x,  num_cell);
			CP_INIT_MEM(dt_F_p_y,  num_cell);
		}

#ifdef MULTIPHASE_BASICS
	CP_INIT_MEM(F_phi, num_cell);
	CV_INIT_MEM(U_phi, num_cell_ghost);
	FV_RESET_MEM(PHI, num_cell_ghost);
	CP_INIT_MEM(F_e_a, num_cell);
	CV_INIT_MEM(U_e_a, num_cell_ghost);
	FV_RESET_MEM(Z_a, num_cell_ghost);
	CP_INIT_MEM(F_gamma, num_cell);
	CV_INIT_MEM(U_gamma, num_cell_ghost);
	FV_RESET_MEM(gamma, num_cell_ghost);
	CP_INIT_MEM(PHI_p, num_cell);
	CP_INIT_MEM(Z_a_p, num_cell);
	CP_INIT_MEM(gamma_p, num_cell);
	if (order > 1)
	    {
		CV_INIT_MEM(gradx_z_a, num_cell_ghost);
		CV_INIT_MEM(grady_z_a, num_cell_ghost);
		CV_INIT_MEM(gradx_phi, num_cell_ghost);
		CV_INIT_MEM(grady_phi, num_cell_ghost);
		if ((_Bool)config[60])
		    {
			CV_INIT_MEM(gradx_gamma, num_cell_ghost);
			CV_INIT_MEM(grady_gamma, num_cell_ghost);
		    }
	    }
#endif

#ifdef MULTIPHASE_BASICS
	CP_INIT_MEM(P_star, num_cell);
	CP_INIT_MEM(U_qt_star, num_cell);
	CP_INIT_MEM(V_qt_star, num_cell);
	CP_INIT_MEM(U_qt_add_c, num_cell);
	CP_INIT_MEM(V_qt_add_c, num_cell);
#endif
#ifdef MULTIPHASE_BASICS_abandoned
	CP_INIT_MEM(RHO_star, num_cell);
	CP_INIT_MEM(gamma_star, num_cell);

	CP_INIT_MEM(RHO_minus_c, num_cell);
	CP_INIT_MEM(P_minus_c, num_cell);
	CP_INIT_MEM(U_qt_minus_c, num_cell);
	CP_INIT_MEM(V_qt_minus_c, num_cell);
	CP_INIT_MEM(gamma_minus_c, num_cell);

	CP_INIT_MEM(RHO_add_c, num_cell);
	CP_INIT_MEM(P_add_c, num_cell);
	CP_INIT_MEM(gamma_add_c, num_cell);

	CP_INIT_MEM(u_star, num_cell);
	CP_INIT_MEM(u_minus_c, num_cell);
	CP_INIT_MEM(u_add_c, num_cell);
#endif
	return cv;
	
 return_NULL:
	exit(5);
}


//Calculate volume.
void vol_comp(const struct cell_var * cv, const struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];
	int **cp = mv->cell_pt;
	
	int p_p, p_n;

	for(int k = 0; k < num_cell; k++)
	    {			
		cv->vol[k] = 0.0;
		for(int j = 0; j < cp[k][0]; j++)
		    {
			if(j == cp[k][0]-1) 
			    {
				p_p = cp[k][1];
				p_n = cp[k][j+1];
			    }				  
			else
			    {
				p_p = cp[k][j+2];
				p_n = cp[k][j+1];
			    } 
			cv->vol[k] = cv->vol[k] + 0.5 * (mv->X[p_n]*mv->Y[p_p] - mv->Y[p_n]*mv->X[p_p]);
		    }
	    }
}



//Determine the normal direction and relationship between cells.
void cell_rel(const struct cell_var * cv, const struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];
	
	int **cp = mv->cell_pt;
	int p_p, p_n, p2_p, p2_n;

	int cell_rec, n_border;
	int i, l, ts;
	double length;

	for(int k = 0; k < num_cell; k++)
	    { 						
		for(int j = 0; j < cp[k][0]; j++)
		    {
			if(j == cp[k][0]-1) 
			    {
				p_p = cp[k][1];
				p_n = cp[k][j+1];
			    }				  
			else
			    {
				p_p = cp[k][j+2];
				p_n = cp[k][j+1];
			    }
			length = sqrt((mv->Y[p_p] - mv->Y[p_n])*(mv->Y[p_p] - mv->Y[p_n])+(mv->X[p_n] - mv->X[p_p])*(mv->X[p_n] - mv->X[p_p]));
			cv->n_x[k][j] = (mv->Y[p_p] - mv->Y[p_n]) / length;
			cv->n_y[k][j] = (mv->X[p_n] - mv->X[p_p]) / length;
			//Inner normal 

			      cell_rec = 0;							   		
			ts = 1;
			while (ts <= MAX(num_cell-k-1, k))
			    {
				// seek in two side
				i = k + ts;
				if (ts > 0)
				    ts = -ts;
				else
				    ts = -ts + 1;
				if (i < 0 || i >= num_cell)
				    continue;
								
				for(l = 0; l < cp[i][0]; l++)
				    {
					if(l == cp[i][0]-1) 
					    {
						p2_p = cp[i][1];
						p2_n = cp[i][l+1];
					    }				  
					else
					    {
						p2_p = cp[i][l+2];
						p2_n = cp[i][l+1];
					    }
					if((p_p == p2_n) && (p2_p == p_n))
					    {
						cv->cell_cell[k][j] = i;
						cell_rec = 1;;
						break;
					    }
				    }
				if (cell_rec)
				    break;
			    }
				
			if (cell_rec)
			    continue;								
					
			for(l = 1, n_border = -1; l <= mv->num_border[0]; l++)
			    {
				n_border += mv->num_border[l] + 1; 
				for(i = n_border-mv->num_border[l]; i < n_border; i++)
				    {				
					p2_p = mv->border_pt[i+1];
					p2_n = mv->border_pt[i];
					if((p_p == p2_p && p_n == p2_n) || (p_p == p2_n && p_n == p2_p))
					    {
						cv->cell_cell[k][j] = mv->border_cond[i];
						cell_rec = 1;
						break;
					    }
				    }							
				if (cell_rec)
				    break;
			    }
			
			if(!cell_rec && k < (int)config[3])
			    {
				fprintf(stderr, "Ther are some wrong cell relationships!\n");
				exit(2);
			    }
		    }							
	    }
}


void cell_centroid(const struct cell_var * cv, const struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];
	const double *X = mv->X, *Y = mv->Y;
	int **cp = mv->cell_pt;
	
	double S, S_tri;

	for(int k = 0; k < num_cell; ++k)
	    {
		S = 0.0;
		cv->X_c[k] = 0.0;
		cv->Y_c[k] = 0.0;

		for(int j = 2; j < cp[k][0]; j++)
		    {
			S_tri = X[cp[k][1]]*Y[cp[k][j]] + X[cp[k][j+1]]*Y[cp[k][1]] + X[cp[k][j]]*Y[cp[k][j+1]] - X[cp[k][j+1]]*Y[cp[k][j]] - X[cp[k][1]]*Y[cp[k][j+1]] - X[cp[k][j]]*Y[cp[k][1]];
			cv->X_c[k] += (X[cp[k][1]] + X[cp[k][j]] + X[cp[k][j+1]]) * S_tri;
			cv->Y_c[k] += (Y[cp[k][1]] + Y[cp[k][j]] + Y[cp[k][j+1]]) * S_tri;
			S += S_tri;
		    }			 
		cv->X_c[k] /= S*3.0;
		cv->Y_c[k] /= S*3.0;
	    }
}
