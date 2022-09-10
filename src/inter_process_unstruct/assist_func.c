#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/inter_process_unstruct.h"


int fluid_var_update(struct flu_var * FV, struct cell_var * cv)
{
	const int num_cell = (int)config[3];
	struct i_f_var ifv;

	for(int k = 0; k < num_cell; ++k)
		{
			cons_qty_copy_cv2ifv(&ifv, cv, k);			

			if(cons2prim(&ifv) == 0)
				{
					fprintf(stderr, "Wrong in copying cons_var to prim_var!\n");
					return 0;
				}
			prim_var_copy_ifv2FV(&ifv, FV, k);

			cons_qty_copy_ifv2cv(&ifv, cv, k);			
		}

	return 1;
}


static int order2_i_f_var_init(const struct cell_var * cv, struct i_f_var * ifv, const int k)
{
	const double n_x = ifv->n_x, n_y = ifv->n_y;
	const double delta_x = ifv->delta_x, delta_y = ifv->delta_y;

	ifv->d_rho  =  cv->gradx_rho[k]*n_x + cv->grady_rho[k]*n_y;
	ifv->d_e    =  cv->gradx_e[k]  *n_x + cv->grady_e[k]  *n_y;
	ifv->d_u    =  cv->gradx_u[k]  *n_x + cv->grady_u[k]  *n_y;
	ifv->d_v    =  cv->gradx_v[k]  *n_x + cv->grady_v[k]  *n_y;
	ifv->t_rho  = -cv->gradx_rho[k]*n_y + cv->grady_rho[k]*n_x;
	ifv->t_e    = -cv->gradx_e[k]  *n_y + cv->grady_e[k]  *n_x;
	ifv->t_u    = -cv->gradx_u[k]  *n_y + cv->grady_u[k]  *n_x;
	ifv->t_v    = -cv->gradx_v[k]  *n_y + cv->grady_v[k]  *n_x;
#ifdef MULTIFLUID_BASICS
	ifv->d_z_a  =  cv->gradx_z_a[k]*n_x + cv->grady_z_a[k]*n_y;
	ifv->t_z_a  = -cv->gradx_z_a[k]*n_y + cv->grady_z_a[k]*n_x;
	ifv->d_phi  =  cv->gradx_phi[k]*n_x + cv->grady_phi[k]*n_y;
	ifv->t_phi  = -cv->gradx_phi[k]*n_y + cv->grady_phi[k]*n_x;
#endif

	if (cons2prim(ifv) == 0)
		{
			fprintf(stderr, "Error happens on primitive variable!\n");
			return 0;
		}

	if ((int)config[31] == 0)
		{
			ifv->d_p = ifv->d_e;
			ifv->t_p = ifv->t_e;

			ifv->RHO += cv->gradx_rho[k]*delta_x + cv->grady_rho[k]*delta_y;
			ifv->P   += cv->gradx_e[k]  *delta_x + cv->grady_e[k]  *delta_y;
			ifv->U   += cv->gradx_u[k]  *delta_x + cv->grady_u[k]  *delta_y;
			ifv->V   += cv->gradx_v[k]  *delta_x + cv->grady_v[k]  *delta_y;
#ifdef MULTIFLUID_BASICS
			ifv->Z_a += cv->gradx_z_a[k]*delta_x + cv->grady_z_a[k]*delta_y;
			ifv->PHI += cv->gradx_phi[k]*delta_x + cv->grady_phi[k]*delta_y;
			ifv->gamma = 1.0+1.0/(ifv->Z_a/(config[6]-1.0)+(1.0-ifv->Z_a)/(config[106]-1.0));
#endif
		}
	else if ((int)config[31] == 1)
		{
			ifv->d_u = (ifv->d_u - ifv->U*ifv->d_rho)/ifv->RHO;
			ifv->d_v = (ifv->d_v - ifv->V*ifv->d_rho)/ifv->RHO;
			ifv->d_p = (ifv->d_e - 0.5*ifv->d_rho*ifv->U*ifv->U - ifv->RHO*ifv->U*ifv->d_u) * (ifv->gamma-1.0);	
			ifv->d_p +=         (- 0.5*ifv->d_rho*ifv->V*ifv->V - ifv->RHO*ifv->V*ifv->d_v) * (ifv->gamma-1.0);

			ifv->U_rho += cv->gradx_rho[k]*delta_x + cv->grady_rho[k]*delta_y;
			ifv->U_e   += cv->gradx_e[k]  *delta_x + cv->grady_e[k]  *delta_y;
			ifv->U_u   += cv->gradx_u[k]  *delta_x + cv->grady_u[k]  *delta_y;
			ifv->U_v   += cv->gradx_v[k]  *delta_x + cv->grady_v[k]  *delta_y;
#ifdef MULTIFLUID_BASICS
			ifv->d_phi = (ifv->d_phi - ifv->PHI*ifv->d_rho)/ifv->RHO;
			if ((_Bool)config[60])
				ifv->d_gamma = (ifv->d_gamma - ifv->gamma*ifv->d_rho)/ifv->RHO;
			ifv->U_phi += cv->gradx_phi[k]*delta_x + cv->grady_phi[k]*delta_y;
#endif
			if(cons2prim(ifv) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return 0;
				}
		}
	return 1;
}


static int order2_i_f_var0(struct i_f_var * ifv)
{		
	ifv->d_rho = 0.0;
	ifv->d_e   = 0.0;
	ifv->d_u   = 0.0;
	ifv->d_v = 0.0;

	ifv->t_rho = 0.0;
	ifv->t_e   = 0.0;
	ifv->t_u   = 0.0;
	ifv->t_v = 0.0;

#ifdef MULTIFLUID_BASICS
	ifv->d_z_a = 0.0;
	ifv->t_z_a = 0.0;
	ifv->d_phi = 0.0;
	ifv->t_phi = 0.0;
#endif

	if(cons2prim(ifv) == 0)
		{
			fprintf(stderr, "Error happens on primitive variable!\n");
			return 0;
		}
		
	return 1;
}

	  
int interface_var_init(const struct cell_var * cv, const struct mesh_var * mv,
					   struct i_f_var * ifv, struct i_f_var * ifv_R,
					   const int k, const int j, const int i, const double gauss)
{
	const int order = (int)config[9];
	int **cc = cv->cell_cell;
	int **cp = mv->cell_pt;

	int p_p, p_n;
	if(j == cp[k][0]-1) 
		{
			p_p=cp[k][1];
			p_n=cp[k][j+1];
		}				  
	else
		{
			p_p=cp[k][j+2];
			p_n=cp[k][j+1];
		}

	ifv->n_x = cv->n_x[k][j];
	ifv->n_y = cv->n_y[k][j];
	ifv->length = sqrt((mv->X[p_p] - mv->X[p_n])*(mv->X[p_p] - mv->X[p_n]) + (mv->Y[p_p] - mv->Y[p_n])*(mv->Y[p_p] - mv->Y[p_n]));

	cons_qty_copy_cv2ifv(ifv, cv, k);
   
	if (order == 2)
		{
			ifv->delta_x = 0.5*(mv->X[p_p]*(1.0+gauss) + mv->X[p_n]*(1.0-gauss)) - cv->X_c[k];
			ifv->delta_y = 0.5*(mv->Y[p_p]*(1.0+gauss) + mv->Y[p_n]*(1.0-gauss)) - cv->Y_c[k];
			if(order2_i_f_var_init(cv, ifv, k) == 0)			
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return 0;
				}
		}
		
	ifv_R->n_x = ifv->n_x;
	ifv_R->n_y = ifv->n_y;
	ifv_R->length = ifv->length;
	
	int cR; //cell_right	
	if (cc[k][j] >= 0)
		{
			cR = cc[k][j];
			cons_qty_copy_cv2ifv(ifv_R, cv, cR);			

			if (order == 2)
				{
					ifv_R->delta_x = 0.5*(mv->X[p_p]*(1.0+gauss) + mv->X[p_n]*(1.0-gauss)) - cv->X_c[cR];
					ifv_R->delta_y = 0.5*(mv->Y[p_p]*(1.0+gauss) + mv->Y[p_n]*(1.0-gauss)) - cv->Y_c[cR];
					if(order2_i_f_var_init(cv, ifv_R, cR) == 0)
						{
							fprintf(stderr, "Error happens on primitive variable!\n");
							return 0;
						}
				}
		}
	else if (cc[k][j] == -1)//initial boundary condition.
		{
			if (i > 1)
				return -1;
			cons_qty_copy_cv2ifv(ifv_R, cv, k);

			if (order == 2)
				if(order2_i_f_var0(ifv_R) == 0)
					{
						fprintf(stderr, "Error happens on primitive variable!\n");
						return 0;
					}
		}
	else if (cc[k][j] == -3)//prescribed boundary condition.
		{
			cons_qty_copy_cv2ifv(ifv_R, cv, k);			

			if (order == 2)
				if(order2_i_f_var0(ifv_R) == 0)
					{
						fprintf(stderr, "Error happens on primitive variable!\n");
						return 0;
					}
		}		
	else if (cc[k][j] != -2&&cc[k][j] != -4)
		{
			printf("No suitable boundary!cc = %d!\n",cc[k][j]);
			return 0;
		}

	if (order == 1)
		{
			if(cons2prim(ifv) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return 0;
				}
			if (cc[k][j] != -2&&cc[k][j] != -4)
				if(cons2prim(ifv_R) == 0)
					{
						fprintf(stderr, "Error happens on primitive variable!\n");
						return 0;
					}
		}

	double u_R, v_R;
	if (cc[k][j] == -2)//reflecting boundary condition.
		{
			*ifv_R = *ifv;
			u_R =  ifv_R->U*ifv_R->n_x + ifv_R->V*ifv_R->n_y;
			v_R = -ifv_R->U*ifv_R->n_y + ifv_R->V*ifv_R->n_x;
			u_R = -u_R;
			ifv_R->U = u_R*ifv_R->n_x - v_R*ifv_R->n_y;
			ifv_R->V = u_R*ifv_R->n_y + v_R*ifv_R->n_x;
		}
	else if (cc[k][j] == -4)//symmetry boundary condition.
		*ifv_R = *ifv;

	return 1;
}


double tau_calc(const struct cell_var * cv, const struct mesh_var * mv)
{
	const double CFL = config[7];
	if (CFL < 0.0)
		return -CFL;
	const int num_cell = (int)config[3];
	int ** cp = mv->cell_pt;
	
	double tau = config[1];
	struct i_f_var ifv, ifv_R;
	double cum, lambda_max;
	int ivi;
	
	double qn, qn_R;
	double c, c_R;	
	
	for(int k = 0; k < num_cell; ++k)
		{
			cum = 0.0;
			
			for(int j = 0; j < cp[k][0]; ++j)
				{
					ivi = interface_var_init(cv, mv, &ifv, &ifv_R, k, j, 0, 0.0);
					if (ivi < 0)
						;
					else if(ivi == 0)
						return -1.0;
					else
						{
							qn = ifv.U*ifv.n_x + ifv.V*ifv.n_y; 
							qn_R = ifv_R.U*ifv_R.n_x + ifv_R.V*ifv_R.n_y;
							c = sqrt(ifv.gamma * ifv.P / ifv.RHO);
							c_R = sqrt(ifv_R.gamma * ifv_R.P / ifv_R.RHO);
							lambda_max = fmax(c+fabs(qn), c_R+fabs(qn_R));
							cum += 0.5*lambda_max * ifv.length;
						}
				}
			tau = fmin(tau, cv->vol[k]/cum * CFL);
		} //To decide tau.
	return tau;
}
