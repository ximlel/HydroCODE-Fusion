#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"


/**
 * @brief Initialize conservative variables on grid cells in struct 'cv' according to primitive variables in struct 'FV'.
 * @param[in] cv:     Structure of grid variable data in computational grid cells.
 * @param[in] FV:     Structure of initial fluid variable data array pointer.
 */
void cons_qty_init(const struct cell_var * cv, const struct flu_var * FV)
{
	const int num_cell = (int)config[3];
	for(int k = 0; k < num_cell; k++)
		{
			cv->U_rho[k] = FV->RHO[k];
			cv->U_u[k]   = FV->RHO[k] * FV->U[k];	
			cv->U_e[k]   = FV->RHO[k] * FV->U[k] * FV->U[k] * 0.5;
			cv->U_v[k]   = FV->RHO[k] * FV->V[k];
			cv->U_e[k]  += FV->RHO[k] * FV->V[k] * FV->V[k] * 0.5;
#ifdef MULTIFLUID_BASICS
			cv->U_phi[k]   = FV->RHO[k] * FV->PHI[k];			
			cv->U_gamma[k] = FV->RHO[k] * FV->gamma[k];
			cv->U_e_a[k]   = FV->Z_a[k] * FV->P[k]/(config[6]-1.0) + 0.5*cv->U_phi[k]*(FV->U[k]*FV->U[k]+FV->V[k]*FV->V[k]);
			cv->U_e[k]    += FV->P[k]/(FV->gamma[k]-1.0);
#else
			cv->U_e[k]    += FV->P[k]/(config[6]-1.0);
#endif
		}
}


int cons2prim(struct i_f_var * ifv)
{
	const double eps = config[4];
	
	ifv->RHO   = ifv->U_rho;
	ifv->U     = ifv->U_u/ifv->U_rho;
	ifv->V  = ifv->U_v/ifv->U_rho;
	if (isnan(ifv->V) || isinf(ifv->V))									
		return 0;
#ifdef MULTIFLUID_BASICS
	double phi_e_a, phi_e_b;
	ifv->PHI = ifv->U_phi/ifv->U_rho;
	phi_e_a  = ifv->U_e_a-0.5*ifv->U_phi*(ifv->U*ifv->U+ifv->V*ifv->V);
	// phi_e_b  = (ifv->U_e-ifv->U_e_a)-0.5*(ifv->U_rho-ifv->U_phi)*(ifv->U*ifv->U+ifv->V*ifv->V);
	phi_e_b  = ifv->U_e-0.5*(ifv->U_rho)*(ifv->U*ifv->U+ifv->V*ifv->V)-phi_e_a;
	ifv->Z_a = phi_e_a*(config[6]-1.0)/(phi_e_a*(config[6]-1.0) + phi_e_b*(config[106]-1.0));
	if (isnan(ifv->Z_a)||isnan(ifv->PHI)||ifv->PHI<(-0.01)||ifv->PHI>(1.0+0.01))
	    return 0;
	else if (ifv->Z_a<(-0.001))
	    {
		printf("Z_a=%.10lf,phi_a=%.10lf\n",ifv->Z_a,ifv->PHI);
		ifv->Z_a = 0.0;
		ifv->U_e_a = 0.5*ifv->U_phi*(ifv->U*ifv->U+ifv->V*ifv->V);
	    }
	else if (ifv->Z_a>(1.0+0.001))
	    {
		printf("Z_a=%.10lf,phi_a=%.10lf\n",ifv->Z_a,ifv->PHI);
		ifv->Z_a = 1.0;
		ifv->U_e_a = ifv->U_e-0.5*(ifv->U_rho-ifv->U_phi)*(ifv->U*ifv->U+ifv->V*ifv->V);
	    }
	else if (ifv->PHI<(-0.001)||ifv->PHI>(1.0+0.001))
	    printf("Z_a=%.10lf,phi_a=%.10lf\n",ifv->Z_a,ifv->PHI);
	// ifv->gamma = (phi_e_a*config[6] + phi_e_b*config[106])/(phi_e_a+phi_e_b);
	ifv->gamma = 1.0/(ifv->Z_a/(config[6]-1.0)+(1.0-ifv->Z_a)/(config[106]-1.0))+1.0;
#endif
	ifv->P  = (ifv->U_e - 0.5*(ifv->U_u*ifv->U_u)/ifv->U_rho) * (ifv->gamma-1.0);	
	ifv->P -= (0.5*(ifv->U_v*ifv->U_v)/ifv->U_rho) * (ifv->gamma-1.0);

	
	if (isnan(ifv->RHO + ifv->U + ifv->P) || isinf(ifv->RHO + ifv->U + ifv->P) || ifv->RHO < eps)
		return 0;
	else if (ifv->P < eps)
		{
			printf("P=%.10lf\n",ifv->P);		
			ifv->U_e = 0.5*(ifv->U_u*ifv->U_u)/ifv->U_rho + eps/(ifv->gamma-1.0);
			ifv->U_e += 0.5*(ifv->U_v*ifv->U_v)/ifv->U_rho;
#ifdef MULTIFLUID_BASICS
			ifv->U_e_a = 0.5*ifv->U_phi*(ifv->U*ifv->U+ifv->V*ifv->V) + ifv->Z_a*eps/(config[6]-1.0);
#endif
		}
	
	return 1;
}


int cons_qty_update(const struct cell_var * cv, const struct mesh_var * mv,
					const struct flu_var  * FV, const double tau)
{
	const int num_cell = (int)config[3];
	int ** cp = mv->cell_pt;
	
	int p_p, p_n;
	double length, flux_v_fix;
	// double gamma, length2;
	int i;
	for(int k = 0; k < num_cell; ++k)
		{
			flux_v_fix = 0.0;				
			//if((_Bool)config[60])
			//	gamma = cv->U_gamma[k]/cv->U_rho[k];
			for(int j = 0; j < cp[k][0]; j++)
				{
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
					length = sqrt((mv->X[p_p] - mv->X[p_n])*(mv->X[p_p]-mv->X[p_n]) + (mv->Y[p_p] - mv->Y[p_n])*(mv->Y[p_p]-mv->Y[p_n]));
					
					cv->U_rho[k] += - tau*cv->F_rho[k][j] * length / cv->vol[k];
					cv->U_e[k]   += - tau*cv->F_e[k][j]   * length / cv->vol[k];	
					cv->U_u[k]   += - tau*cv->F_u[k][j]   * length / cv->vol[k];
					cv->U_v[k] += - tau*cv->F_v[k][j] * length / cv->vol[k];
#ifdef MULTIFLUID_BASICS
					cv->U_phi[k] += - tau*cv->F_phi[k][j] * length / cv->vol[k];
					if((_Bool)config[60])
						cv->U_gamma[k] += - tau*cv->F_gamma[k][j] * length / cv->vol[k];
#endif
					if ((int)config[61] == 1)
						{													
							//							flux_v_fix += tau*length/cv->vol[k]*FV->RHO[k]*cv->RHO_p[k][j]*(cv->U_p[k][j]*cv->n_x[k][j]+cv->V_p[k][j]*cv->n_y[k][j])*((cv->U_p[k][j]-FV->U[k])*(cv->U_p[k][j]-FV->U[k])+(cv->V_p[k][j]-FV->V[k])*(cv->V_p[k][j]-FV->V[k]));
							for (i = j+1; i < cp[k][0]; i++)
								{
									if(i == cp[k][0]-1) 
										{
											p_p=cp[k][1];
											p_n=cp[k][i+1];
										}				  
									else
										{
											p_p=cp[k][i+2];
											p_n=cp[k][i+1];
										}
									// length2 = sqrt((mv->X[p_p] - mv->X[p_n])*(mv->X[p_p]-mv->X[p_n]) + (mv->Y[p_p] - mv->Y[p_n])*(mv->Y[p_p]-mv->Y[p_n]));
									// flux_v_fix -= tau*length/cv->vol[k]*cv->RHO_p[k][j]*(cv->U_p[k][j]*cv->n_x[k][j]+cv->V_p[k][j]*cv->n_y[k][j])*tau*length2/cv->vol[k]*cv->RHO_p[k][i]*(cv->U_p[k][i]*cv->n_x[k][i]+cv->V_p[k][i]*cv->n_y[k][i])*((cv->U_p[k][j]-cv->U_p[k][i])*(cv->U_p[k][j]-cv->U_p[k][i])+(cv->V_p[k][j]-cv->V_p[k][i])*(cv->V_p[k][j]-cv->V_p[k][i]));
								}

							flux_v_fix += tau*length/cv->vol[k]*cv->RHO_p[k][j]*(cv->U_p[k][j]*cv->n_x[k][j]+cv->V_p[k][j]*cv->n_y[k][j])*((cv->U_p[k][j]-FV->U[k])*(cv->U_p[k][j]-FV->U[k])+(cv->V_p[k][j]-FV->V[k])*(cv->V_p[k][j]-FV->V[k]))/2.0;

						}
				}

			/*			
						if(!(isinf(config[106]) || isinf(config[110]) || isinf(config[111])))
						{
						gamma  = cv->U_phi[k]/cv->U_rho[k]*config[6]*config[110] + (1.0-cv->U_phi[k]/cv->U_rho[k])*config[106]*config[111];
						gamma /= cv->U_phi[k]/cv->U_rho[k]*config[110] + (1.0-cv->U_phi[k]/cv->U_rho[k])*config[111];
						cv->U_gamma[k] = gamma*cv->U_rho[k];
						}
						else if((_Bool)config[60])
						cv->U_gamma[k] = gamma*cv->U_rho[k];
			*/
			
			cv->U_e[k] += flux_v_fix;///cv->U_rho[k]/2.0;
		}
	return 0;
}
