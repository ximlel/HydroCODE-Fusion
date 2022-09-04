#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "../include/var_struc.h"
#include "../include/tools.h"


static inline double mu_BJ(double x)
{
	return (x<1.0?x:1.0);
}
static inline double mu_Ven(double x)
{
	return ((x*x+2.0*x)/(x*x+x+2.0));
}

static void lsq_limiter(const struct cell_var * cv, const struct mesh_var * mv, 
						double * grad_W_x, double * grad_W_y, double * W)
{
	const double eps = config[4];
	const int num_cell = (int)config[3];
	const int lim = (int)config[40]; //limiter
	double (*mu[])(double) = { mu_BJ, mu_Ven };
	
	int **cp = mv->cell_pt;
	int **cc = cv->cell_cell;
	const double *X_c = cv->X_c;
	const double *Y_c = cv->Y_c;	
	const double *X = mv->X;
	const double *Y = mv->Y;

	int cell_R;
	double tmp_x, tmp_y;
	double M_c[2][2];
	double W_c_min, W_c_max, W_c_x_p;
	double fai_W;
	int p_p,p_n;
	
	for(int k = 0; k < num_cell; ++k)
		{
			M_c[0][0] = 0.0;  M_c[0][1] = 0.0;
			M_c[1][0] = 0.0;  M_c[1][1] = 0.0;
			grad_W_x[k] = 0.0;
			grad_W_y[k] = 0.0;

			for(int j = 0; j < cp[k][0]; j++)
				{
					if (cc[k][j] >= 0)
						cell_R = cc[k][j];
					else if (cc[k][j] == -1 || cc[k][j] == -2 || cc[k][j] == -3 || cc[k][j] == -4)
						continue;
					else
						{
							fprintf(stderr, "No suitable boundary!\n");
							exit(2);
						}
							
					M_c[0][0] += (X_c[cell_R] - X_c[k]) * (X_c[cell_R] - X_c[k]);
					M_c[0][1] += (X_c[cell_R] - X_c[k]) * (Y_c[cell_R] - Y_c[k]);
					M_c[1][0] += (Y_c[cell_R] - Y_c[k]) * (X_c[cell_R] - X_c[k]);
					M_c[1][1] += (Y_c[cell_R] - Y_c[k]) * (Y_c[cell_R] - Y_c[k]);
					grad_W_x[k] += (W[cell_R] - W[k]) * (X_c[cell_R] - X_c[k]);
					grad_W_y[k] += (W[cell_R] - W[k]) * (Y_c[cell_R] - Y_c[k]);
				}		
			//inverse
			if(rinv(M_c[0], 2) == 0)
				exit(3);

			tmp_x = M_c[0][0] * grad_W_x[k] + M_c[0][1] * grad_W_y[k];
			tmp_y = M_c[1][0] * grad_W_x[k] + M_c[1][1] * grad_W_y[k];
			grad_W_x[k] = tmp_x;
			grad_W_y[k] = tmp_y;
		}

	for(int k = 0; k < num_cell; ++k)
		{
			W_c_min = W[k];
			W_c_max = W[k];
			for(int j = 0; j < cp[k][0]; ++j)
				{
					if (cc[k][j] >= 0)
						cell_R = cc[k][j];
					else if (cc[k][j] == -1 || cc[k][j] == -2 || cc[k][j] == -3 || cc[k][j] == -4)
						continue;
					else
						{
							printf("No suitable boundary!\n");
							exit(2);
						}
							
					if(W[cell_R] < W_c_min)
						W_c_min = W[cell_R];
					else if(W[cell_R] > W_c_max)
						W_c_max = W[cell_R];
				}
			fai_W = 1.0;
			for(int j = 0; j < cp[k][0]; ++j)
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
					//					W_c_x_p = W[k] + grad_W_x[k] * (X[cp[k][j+1]] - X_c[k]) + grad_W_y[k] * (Y[cp[k][j+1]] - Y_c[k]);
					W_c_x_p = W[k] + grad_W_x[k] * (0.5*(X[p_p]+X[p_n]) - X_c[k]) + grad_W_y[k] * (0.5*(Y[p_p]+Y[p_n]) - Y_c[k]);	
					if (fabs(W_c_x_p - W[k]) < eps)
						;	
					else if((W_c_x_p - W[k]) > 0.0)
						fai_W = fmin(fai_W, mu[lim]((W_c_max - W[k])/(W_c_x_p - W[k])));
					else
						fai_W = fmin(fai_W, mu[lim]((W_c_min - W[k])/(W_c_x_p - W[k])));
				}
			grad_W_x[k] = grad_W_x[k] * fai_W;
			grad_W_y[k] = grad_W_y[k] * fai_W;
		}
}

static void minmod_limiter(const struct cell_var * cv, const struct mesh_var * mv, 
						   double * grad_W, double * W)
{
	// const double eps = config[4];
	const int num_cell = (int)config[3];
	double alpha = config[41];
	int **cc = cv->cell_cell;

	int cell_R;
	double grad_W_tmp;
	for(int k = 0; k < num_cell; k++)
		{
			for(int j = 0; j < 2; j++)
				{
					if (cc[k][j] >= 0)
						{
							cell_R = cc[k][j];
							grad_W_tmp = alpha*(W[cell_R] - W[k]) / (cv->X_c[cell_R] - cv->X_c[k]);
						}
					else if (cc[k][j] == -1 || cc[k][j] == -3 || cc[k][j] == -4)
						grad_W_tmp = 0.0;
					else if (cc[k][j] == -2)
						continue;
					else
						{
							fprintf(stderr, "No suitable boundary!\n");
							exit(2);
						}
					/*
					  if (grad_W_tmp * grad_W[k] < eps*eps)
					  {
					  grad_W[k] = 0.0;
					  break;
					  }
					  else if (grad_W[k] > 0.0)
					  grad_W[k] = fmin(grad_W[k], grad_W_tmp);
					  else
					  grad_W[k] = fmax(grad_W[k], grad_W_tmp);
					*/
					if (grad_W_tmp > 0.0 && grad_W[k] > 0.0)
						grad_W[k] = fmin(grad_W[k], grad_W_tmp);
					else if (grad_W_tmp < 0.0 && grad_W[k] < 0.0)
						grad_W[k] = fmax(grad_W[k], grad_W_tmp);
					else
						grad_W[k] = 0.0;
				}
		}
}

static void minmod_limiter_2D(const struct cell_var * cv, const struct mesh_var * mv, 
							  double * gradx_W, double * grady_W, const double * W,
							  const int isUorV)//isUorV(U:1,V:-1,NO:0)
{
	// const double eps = config[4];
	const int num_cell = (int)config[3];
	double alpha = config[41];
	int **cc = cv->cell_cell;
	int **cp = mv->cell_pt;

	int cell_R, p_p, p_n;
	double grad_W_tmp;
	for(int k = 0; k < num_cell; k++)
		{
			for(int j = 0; j < 4 ;j++)
				{
					if (j == 1 || j == 3)
						{
							if (cc[k][j] >= 0)
								{
									cell_R = cc[k][j];
									grad_W_tmp = alpha*(W[cell_R] - W[k]) / (cv->X_c[cell_R] - cv->X_c[k]);
								}
							else if (cc[k][j] == -1 || cc[k][j] == -3 || cc[k][j] == -4)
								grad_W_tmp = 0.0;
							else if (cc[k][j] == -2)
								{
									if (isUorV == 0 || isUorV == -1)
										grad_W_tmp = 0.0;
									else if (isUorV == 1)
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
											grad_W_tmp = alpha*W[k] / (cv->X_c[k] - 0.5*(mv->X[p_p]+mv->X[p_n]));
										}
								}
							else
								{
									fprintf(stderr, "No suitable boundary!cc = %d,%d,%d\n",cc[k][j],k,j);
									//exit(2);
								}
							if (isinf(gradx_W[k]))
								gradx_W[k] = grad_W_tmp;
							else
								{
									/*
									  if (grad_W_tmp * gradx_W[k] < eps*eps)
									  gradx_W[k] = 0.0;
									  else if (gradx_W[k] > 0.0)
									  gradx_W[k] = fmin(gradx_W[k], grad_W_tmp);
									  else
									  gradx_W[k] = fmax(gradx_W[k], grad_W_tmp);
									*/
									if (grad_W_tmp > 0.0 && gradx_W[k] > 0.0)
										gradx_W[k] = fmin(gradx_W[k], grad_W_tmp);
									else if (grad_W_tmp < 0.0 && gradx_W[k] < 0.0)
										gradx_W[k] = fmax(gradx_W[k], grad_W_tmp);
									else
										gradx_W[k] = 0.0;
								}
						}
					else
						{
							if (cc[k][j] >= 0)
								{
									cell_R = cc[k][j];
									grad_W_tmp = alpha*(W[cell_R] - W[k]) / (cv->Y_c[cell_R] - cv->Y_c[k]);
								}
							else if (cc[k][j] == -1 || cc[k][j] == -3 || cc[k][j] == -4)
								grad_W_tmp = 0.0;
							else if (cc[k][j] == -2)
								{
									if (isUorV == 0 || isUorV == 1)
										grad_W_tmp = 0.0;
									else if (isUorV == -1)
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
											grad_W_tmp = alpha*W[k] / (cv->Y_c[k] - 0.5*(mv->Y[p_p]+mv->Y[p_n]));
										}
								}
							else
								{
									fprintf(stderr, "No suitable boundary!cc = %d,%d,%d\n",cc[k][j],k,j);
									//exit(2);
								}
							if (isinf(grady_W[k]))
								grady_W[k] = grad_W_tmp;
							else
								{
									/*
									  if (grad_W_tmp * grady_W[k] < eps*eps)
									  grady_W[k] = 0.0;
									  else if (grady_W[k] > 0.0)
									  grady_W[k] = fmin(grady_W[k], grad_W_tmp);
									  else
									  grady_W[k] = fmax(grady_W[k], grad_W_tmp);
									*/
									if (grad_W_tmp > 0.0 && grady_W[k] > 0.0)
										grady_W[k] = fmin(grady_W[k], grad_W_tmp);
									else if (grad_W_tmp < 0.0 && grady_W[k] < 0.0)
										grady_W[k] = fmax(grady_W[k], grad_W_tmp);
									else
										grady_W[k] = 0.0;
								}
						}
				}
		}
}


void slope_limiter(const struct cell_var * cv,const struct mesh_var * mv, const struct flu_var * FV)
{
	const int num_cell = (int)config[3];
	// int **cc = cv->cell_cell;

	if ((int)config[31] == 0)
	    {
		/*
		  lsq_limiter(cv, mv, cv->gradx_rho, cv->grady_rho, FV->RHO);
		  lsq_limiter(cv, mv, cv->gradx_e, cv->grady_e, FV->P);
		  lsq_limiter(cv, mv, cv->gradx_u, cv->grady_u, FV->U);
		  lsq_limiter(cv, mv, cv->gradx_v, cv->grady_v, FV->V);
#ifdef MULTIFLUID_BASICS
		  lsq_limiter(cv, mv, cv->gradx_phi, cv->grady_phi, FV->PHI);
		  lsq_limiter(cv, mv, cv->gradx_z_a, cv->grady_z_a, FV->Z_a);
#endif
		*/
		for(int k = 0; k < num_cell; k++)
		    {
			cv->gradx_rho[k] = (cv->RHO_p[k][1] - cv->RHO_p[k][3])/config[10];
			cv->grady_rho[k] = (cv->RHO_p[k][2] - cv->RHO_p[k][0])/config[11];
			cv->gradx_u[k]   = (cv->U_p[k][1]   - cv->U_p[k][3])  /config[10];
			cv->grady_u[k]   = (cv->U_p[k][2]   - cv->U_p[k][0])  /config[11];
			cv->gradx_v[k]   = (cv->V_p[k][1]   - cv->V_p[k][3])  /config[10];
			cv->grady_v[k]   = (cv->V_p[k][2]   - cv->V_p[k][0])  /config[11];
			cv->gradx_e[k]   = (cv->P_p[k][1]   - cv->P_p[k][3])  /config[10];
			cv->grady_e[k]   = (cv->P_p[k][2]   - cv->P_p[k][0])  /config[11];
#ifdef MULTIFLUID_BASICS
			cv->gradx_phi[k] = (cv->PHI_p[k][1] - cv->PHI_p[k][3])/config[10];
			cv->grady_phi[k] = (cv->PHI_p[k][2] - cv->PHI_p[k][0])/config[11];
			cv->gradx_z_a[k] = (cv->Z_a_p[k][1] - cv->Z_a_p[k][3])/config[10];
			cv->grady_z_a[k] = (cv->Z_a_p[k][2] - cv->Z_a_p[k][0])/config[11];	
#endif
			/*
			  if (cc[k][0] < 0 || cc[k][1] < 0 || cc[k][2] < 0 || cc[k][3] < 0)
			  {
			  cv->gradx_u[k]   = 1.0/0.0;
			  cv->grady_u[k]   = 1.0/0.0;
			  cv->gradx_v[k]   = 1.0/0.0;
			  cv->grady_v[k]   = 1.0/0.0;
			  cv->gradx_e[k]   = 1.0/0.0;
			  cv->grady_e[k]   = 1.0/0.0;
#ifdef MULTIFLUID_BASICS
			  cv->gradx_phi[k] = 1.0/0.0;
			  cv->grady_phi[k] = 1.0/0.0;
			  cv->gradx_z_a[k] = 1.0/0.0;
			  cv->grady_z_a[k] = 1.0/0.0;	
#endif
			  }
			  else
			  {
			  cv->gradx_rho[k] = (FV->RHO[cc[k][1]]-FV->RHO[cc[k][3]]) /config[10]/2.0;
			  cv->grady_rho[k] = (FV->RHO[cc[k][2]]-FV->RHO[cc[k][0]]) /config[11]/2.0;
			  cv->gradx_u[k]   = (FV->U[cc[k][1]] - FV->U[cc[k][3]]) /config[10]/2.0;
			  cv->grady_u[k]   = (FV->U[cc[k][2]] - FV->U[cc[k][0]]) /config[11]/2.0;
			  cv->gradx_v[k]   = (FV->V[cc[k][1]] - FV->V[cc[k][3]]) /config[10]/2.0;
			  cv->grady_v[k]   = (FV->V[cc[k][2]] - FV->V[cc[k][0]]) /config[11]/2.0;
			  cv->gradx_e[k]   = (FV->P[cc[k][1]] - FV->P[cc[k][3]]) /config[10]/2.0;
			  cv->grady_e[k]   = (FV->P[cc[k][2]] - FV->P[cc[k][0]]) /config[11]/2.0;
#ifdef MULTIFLUID_BASICS
			  cv->gradx_phi[k] = (FV->PHI[cc[k][1]] - FV->PHI[cc[k][3]])/config[10];
			  cv->grady_phi[k] = (FV->PHI[cc[k][2]] - FV->PHI[cc[k][0]])/config[11];
			  cv->gradx_z_a[k] = (FV->Z_a[cc[k][1]] - FV->Z_a[cc[k][3]])/config[10];
			  cv->grady_z_a[k] = (FV->Z_a[cc[k][2]] - FV->Z_a[cc[k][0]])/config[11];
#endif
			  }

			  if (fabs(cv->gradx_rho[k]-(cv->RHO_p[k][1] - cv->RHO_p[k][3])/config[10])>0.00001)
			  {
			  printf("X:%.10lf,%.10lf,%d\n",cv->gradx_rho[k],(cv->RHO_p[k][1] - cv->RHO_p[k][3])/config[10],k);
			  printf("Y:%.10lf,%.10lf,%d\n",cv->grady_rho[k],(cv->RHO_p[k][2] - cv->RHO_p[k][0])/config[11],k);
			  }
			*/
		    }
		/*
		  minmod_limiter_2D(cv, mv, cv->gradx_rho, cv->grady_rho, FV->RHO, 0);
		  minmod_limiter_2D(cv, mv, cv->gradx_e,   cv->grady_e,   FV->P, 0);
		  minmod_limiter_2D(cv, mv, cv->gradx_u,   cv->grady_u,   FV->U, 1);
		  minmod_limiter_2D(cv, mv, cv->gradx_v,   cv->grady_v,   FV->V,-1);
#ifdef MULTIFLUID_BASICS
		  minmod_limiter_2D(cv, mv, cv->gradx_phi, cv->grady_phi, FV->PHI, 0);
		  minmod_limiter_2D(cv, mv, cv->gradx_z_a, cv->grady_z_a, FV->Z_a, 0);
#endif
		  }
		*/
	    }
}

/*
void slope_limiter2(const struct cell_var * cv,const struct mesh_var * mv, const struct flu_var * FV)
{
	const int num_cell = (int)config[3];
	int **cp = mv->cell_pt;
	int **cc = cv->cell_cell;
	const double *X = mv->X;
	const double *Y = mv->Y;
	int p_p,p_n;
	double Xu,Xd,Xl,Xr, Yu,Yd,Yl,Yr;

	if ((int)config[31] == 0)
	    {
		lsq_limiter(cv, mv, cv->gradx_rho, cv->grady_rho, FV->RHO);
		lsq_limiter(cv, mv, cv->gradx_e, cv->grady_e, FV->P);
		lsq_limiter(cv, mv, cv->gradx_u, cv->grady_u, FV->U);
		lsq_limiter(cv, mv, cv->gradx_v, cv->grady_v, FV->V);
#ifdef MULTIFLUID_BASICS
		lsq_limiter(cv, mv, cv->gradx_phi, cv->grady_phi, FV->PHI);
		lsq_limiter(cv, mv, cv->gradx_z_a, cv->grady_z_a, FV->Z_a);
#endif

		for(int k = 0; k < num_cell; k++)
		    {
			for(int j = 0; j < cp[k][0]; ++j)
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
				(0.5*(X[p_p]+X[p_n]) - X_c[k]) (0.5*(Y[p_p]+Y[p_n]);
				}
				cv->gradx_rho[k] = (cv->RHO_p[k][1] - cv->RHO_p[k][3])/config[10];
				cv->grady_rho[k] = (cv->RHO_p[k][2] - cv->RHO_p[k][0])/config[11];
				cv->gradx_u[k]   = (cv->U_p[k][1]   - cv->U_p[k][3])  /config[10];
				cv->grady_u[k]   = (cv->U_p[k][2]   - cv->U_p[k][0])  /config[11];
				cv->gradx_v[k]   = (cv->V_p[k][1]   - cv->V_p[k][3])  /config[10];
				cv->grady_v[k]   = (cv->V_p[k][2]   - cv->V_p[k][0])  /config[11];
				cv->gradx_e[k]   = (cv->P_p[k][1]   - cv->P_p[k][3])  /config[10];
				cv->grady_e[k]   = (cv->P_p[k][2]   - cv->P_p[k][0])  /config[11];
#ifdef MULTIFLUID_BASICS
					cv->gradx_phi[k] = (cv->PHI_p[k][1] - cv->PHI_p[k][3])/config[10];
					cv->grady_phi[k] = (cv->PHI_p[k][2] - cv->PHI_p[k][0])/config[11];
					cv->gradx_z_a[k] = (cv->Z_a_p[k][1] - cv->Z_a_p[k][3])/config[10];
					cv->grady_z_a[k] = (cv->Z_a_p[k][2] - cv->Z_a_p[k][0])/config[11];	
#endif
			    }
			minmod_limiter_2D(cv, mv, cv->gradx_rho, cv->grady_rho, FV->RHO, 0);
			minmod_limiter_2D(cv, mv, cv->gradx_e,   cv->grady_e,   FV->P, 0);
			minmod_limiter_2D(cv, mv, cv->gradx_u,   cv->grady_u,   FV->U, 1);
			minmod_limiter_2D(cv, mv, cv->gradx_v,   cv->grady_v,   FV->V,-1);
#ifdef MULTIFLUID_BASICS
			minmod_limiter_2D(cv, mv, cv->gradx_phi, cv->grady_phi, FV->PHI, 0);
			minmod_limiter_2D(cv, mv, cv->gradx_z_a, cv->grady_z_a, FV->Z_a, 0);
#endif
			    }
		    }
}
*/
