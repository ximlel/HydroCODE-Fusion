#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/finite_volume.h"


#include "../include/file_io.h"

void finite_volume_scheme(struct flu_var * FV, const struct mesh_var * mv, const char * scheme, const char * problem)
{
	clock_t start_clock;
	double cpu_time = 0.0;

config[53] = 1.0;//R-K

	const int dim = (int)config[0];
	const int order = (int)config[9];
	const int el = isinf(config[8]) ? 0 : (int)config[8];
	const int num_cell = (int)config[3];
	const double eps = config[4];
	int ** cp = mv->cell_pt;

	struct cell_var cv = cell_mem_init(mv, FV);

	cons_qty_init(&cv, FV);

	vol_comp(&cv, mv);

	cell_rel(&cv, mv);

	if (order > 1)
		cell_centroid(&cv, mv);

	printf("Grid has been constructed.\n");

	double tau; // the length of the time step
	double t_all = 0.0;
	struct i_f_var ifv, ifv_R, ifv_tmp;

	const double delta_plot_t = 0.05;
	double plot_t = 0.0;

	int k, j, ivi, stop_step = 0, stop_t = 0;
	int i;
	for(i=0; i < (int)config[5] && stop_step != 1; ++i)
		{
			/*
			if (t_all >= plot_t)
				{				
					file_write_TEC(*FV, *mv, problem, plot_t, dim);
					plot_t += delta_plot_t;
				}
			*/
			  			
			start_clock = clock();

			fluid_var_update(FV, &cv);

			if (order > 1)
				{
					if (el != 0 && i > 0)
						cell_centroid(&cv, mv);
					if (mv->bc != NULL)
						mv->bc(&cv, *mv, FV, t_all);
					//if (!(dim == 1 && i == 0) && i > 50)
//					if (i==0)
//					   slope_limiter2(&cv, mv, FV);
//					else
						slope_limiter(&cv, mv, FV);
				}

			if (mv->bc != NULL)
				mv->bc(&cv, *mv, FV, t_all);
			if(stop_step != 2)
				tau = tau_calc(&cv, mv);					

			if(tau < eps)
				{
					printf("\nThe length of the time step is so small at step %d, t_all=%lf, tau=%lf.\n",i,t_all,tau);
					stop_t = 1;
				}

			if((t_all + tau) > config[1])
				{
					printf("\nThe time is enough at step %d.\n",i);				
					tau = config[1] - t_all;
					if (isinf(config[53]) || stop_step == 2)
						stop_t = 1;
				} // Time

			if (isinf(config[53]) || stop_step == 2)
				t_all += tau;
			config[16] = tau;

			for(k = 0; k < num_cell; k++)
				{
					if (stop_step == 1)
						break;
					for(j = 0; j < cp[k][0]; j++)
						{
							ivi = interface_var_init(&cv, mv, &ifv, &ifv_R, k, j, i, 0.0);
							//ivi = interface_var_init(&cv, mv, &ifv, &ifv_R, k, j, i, 1.0/sqrt(3));
							if(ivi == 0)
								{
									stop_step = 1;
									break;
								}
							else if (ivi == 1)
								{
									if (dim == 1 && cv.n_x[k][j] < 0.0)
										{
											ifv_tmp = ifv;
											ifv     = ifv_R;
											ifv_R   = ifv_tmp;
										}

									if (order == 1)
										{
											if (strcmp(scheme,"Roe") == 0)
												Roe_scheme(&ifv, &ifv_R);
											else if (strcmp(scheme,"HLL") == 0)
												HLL_scheme(&ifv, &ifv_R);
											else if(strcmp(scheme,"Riemann_exact") == 0)
												Riemann_exact_scheme(&ifv, &ifv_R);
											else
												{
													printf("No Riemann solver!\n");
													exit(2);
												}
										}
									else if (order == 2)
										{
											if(strcmp(scheme,"GRP") == 0)					
												GRP_scheme(&ifv, &ifv_R, tau);
											else if(strcmp(scheme,"GRP_2D") == 0)
												GRP_2D_scheme(&ifv, &ifv_R, tau);
											else if(strcmp(scheme,"Riemann_exact") == 0)
												Riemann_exact_scheme(&ifv, &ifv_R);
											else
												{
													printf("No Riemann solver!\n");
													exit(2);
												}
										}
								}
							if (ivi != -1)
								flux_copy_ifv2cv(&ifv, &cv, k ,j);
/*
							ivi = interface_var_init(&cv, mv, &ifv, &ifv_R, k, j, i, 0.0);
							//ivi = interface_var_init(&cv, mv, &ifv, &ifv_R, k, j, i, -1.0/sqrt(3));
							if(ivi == 0)
								{
									stop_step = 1;
									break;
								}
							else if (order == 2 && ivi == 1)
								{
									if(strcmp(scheme,"GRP") == 0)					
										GRP_scheme(&ifv, &ifv_R, tau);
									else if(strcmp(scheme,"GRP_2D") == 0)
										GRP_2D_scheme(&ifv, &ifv_R, tau);
									else if(strcmp(scheme,"Riemann_exact") == 0)
										Riemann_exact_scheme(&ifv, &ifv_R);
									else
										{
											printf("No Riemann solver!\n");
											exit(2);
										}
								}
							flux_add_ifv2cv(&ifv, &cv, k ,j);
*/
						}
				}
			start_clock = clock();

			if (stop_step != 1)
				{
					//	cons_qty_update(&cv, mv, *FV, tau);
					if (cons_qty_update_corr_ave_P(&cv, mv, FV, tau, stop_step) == 0)
						stop_step = 1;
				}

			if(!isinf(config[53]))
				{
					if (stop_step == 0)
						stop_step = 2;
					else if (stop_step == 2)
						stop_step = 0;
				}


			DispPro(t_all*100.0/config[1], i);

			if (stop_step == 1 || stop_t == 1)
				break;

			cpu_time += (clock() - start_clock) / (double)CLOCKS_PER_SEC;
		}
	fluid_var_update(FV, &cv);

//	for (int j = 0; j < num_cell; j++)
//internal energy
//		FV->RHO[j] = (cv.U_e[j]-0.5*cv.U_u[j]*cv.U_u[j]/cv.U_rho[j]-0.5*cv.U_v[j]*cv.U_v[j]/cv.U_rho[j])/cv.U_rho[j];
//entropy
//		FV->V[j] = FV->P[j]/pow(FV->RHO[j],FV->gamma[j]);


	printf("\nThe cost of CPU time for the Eulerian method is %g seconds.\n", cpu_time);
}
