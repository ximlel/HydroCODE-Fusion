#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/file_io.h"
#include "../include/inter_process_unstruct.h"
#include "../include/flux_calc.h"


/**
 * @brief  This function use various finite volume schemes to solve (augmented) Euler equations for single-/two-component fluid
 *         motion on unstructured grids in Eulerian coordinate.
 * @param[in,out] FV:     Structure of initial fluid variable data array pointer.
 * @param[in] mv:         Structure of meshing variable data.
 * @param[in] scheme:     Scheme name.
 * @param[in] problem:    Name of the numerical results for the test problem.
 * @param[in] N_plot:     Number of time steps for plotting.
 * @param[in] time_plot:  Array of the plotting time recording.
 */
void finite_volume_scheme_unstruct(struct flu_var * FV, const struct mesh_var * mv, const char * scheme,
				   const char * problem, int * N_plot, double time_plot[])
{
	int const num_cell = (int)config[3];  // Total grid cell number
	int const N        = (int)config[5];  // the maximum number of time steps
	int const el       = (int)config[8];  // coordinate framework
	int const order    = (int)config[9];  // order of the scheme
	double const t_all =      config[1];  // the total time
	double const eps   =      config[4];  // the largest value could be seen as zero
	double       tau   =      config[16]; // the length of the time step

	clock_t start_clock;
	double cpu_time = 0.0;

	int ** cp = mv->cell_pt;

	struct cell_var cv;
	cell_mem_init_free(&cv, mv, FV, 1); // Initialize memory

	cons_qty_init(&cv, FV);

	vol_comp(&cv, mv);

	cell_rel(&cv, mv);

	if (order > 1)
		cell_centroid(&cv, mv);

	printf("Unstructured grid has been constructed.\n");

	struct i_f_var ifv, ifv_R;
	double time_c = 0.0;
	_Bool stop_t = false;
	int i, ivi, RK = 0, N_count = 0;
	for(i = 1; i <= N; ++i)
		{
			start_clock = clock();
			if (time_c >= time_plot[N_count] && N_count < (*N_plot-1))
				{
					file_write_2D_BLOCK_TEC(*FV, *mv, problem, time_plot[N_count]);
					N_count++;
				}

			fluid_var_update(FV, &cv);

			if (order > 1)
				{
					if (el != 0 && i > 1) // @todo ALE grid movement
						cell_centroid(&cv, mv);
					if (mv->bc != NULL)
						mv->bc(&cv, mv, FV, time_c);
					if (!(int)config[31])
						slope_limiter_prim(&cv, mv, FV);
				}
			if (mv->bc != NULL)
				mv->bc(&cv, mv, FV, time_c);

			if(isfinite(t_all) || !isfinite(config[16]) || config[16] <= 0.0 || !RK)
			    {
				tau = tau_calc(&cv, mv);
				if(tau < eps)
				    {
					printf("\nThe length of the time step is so small on [%d, %g, %g] (t_n, time_c, tau)\n", i, time_c, tau);
					stop_t = true;
				    }
				else if((time_c + tau) > (t_all - eps))
				    {
					printf("\nThe time is enough at step %d.\n",i);				
					tau = t_all - time_c;
				    }
				else if(!isfinite(tau))
				    {
					printf("NAN or INFinite error on [%d, %g, %g] (t_n, time_c, tau) - CFL\n", i, time_c, tau); 
					goto return_NULL;
				    }
			    }

			for(int k = 0; k < num_cell; k++)
				{
					for(int j = 0; j < cp[k][0]; j++)
						{
							ivi = interface_var_init(&cv, mv, &ifv, &ifv_R, k, j, i, 0.0);
							// ivi = interface_var_init(&cv, mv, &ifv, &ifv_R, k, j, i, 1.0/sqrt(3));
							if(ivi == 0)
								stop_t = true;
							else if (ivi == 1)
								{
									if (order == 1)
										{
											if (strcmp(scheme,"Roe") == 0)
												Roe_flux(&ifv, &ifv_R);
											else if (strcmp(scheme,"HLL") == 0)
												HLL_flux(&ifv, &ifv_R);
											else if(strcmp(scheme,"Riemann_exact") == 0)
												Riemann_exact_flux(&ifv, &ifv_R);
											else
												{
													printf("No Riemann solver!\n");
													exit(4);
												}
										}
									else if (order == 2)
										{
											if(strcmp(scheme,"GRP_2D") == 0)
												GRP_2D_flux(&ifv, &ifv_R, tau);
											else
												{
													printf("No Riemann solver!\n");
													exit(4);
												}
										}
								}
							if (ivi != -1)
								flux_copy_ifv2cv(&ifv, &cv, k ,j);
/*
							ivi = interface_var_init(&cv, mv, &ifv, &ifv_R, k, j, i, 0.0);
							//ivi = interface_var_init(&cv, mv, &ifv, &ifv_R, k, j, i, -1.0/sqrt(3));
							if(ivi == 0)
								stop_t = true;
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

			// cons_qty_update(&cv, mv, *FV, tau);
			if (cons_qty_update_corr_ave_P(&cv, mv, FV, tau, RK) == 0)
			    stop_t = true;

			if((_Bool)config[53])
			    RK = RK ? 0 : 1;
			if(!(_Bool)config[53] || RK == 1)
			    time_c += tau;
			if(isfinite(t_all))
			    DispPro(time_c*100.0/t_all, i);
			else
			    DispPro(i*100.0/N, i);
			if(stop_t || time_c > (t_all - eps) || !isfinite(time_c))
			    break;

			cpu_time += (clock() - start_clock) / (double)CLOCKS_PER_SEC;
		}
	printf("\nTime is up at time step %d.\n", i);
	printf("\nThe cost of CPU time for the finite volume scheme on unstructured grids in Eulerian coordinate is %g seconds.\n", cpu_time);

return_NULL:
	config[5] = (double)i;
  *N_plot = N_count+1;
  if(isfinite(time_c))
      time_plot[N_count] = time_c;
  else if(isfinite(t_all))
      time_plot[N_count] = t_all;
  else if(isfinite(tau))
      time_plot[N_count] = i*tau;

	fluid_var_update(FV, &cv);
	cell_mem_init_free(&cv, mv, FV, 0);
}
