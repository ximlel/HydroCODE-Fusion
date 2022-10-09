#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"
#include "../include/finite_volume.h"
#include "../include/meshing.h"


double config[N_CONF]; //!< Initial configuration data array.

#define CV_INIT_FV_RESET_MEM(v, N)					\
    do {								\
    	CV.v = (double **)malloc(N * sizeof(double *));			\
	if(CV.v == NULL)						\
	    {								\
		printf("NOT enough memory! %s\n", #v);			\
		retval = 5;						\
		goto return_NULL;					\
	    }								\
	for(k = 0; k < N; ++k)						\
	    {								\
		CV.v[k] = (double *)malloc(Md * sizeof(double));	\
		if(CV.v[k] == NULL)					\
		    {							\
			printf("NOT enough memory! %s[%d]\n", #v, k);	\
			retval = 5;					\
			goto return_NULL;				\
		    }							\
	    }								\
	memmove(CV.v[0]+1, FV0.v, Ncell * sizeof(double));		\
	free(FV0.v);							\
	FV0.v = NULL;							\
    } while(0)

/**
 * @brief This is the main function which constructs the
 *        main structure of the Eulerian hydrocode.
 * @param[in] argc: ARGument Counter.
 * @param[in] argv: ARGument Values.
 *            - argv[1]: Folder name of test example (input path).
 *            - argv[2]: Folder name of numerical results (output path).
 *            - argv[3]: Order of numerical scheme[_scheme name] (= 1[_Riemann_exact] or 2[_GRP]).
 *            - argv[4]: Spatial dimension number for radially symmetric flow.
 *              - M=1: planar flow.
 *              - M=2: cylindrical flow.
 *              - M=3: spherical flow.
 *            - argv[5,6,…]: Configuration supplement config[n]=(double)C (= n=C).
 * @return Program exit status code.
 */
int main(int argc, char *argv[])
{
  int k, j, retval = 0;
  // Initialize configuration data array
  for(k = 1; k < N_CONF; k++)
      config[k] = INFINITY;
  char * scheme = NULL; // Riemann_exact(Godunov), GRP
  arg_preprocess(4, argc, argv, scheme);

  // Set dimension.
  config[0] = (double)1; // Dimension of input data = 1

  // The number of times steps of the fluid data stored for plotting.
  int N;
  double *time_plot;
  /* 
   * We read the initial data files.
   * The function initialize return a point pointing to the position
   * of a block of memory consisting Ncell variables of type double.
   * The Ncell variables are the initial value.
   */
  struct flu_var FV0 = initialize_1D(argv[1], &N, &time_plot); // Structure of initial data array pointer.
  const int Ncell = (int)config[3]; // Number of computing cells in r direction
  const int Md    = Ncell+2;        // max vector dimension
  const int order = (int)config[9];
  double gamma = config[6];
  int M = atoi(argv[4]); // M=1 planar; M=2 cylindrical; M=3 spherical
  if(M != 1 && M != 2 && M != 3)
      {
	  printf("Wrong spatial dimension number!\n");
	  exit(4);
      }

  struct radial_mesh_var smv = radial_mesh_init(argv[1]);
  radial_mesh_update(&smv);

  struct cell_var_stru CV = {NULL}; // Structure of fluid variables in computational cells array pointer.
  double ** R = NULL;
  double * cpu_time = (double *)malloc(N * sizeof(double));
  R = (double **)malloc(N * sizeof(double *));
  if(cpu_time == NULL)
      {
	  printf("NOT enough memory! CPU_time\n");
	  retval = 5;
	  goto return_NULL;
      }
  if(R == NULL)
      {
	  printf("NOT enough memory! R\n");
	  retval = 5;
	  goto return_NULL;
      }
  for(k = 0; k < N; ++k)
  {
    R[k] = (double *)malloc(Md * sizeof(double));
    if(R[k] == NULL)
    {
      printf("NOT enough memory! R[%d]\n", k);
      retval = 5;
      goto return_NULL;
    }
  }
  memmove(R[0], smv.RR, Md * sizeof(double));

  CV_INIT_FV_RESET_MEM(U, N);
  CV_INIT_FV_RESET_MEM(P, N);
  CV_INIT_FV_RESET_MEM(RHO, N);
#ifdef MULTIFLUID_BASICS
  CV_INIT_FV_RESET_MEM(gamma, N);
  FV0.gamma = CV.gamma[0];
  for(k = 1; k < N; ++k)
      {
	  free(CV.gamma[k]);
	  CV.gamma[k] = CV.gamma[0];
      }
#endif
  CV.E = (double **)malloc(N * sizeof(double *));
  if(CV.E == NULL)
      {
	  printf("NOT enough memory! E\n");
	  retval = 5;
	  goto return_NULL;
      }
  for(k = 0; k < N; ++k)
  {
    CV.E[k] = (double *)malloc(Md * sizeof(double));
    if(CV.E[k] == NULL)
    {
      printf("NOT enough memory! E[%d]\n", k);
      retval = 5;
      goto return_NULL;
    }
  }
  for(j = 1; j <= Ncell; ++j)
      {
#ifdef MULTIFLUID_BASICS
	  gamma = CV.gamma[0][j];
#endif
	  CV.E[0][j] = 0.5*CV.U[0][j]*CV.U[0][j] + CV.P[0][j]/(gamma-1.0)/CV.RHO[0][j];
      }

  // Use GRP/Godunov scheme to solve it on Lagrangian coordinate.
  config[8] = (double)1;
  switch(order)
      {
      case 1:
	  config[41] = 0.0; // alpha = 0.0
      case 2:
	  grp_solver_radial_LAG_source(CV, &smv, R, M, argv[2], cpu_time, &N, time_plot);
	  break;
      default:
	  printf("NOT appropriate order of the scheme! The order is %d.\n", order);
	  retval = 4;
	  goto return_NULL;
      }

#ifndef NODATPLOT
  file_1D_write(Ncell+1, N, CV, R, cpu_time, argv[2], time_plot);
#endif
#ifdef HDF5PLOT
  file_1D_write_HDF5(Ncell+1, N, CV, R, cpu_time, argv[2], time_plot);
#endif
#ifndef NOTECPLOT
  FV0.RHO = CV.RHO[N-1];
  FV0.U   = CV.U[N-1];
  FV0.P   = CV.P[N-1];
  file_radial_write_TEC(FV0, smv, argv[2], time_plot[N-1]);
#endif

return_NULL:
  radial_mesh_mem_free(&smv);

  FV0.RHO   = NULL;
  FV0.U     = NULL;
  FV0.P     = NULL;
  for(k = 0; k < N; ++k)
  {
    free(CV.E[k]);
    free(CV.RHO[k]);
    free(CV.U[k]);
    free(CV.P[k]);
    free(R[k]);
    CV.E[k]     = NULL;
    CV.RHO[k]   = NULL;
    CV.U[k]     = NULL;
    CV.P[k]     = NULL;
    R[k]        = NULL;
  }
  free(CV.E);
  free(CV.RHO);
  free(CV.U);
  free(CV.P);
  CV.E     = NULL;
  CV.RHO   = NULL;
  CV.U     = NULL;
  CV.P     = NULL;
#ifdef MULTIFLUID_BASICS
  FV0.gamma = NULL;
  free(CV.gamma[0]);
  for(k = 0; k < N; ++k)
      CV.gamma[k] = NULL;
  free(CV.gamma);
  CV.gamma = NULL;
#endif
  free(R);
  R = NULL;
  free(cpu_time);
  cpu_time = NULL;

  return retval;
}
