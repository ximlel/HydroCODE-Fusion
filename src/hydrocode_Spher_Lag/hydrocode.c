#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"
#include "../include/finite_volume.h"
#include "../include/meshing.h"


/**
 * @brief This is the main function which constructs the
 *        main structure of the Eulerian hydrocode.
 * @param[in] argc: ARGument Counter.
 * @param[in] argv: ARGument Values.
 *            - argv[1]: Folder name of test example (input path).
 *            - argv[2]: Folder name of numerical results (output path).
 *            - argv[3]: Dimension of input data (= 1).
 *            - argv[4]: Order of numerical scheme[_scheme name] (= 1[_Riemann_exact] or 2[_GRP]).
 *            - argv[5]: Spatial dimension number for radially symmetric flow.
 *              - M=1: planar flow.
 *              - M=2: cylindrical flow.
 *              - M=3: spherical flow.
 *            - argv[6,7,…]: Configuration supplement config[n]=(double)C (= n=C).
 * @return Program exit status code.
 */
int main(int argc, char *argv[])
{
  int k, retval = 0;
  // Initialize configuration data array
  for(k = 1; k < N_CONF; k++)
      config[k] = INFINITY;
  char * scheme = NULL; // Riemann_exact(Godunov), GRP
  arg_preprocess(5, argc, argv, scheme);
  if ((int)config[0] != 1)
	{
	    printf("No appropriate dimension was entered!\n");
	    return 4;
	}

    // The number of times steps of the fluid data stored for plotting.
	int N; // (int)(config[5]) + 1;
	double *time_plot;
    /* 
     * We read the initial data files.
     * The function initialize return a point pointing to the position
     * of a block of memory consisting (m+1) variables of type double.
     * The value of first array element of these variables is m.
     * The following m variables are the initial value.
     */
	struct flu_var FV0 = initialize_1D(argv[1], &N, &time_plot); // Structure of initial data array pointer.
	int M = atoi(argv[5]); // m=1 planar; m=2 cylindrical; m=3 spherical
    if(M !=1 && M != 2 && M != 3)
	{
	    printf("Wrong spatial dimension number!\n");
	    exit(4);
	}

    struct spher_mesh_var smv = spher_mesh_init(argv[1]);

  double * cpu_time = (double *)malloc(N * sizeof(double));
  if(cpu_time == NULL)
      {
	  printf("NOT enough memory! CPU_time\n");
	  retval = 5;
	  goto return_NULL;
      }

	grp_solver_spher_LAG_source(&FV0, &smv, M, cpu_time, N, time_plot);

	// file_spher_write_TEC(time)；

return_NULL:
    spher_mesh_mem_free(&smv);
	return retval;
}
