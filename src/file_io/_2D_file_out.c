/**
 * @file  _2D_file_out.c
 * @brief This is a set of functions which control the readout of two-dimensional data.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"


/**
 * @brief Print out fluid variable 'v' with array data element 'v_print'.
 */
#define PRINT_NC(v, v_print)						\
    do {								\
    strcpy(file_data, add_out);						\
    strcat(file_data, "/");						\
    strcat(file_data, #v);						\
    strcat(file_data, ".dat");						\
    if((fp_write = fopen(file_data, "w")) == NULL)			\
	{								\
	    printf("Cannot open solution output file: %s!\n", #v);	\
	    exit(1);							\
	}								\
    for(k = 0; k <= N; ++k)						\
	{								\
	    for(i = 0; i < n_y; ++i)					\
		{							\
		    for(j = 0; j < n_x; ++j)				\
			fprintf(fp_write, "%.10g\t", (v_print));	\
		    fprintf(fp_write, "\n");				\
		}							\
	    fprintf(fp_write, "\n\n");					\
	}								\
    fclose(fp_write);							\
    } while (0)

void _2D_file_write(const int n_x, const int n_y, const int N, struct cell_var_stru * CV,
		    double * X[], double * Y[], double * cpu_time, const char * name)
{
    char add_out[FILENAME_MAX+40];
    // Get the address of the output data folder of the test example.
    example_io(name, add_out, 0);
    
    char file_data[FILENAME_MAX+40] = "";
    FILE * fp_write;

//===================Write Solution File=========================

    int k, i, j;
    PRINT_NC(RHO, (CV+k)->RHO[j][i]);
    PRINT_NC(U, (CV+k)->U[j][i]);
    PRINT_NC(V, (CV+k)->V[j][i]);
    PRINT_NC(P, (CV+k)->P[j][i]);
    PRINT_NC(E, (CV+k)->E[j][i]);
    PRINT_NC(X, 0.25 * (X[j][i] + X[j][i+1] + X[j+1][i] + X[j+1][i+1]));
    PRINT_NC(Y, 0.25 * (Y[j][i] + Y[j][i+1] + Y[j+1][i] + Y[j+1][i+1]));

//===================Write LOG File=========================
  strcpy(file_data, add_out);
  strcat(file_data, "/log");
  strcat(file_data, ".dat");
  if((fp_write = fopen(file_data, "w")) == NULL)
  {
    printf("Cannot open log output file!\n");
    exit(1);
  }

  fprintf(fp_write, "%s is initialized with %d grids.\n\n", name, n_x * n_y);
  fprintf(fp_write, "Configurated:\n");
  fprintf(fp_write, "dim\t\t= %d\n", (int)config[0]);
  if(isfinite(config[1]))
      fprintf(fp_write, "t_all\t= %d\n", (int)config[1]);
  else if(isfinite(config[16]))
      fprintf(fp_write, "tau\t\t= %g\n", config[16]);
  fprintf(fp_write, "eps\t\t= %g\n", config[4]);
  fprintf(fp_write, "gamma\t= %g\n", config[6]);
  fprintf(fp_write, "CFL\t\t= %g\n", config[7]);
  fprintf(fp_write, "n_x\t= %g\n", config[10]);
  fprintf(fp_write, "n_y\t= %g\n", config[11]);
  fprintf(fp_write, "bond\t= %d\n", (int)config[17]);
  fprintf(fp_write, "bond_y\t= %d\n", (int)config[18]);
  fprintf(fp_write, "\nA total of %d time steps are computed.\n", (int)config[5]);

  fclose(fp_write);
}
