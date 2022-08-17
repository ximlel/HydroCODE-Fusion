/**
 * @file  _1D_file_out.c
 * @brief This is a set of functions which control the readout of one-dimensional data.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"

/**
 * @brief This function write the solution into output files.
 * @note  It is quite simple so there will be no more comments.
 * @param[in] m: The number of spatial points in the output data.
 * @param[in] N: The number of time steps in the output data.
 * @param[in,out] CV:   structural body of grid variable data.
 * @param[in,out] X[]:  Array of the coordinate data.
 * @param[in] cpu_time: Array of the CPU time recording.
 * @param[in] name:     Name of the test example.
 */

#define PRINT_NP(v, v_print)						\
    do {								\
	strcpy(file_data, add_out);					\
	strcat(file_data, "/");						\
	strcat(file_data, #v);						\
	strcat(file_data, ".dat");					\
	if((fp_write = fopen(file_data, "w")) == NULL)			\
	    {								\
		printf("Cannot open solution output file: %s!\n", #v);	\
		exit(1);						\
	    }								\
	for(n = 0; n < N; ++n)						\
	    {								\
		for(j = 0; j < m; ++j)					\
		    fprintf(fp_write, "%.10g\t", (v_print));		\
		fprintf(fp_write, "\n");				\
	    }								\
	fclose(fp_write);						\
    } while (0)

void _1D_file_write(const int m, const int N, struct cell_var CV, double * X[], 
                    const double * cpu_time, const char * name)
{
  // Records the time when the program is running.
  /*
  struct tm * local_time;
  time_t t;
  t=time(NULL);
  local_time=localtime(&t);
  char str_time[100];
  sprintf(str_time, "_%02d%02d%02d%02d%02d%02d", local_time->tm_year-100, local_time->tm_mon+1, local_time->tm_mday, local_time->tm_hour, local_time->tm_min, local_time->tm_sec);
  */

    char add_out[FILENAME_MAX+40];
    // Get the address of the output data folder of the test example.
    example_io(name, add_out, 0);
    
    char file_data[FILENAME_MAX+40];
    FILE * fp_write;

//===================Write Output Data File=========================

    int j, n;
    PRINT_NP(RHO, CV.RHO[n][j]);
    PRINT_NP(U, CV.U[n][j]);
    PRINT_NP(P, CV.P[n][j]);
    PRINT_NP(E, CV.E[n][j]);
    PRINT_NP(X, 0.5 * (X[n][j] + X[n][j+1]));

//======================Write Log File============================
  strcpy(file_data, add_out);
  strcat(file_data, "/log");
  strcat(file_data, ".dat");
  if((fp_write = fopen(file_data, "w")) == NULL)
  {
    printf("Cannot open log output file!\n");
    exit(1);
  }

  fprintf(fp_write, "%s is initialized with %d grids.\n\n", name, m);
  fprintf(fp_write, "Configurated:\n");
  if(!isinf(config[1]))
      fprintf(fp_write, "t_all = %g\n", (int)config[1]);
  else if(!isinf(config[5]))
      {
	  fprintf(fp_write, "N_t   = %d\n", (int)config[5]);
	  fprintf(fp_write, "tau   = %g\n", config[16]);
      }
  fprintf(fp_write, "dim   = %g\n", config[3]);
  fprintf(fp_write, "eps   = %g\n", config[4]);
  fprintf(fp_write, "gamma = %g\n", config[6]);
  fprintf(fp_write, "CFL   = %g\n", config[7]);
  fprintf(fp_write, "h     = %g\n", config[10]);
  fprintf(fp_write, "bond  = %d\n", (int)config[17]);
  fprintf(fp_write, "\n%d time steps computed.\n", (int)config[5]);
  /*
  double* sum = calloc(N, sizeof(double));
  sum[0] = 0.0;
  fprintf(fp_write, "CPU time for each step:");
  for(n = 1; n < N; ++n)
  {
    fprintf(fp_write, "%.18f  ", cpu_time[n]);
    sum[n] = sum[n-1] + cpu_time[n];
  }
  fprintf(fp_write, "\nTotal CPU time at each step:");
  for(n = 1; n < N; ++n)
    fprintf(fp_write, "%.18f  ", sum[n]);
  free(sum);
  */
  fclose(fp_write);
}
