/**
 * @file  _1D_file_out.c
 * @brief This is a set of functions which control the readout of one-dimensional data.
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
	strcpy(file_data, add_out);					\
	strcat(file_data, "/");						\
	strcat(file_data, #v);						\
	strcat(file_data, ".dat");					\
	if((fp_write = fopen(file_data, "w")) == NULL)			\
	    {								\
		printf("Cannot open solution output file: %s!\n", #v);	\
		exit(1);						\
	    }								\
	for(k = 0; k < N; ++k)						\
	    {								\
		for(j = 0; j < m; ++j)					\
		    fprintf(fp_write, "%.10g\t", (v_print));		\
		fprintf(fp_write, "\n");				\
	    }								\
	fclose(fp_write);						\
    } while (0)

/**
 * @brief This function write the 1-D solution into output .dat files.
 * @note  It is quite simple so there will be no more comments.
 * @param[in] m: The number of spatial points in the output data.
 * @param[in] N: The number of time steps in the output data.
 * @param[in] CV:  Structure of grid variable data.
 * @param[in] X[]: Array of the coordinate data.
 * @param[in] cpu_time:  Array of the CPU time recording.
 * @param[in] name:      Name of the numerical results.
 * @param[in] time_plot: Array of the plotting time recording.
 */
void _1D_file_write(const int m, const int N, const struct cell_var_stru CV, 
                    double * X[], const double * cpu_time, const char * name, const double * time_plot)
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
    
    char file_data[FILENAME_MAX+40] = "";
    FILE * fp_write;

//===================Write Output Data File=========================

    int k, j;
    PRINT_NC(RHO, CV.RHO[k][j]);
    PRINT_NC(U, CV.U[k][j]);
    PRINT_NC(P, CV.P[k][j]);
    PRINT_NC(E, CV.E[k][j]);
    PRINT_NC(X, 0.5 * (X[k][j] + X[k][j+1]));

    strcpy(file_data, add_out);
    strcat(file_data, "/time_plot.dat");
    if((fp_write = fopen(file_data, "w")) == NULL)
	{
	    printf("Cannot open solution output file: time_plot!\n");
	    exit(1);
	}
    for(k = 0; k < N; ++k)
	fprintf(fp_write, "%.10g\n", time_plot[k]);
    fclose(fp_write);

//======================Write Log File============================
    config_write(add_out, cpu_time, name);
}
