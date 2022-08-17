/**
 * @file  _1D_f_io.c
 * @brief This is a set of functions which control the reading and readout of
 *        one-dimensional data.
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
 * @param[in] RHO,U,P,Ene,X[]: Array of the density/velocity/pressure/energy/coordinate data.
 * @param[in] cpu_time: Array of the CPU time recording.
 * @param[in] name:     Name of the test example.
 * @param[in] add_out:  Adress of the data output folder of the test example. 
 */
void _1D_file_write(const int m, const int N, struct cell_var CV, double * X[], 
                    const double * cpu_time, const char * name)
{
  // Records the time when the program is running.
  /*
  struct tm * local_time;
  time_t t;
  t=time(NULL);
  local_time=localtime(&t);
  char str_time[FILENAME_MAX];
  sprintf(str_time, "_%02d%02d%02d%02d%02d%02d", local_time->tm_year-100, local_time->tm_mon+1, local_time->tm_mday, local_time->tm_hour, local_time->tm_min, local_time->tm_sec);
  */

  char file_data[FILENAME_MAX];
  FILE * fp_write;

//===================Write Output Data File=========================

  strcpy(file_data, add_out);
  strcat(file_data, "/RHO");
  strcat(file_data, ".dat");
  if((fp_write = fopen(file_data, "w")) == NULL)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  int j = 0, n = 0;
  for(n = 0; n < N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.10g\t", RHO[n][j]);
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);


  strcpy(file_data, add_out);
  strcat(file_data, "/U");
  strcat(file_data, ".dat");
  if((fp_write = fopen(file_data, "w")) == NULL)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(n = 0; n < N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.10g\t", U[n][j]);
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);

  
  strcpy(file_data, add_out);
  strcat(file_data, "/P");
  strcat(file_data, ".dat");
  if((fp_write = fopen(file_data, "w")) == NULL)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(n = 0; n < N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.10g\t", P[n][j]);
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);

  
  strcpy(file_data, add_out);
  strcat(file_data, "/E");
  strcat(file_data, ".dat");
  if((fp_write = fopen(file_data, "w")) == NULL)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(n = 0; n < N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.10g\t", Ene[n][j]);
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);

  
  strcpy(file_data, add_out);
  strcat(file_data, "/X");
  strcat(file_data, ".dat");
  if((fp_write = fopen(file_data, "w")) == NULL)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(n = 0; n < N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.10g\t", 0.5 * (X[n][j] + X[n][j+1]));
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);

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
      fprintf(fp_write, "t_all = %g\n", config[1]);
  else if (!isinf(config[5])
      {
	  fprintf(fp_write, "N_t   = %d\n", (int)config[5]);
	  fprintf(fp_write, "tau   = %g\n", config[16]);
      }
  fprintf(fp_write, "eps   = %g\n", config[4]);
  fprintf(fp_write, "gamma = %g\n", config[6]);
  fprintf(fp_write, "CFL   = %g\n", config[7]);
  fprintf(fp_write, "h     = %g\n", config[10]);
  fprintf(fp_write, "bond  = %d\n", (int)config[17]);
  fprintf(fp_write, "%d time steps computed.\n", (int)config[5]);
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
