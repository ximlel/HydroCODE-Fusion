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


void file_write(int n_x, int n_y, int N, double *rho[][n_x], double *u[][n_x], double *v[][n_x], double *p[][n_x], double * cpu_time, double * config, double tau, char * scheme)
{
  FILE * fp_write;
  char file_data[100] = "";
  char str_time[100] = "";

  struct tm * local_time;
  time_t t;
  t=time(NULL);
  local_time=localtime(&t);

//===================Write solution File=========================

  strcat(file_data, "./data_out/rho_");
  strcat(file_data, scheme);
  sprintf(str_time, "_%02d%02d%02d%02d%02d%02d", local_time->tm_year-100, local_time->tm_mon+1, local_time->tm_mday, local_time->tm_hour, local_time->tm_min, local_time->tm_sec);
  //strcat(file_data, str_time);
  strcat(file_data, ".txt");

  printf("%s\n", str_time);

  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  int j = 0, i = 0, k = 0;
  for(k = 0; k <= N; ++k)
  {
    for(i = 0; i < n_y; ++i)
    {
      //if(i % 2)
	//continue;
      //if(i)
	//continue;
      for(j = 0; j < n_x; ++j)
      {
	//if(j % 2)
	  //continue;
	fprintf(fp_write, "%.18lf\t", rho[k][j][i]);
      }
      fprintf(fp_write, "\n");
    }
  }
  fclose(fp_write);


  file_data[11] = 'u';
  file_data[12] = '_';
  file_data[13] = '_';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(k = 0; k <= N; ++k)
  {
    for(i = 0; i < n_y; ++i)
    {
      //if(i % 2)
	//continue;
      //if(i)
	//continue;
      for(j = 0; j < n_x; ++j)
      {
	//if(j % 2)
	  //continue;
	fprintf(fp_write, "%.18lf\t", u[k][j][i]);
      }
      fprintf(fp_write, "\n");
    }
  }
  fclose(fp_write);


  file_data[11] = 'v';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(k = 0; k <= N; ++k)
  {
    for(i = 0; i < n_y; ++i)
    {
      //if(i % 2)
	//continue;
      //if(i)
	//continue;
      for(j = 0; j < n_x; ++j)
      {
	//if(j % 2)
	  //continue;
	fprintf(fp_write, "%.18lf\t", v[k][j][i]);
      }
      fprintf(fp_write, "\n");
    }
  }
  fclose(fp_write);


  file_data[11] = 'p';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(k = 0; k <= N; ++k)
  {
    for(i = 0; i < n_y; ++i)
    {
      //if(i % 2)
	//continue;
      //if(i)
      //continue;
      for(j = 0; j < n_x; ++j)
      {
	//if(j % 2)
	//continue;
	fprintf(fp_write, "%.18lf\t", p[k][j][i]);
      }
      fprintf(fp_write, "\n");
    }
  }
  fclose(fp_write);
//===================Write LOG File=========================
  file_data[11] = 'l';
  file_data[12] = 'o';
  file_data[13] = 'g';


  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open log output file!\n");
    exit(1);
  }

  double sum[N+1];
  sum[0] = 0.0;

  fprintf(fp_write, "%s initialized with %d x-grids, %d y-grids.\n\n", scheme, n_x, n_y);

  fprintf(fp_write, "Configurated:\n");
  fprintf(fp_write, "gamma = %g\n", config[0]);
  fprintf(fp_write, "CFL   = %g\n", config[1]);
  fprintf(fp_write, "h_x   = %g\n", config[2]);
  fprintf(fp_write, "h_y   = %g\n", config[3]);
  fprintf(fp_write, "tau   = %g\n", tau);
  fprintf(fp_write, "eps   = %g\n", config[4]);
  fprintf(fp_write, "tim   = %d\n\n", (int)config[5]);
  fprintf(fp_write, "%d time steps computed.\n", N);
  fprintf(fp_write, "CPU time for each step:");
  for(k = 1; k <= N; ++k)
  {
    fprintf(fp_write, "%.18lf  ", cpu_time[k]);
    sum[k] = sum[k-1] + cpu_time[k];
  }
  fprintf(fp_write, "\nTotal CPU time at each step:");
  for(k = 1; k <= N; ++k)
    fprintf(fp_write, "%.18lf  ", sum[k]);

  fclose(fp_write);
}





void file_write_noneval(int m, int n, int N, double *rho[m], double cpu_time, double * config, double tau, char * scheme)
{
  FILE * fp_write;
  char file_data[100] = "";
  char str_time[100] = "";

  struct tm * local_time;
  time_t t;
  t=time(NULL);
  local_time=localtime(&t);

//===================Write solution File=========================

  strcat(file_data, "./data_out/rho_");
  strcat(file_data, scheme);
  //sprintf(str_time, "_%02d%02d%02d%02d%02d%02d", local_time->tm_year-100, local_time->tm_mon+1, local_time->tm_mday, local_time->tm_hour, local_time->tm_min, local_time->tm_sec);
  //strcat(file_data, str_time);
  strcat(file_data, ".txt");

  //printf("%s\n", str_time);

  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  int j = 0, i = 0;

  for(i = 0; i < n; ++i)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.18lf  ", rho[j][i]);
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);
  printf("ha\n");
//===================Write LOG File=========================
  file_data[11] = 'l';
  file_data[12] = 'o';
  file_data[13] = 'g';


  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open log output file!\n");
    exit(1);
  }

  fprintf(fp_write, "%s initialized with %d x-grids, %d y-grids.\n\n", scheme, m, n);

  fprintf(fp_write, "Configurated:\n");
  fprintf(fp_write, "gamma = %g\n", config[0]);
  fprintf(fp_write, "CFL   = %g\n", config[1]);
  fprintf(fp_write, "h_x   = %g\n", config[2]);
  fprintf(fp_write, "h_y   = %g\n", config[3]);
  fprintf(fp_write, "tau   = %g\n", tau);
  fprintf(fp_write, "eps   = %g\n", config[4]);
  fprintf(fp_write, "tim   = %d\n\n", (int)config[5]);
  fprintf(fp_write, "%d time steps computed.\n", N);
  fprintf(fp_write, "CPU time: %g", cpu_time);

  fclose(fp_write);
}




void file_write_y(int m, int n, int N, double *rho[][m], double *u[][m], double *v[][m], double *p[][m], double * cpu_time, double * config, double tau, char * scheme)
{
  FILE * fp_write;
  char file_data[100] = "";
  char str_time[100] = "";

  struct tm * local_time;
  time_t t;
  t=time(NULL);
  local_time=localtime(&t);

//===================Write solution File=========================

  strcat(file_data, "./data_out/rho_");
  strcat(file_data, scheme);
  sprintf(str_time, "_%02d%02d%02d%02d%02d%02d", local_time->tm_year-100, local_time->tm_mon+1, local_time->tm_mday, local_time->tm_hour, local_time->tm_min, local_time->tm_sec);
  //strcat(file_data, str_time);
  strcat(file_data, ".txt");

  printf("%s\n", str_time);

  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  int j = 0, i = 0, k = 0;
  for(k = 0; k <= N; ++k)
  {
    for(i = 0; i < n; ++i)
    {
      fprintf(fp_write, "%.18lf\t", rho[k][0][i]);
    }
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);


  file_data[11] = 'u';
  file_data[12] = '_';
  file_data[13] = '_';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(k = 0; k <= N; ++k)
  {
    for(i = 0; i < n; ++i)
    {
      fprintf(fp_write, "%.18lf\t", u[k][0][i]);
    }
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);


  file_data[11] = 'v';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(k = 0; k <= N; ++k)
  {
    for(i = 0; i < n; ++i)
    {
      fprintf(fp_write, "%.18lf\t", v[k][0][i]);
    }
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);


  file_data[11] = 'p';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(k = 0; k <= N; ++k)
  {
    for(i = 0; i < n; ++i)
    {
      fprintf(fp_write, "%.18lf\t", p[k][0][i]);
    }
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);
//===================Write LOG File=========================
  file_data[11] = 'l';
  file_data[12] = 'o';
  file_data[13] = 'g';


  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open log output file!\n");
    exit(1);
  }

  double sum[N+1];
  sum[0] = 0.0;

  fprintf(fp_write, "%s initialized with %d x-grids, %d y-grids.\n\n", scheme, m, n);

  fprintf(fp_write, "Configurated:\n");
  fprintf(fp_write, "gamma = %g\n", config[0]);
  fprintf(fp_write, "CFL   = %g\n", config[1]);
  fprintf(fp_write, "h_x   = %g\n", config[2]);
  fprintf(fp_write, "h_y   = %g\n", config[3]);
  fprintf(fp_write, "tau   = %g\n", tau);
  fprintf(fp_write, "eps   = %g\n", config[4]);
  fprintf(fp_write, "tim   = %d\n\n", (int)config[5]);
  fprintf(fp_write, "%d time steps computed.\n", N);
  fprintf(fp_write, "CPU time for each step:");
  for(k = 1; k <= N; ++k)
  {
    fprintf(fp_write, "%.18lf  ", cpu_time[k]);
    sum[k] = sum[k-1] + cpu_time[k];
  }
  fprintf(fp_write, "\nTotal CPU time at each step:");
  for(k = 1; k <= N; ++k)
    fprintf(fp_write, "%.18lf  ", sum[k]);

  fclose(fp_write);
}
