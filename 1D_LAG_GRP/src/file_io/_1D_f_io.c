#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>


#include "../include/file_io.h"


#ifndef N_CONF
#define N_CONF 5
#endif /* N_CONF */


extern double * U0;
extern double * P0;
extern double * RHO0;



/* this function counts how many numbers are there
 * the initial data file.
 */
int _1D_file_pre_read(FILE * fp)
{
  int num = 0;
  int flag = 0;  /* We need to know how many numbers are there in
		  * the initial data. "flag" helps us to count.
		  * We read characters one by one from the data
		  * file. The value of "flag" is 1 when read a
		  * number-using character (i.e. 0, 1, 2, and so
		  * on and the dot), while is 0 when read a 
		  * non-number-using character. 
		  */
  char ch;

  while((ch = getc(fp)) != EOF)
  {
    //if(((ch < 45) || ((ch > 57) && (ch != 69) && (ch != 101)) || (ch == 47)) && (flag))
    if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (flag))
    {
      ++num;
      flag = 0;
    }
    else if( ((ch == 46)||(ch == 45)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (!flag) )
      flag = 1;
    else if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (!flag))
      continue;
    else if( ((ch == 46)||(ch == 45)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (flag) )
      continue;
    else
    {
      printf("Input contains illigal character(ASCII=%d, flag=%d)!\n", (int)ch, flag);
      if(!U0)
	free(U0);
      exit(1);
    }
  }

  return num;
}


/* This function reads the initial data file. The function 
 * initialize return a pointer pointing to the position of
 * a block of memory consisting (m+1) variables* of type
 * double. The value of first of these variables is m.
 * The following m variables are the initial value.
 */
void _1D_initialize(char * name, char * addRHO, char * addU, char * addP)
{
  FILE * fp_U, * fp_P, * fp_rho;
  int num_U = 0, num_P = 0, num_rho = 0;  
  char  ch;
  int file_read_state;
  //double * U0;


  //open the initial data file
  if((fp_rho = fopen(addRHO, "r")) == 0)
  {
    printf("Cannot open initial data file rho!\n");
    exit(1);
  }
  num_rho = _1D_file_pre_read(fp_rho);
  fseek(fp_rho, 0L, SEEK_SET);

  if((fp_U = fopen(addU, "r")) == 0)
  {
    printf("Cannot open initial data file U!\n");
    exit(1);
  }
  num_U = _1D_file_pre_read(fp_U);
  if(num_U != num_rho)
  {
    printf("U:Unequal! num_U=%d  num_rho=%d\n", num_rho, num_U);
    printf("Unequal!\n");
    exit(3);
  }
  fseek(fp_U, 0L, SEEK_SET);

  if((fp_P = fopen(addP, "r")) == 0)
  {
    printf("Cannot open initial data file P!\n");
    exit(1);
  }
  num_P = _1D_file_pre_read(fp_P);
  if(num_P != num_rho)
  {
    printf("P:Unequal! num_rho=%d  num_P=%d\n", num_rho, num_P);
    exit(3);
  }
  fseek(fp_P, 0L, SEEK_SET);


  RHO0 = (double *)malloc((num_rho + 1) * sizeof(double));
  RHO0[0] = (double)num_rho;
  file_read_state = file_read(fp_rho, RHO0+1, num_rho);
  fclose(fp_rho);
  if(file_read_state)
  {
    free(RHO0);
    if(file_read_state == num_rho)
      printf("Error on file reading!\n");
    else
      printf("\nThe %dth entry in the file RHO0 is not a number.\n", file_read_state);
    exit(2);
  }

  U0 = (double *)malloc((num_U + 1) * sizeof(double));
  U0[0] = (double)num_U;
  file_read_state = file_read(fp_U, U0+1, num_U);
  fclose(fp_U);
  if(file_read_state)
  {
    free(RHO0);
    free(U0);
    if(file_read_state == num_U)
      printf("Error on file reading!\n");
    else
      printf("\nThe %dth entry in the file U0 is not a number.\n", file_read_state);
    exit(2);
  }

  P0 = (double *)malloc((num_P + 1) * sizeof(double));
  P0[0] = (double)num_P;
  file_read_state = file_read(fp_P, P0+1, num_P);
  fclose(fp_P);
  if(file_read_state)
  {
    free(RHO0);
    free(U0);
    free(P0);
    if(file_read_state == num_P)
      printf("Error on file reading!\n");
    else
      printf("\nThe %dth entry in the file P0 is not a number.\n", file_read_state);
    exit(2);
  }


  printf("%s initialized, m=%d\n", name, num_rho);
}


/* This function read the configuration data file,
 * and store the configuration data in the array
 * "config".
 * config[0] is the constant of the perfect gas
 * config[1] is the length of the time step
 * config[2] is the spatial grid size
 * config[3] is the largest value can be seen as zero
 * config[4] is the number of time steps
 */
void _1D_configurate(double * config, char * name, char * add)
{
  FILE * fp_data;

  //open the configuration data file
  //printf("%s will open the cinfiguration data file: ", name);
  //printf("%s\n\n", add);
  if((fp_data = fopen(add, "r")) == 0)
  {
    printf("Cannot open configuration data file!\n");
    free(U0);
    free(P0);
    free(RHO0);
    exit(1);
  }

  int n_conf, state;

  //read the configuration data file
  if((n_conf = _1D_file_pre_read(fp_data)) != N_CONF)
  {
    printf("Configuration data file error, n_config=%d.\n", n_conf);
    free(U0);
    free(P0);
    free(RHO0);
    exit(2);
  }

  fseek(fp_data, 0L, SEEK_SET);
  state = file_read(fp_data, config, n_conf);
  fclose(fp_data);

  if(state)
  {
    printf("Configuration data file error, at %d.\n", state);
    free(U0);
    free(P0);
    free(RHO0);
    exit(2);
  }
  if(config[0] < (1.0 + config[3]))
  {
    printf("The constant of the perfect gas(%lf) should be larger than 1.0.\n", config[0]);
    free(U0);
    exit(2);
  }
  if(config[3] < 0)
  {
    printf("eps(%lf) should be positive.\n", config[3]);
    free(U0);
    free(P0);
    free(RHO0);
    exit(2);
  }

  printf("%s configurated:\n", name);
  printf("gamma = %g\n", config[0]);
  printf("tau   = %g\n", config[1]);
  printf("h     = %g\n", config[2]);
  printf("eps   = %g\n", config[3]);
  printf("time  = %d\n", (int)config[4]);
}


/* This function write the solution into an output file.
 * It is quite simple so there will be no more comments.
 */
void _1D_file_write(int m, int N, double * RHO[], double * U[], double * P[], double * Ene[], double * X[], double * cpu_time, double * config, char * scheme)
{
  FILE * fp_write;
  char file_data[100] = "";
  char str_time[100];

  struct tm * local_time;
  time_t t;
  t=time(NULL);
  local_time=localtime(&t);

  strcat(file_data, "./data_out/RHO_");
  strcat(file_data, scheme);
  sprintf(str_time, "_%02d%02d%02d%02d%02d%02d", local_time->tm_year-100, local_time->tm_mon+1, local_time->tm_mday, local_time->tm_hour, local_time->tm_min, local_time->tm_sec);
  //strcat(file_data, str_time);
  strcat(file_data, ".txt");

  //printf("%s\n", str_time);

  int j = 0, n = 0;


  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(n = 0; n <= N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.18lf\t", RHO[n][j]);
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);



  file_data[11] = 'U';
  file_data[12] = '_';
  file_data[13] = '_';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(n = 0; n <= N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.18lf\t", U[n][j]);
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);

  file_data[11] = 'P';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(n = 0; n <= N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.18lf\t", P[n][j]);
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);


  file_data[11] = 'E';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(n = 0; n <= N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.18lf\t", Ene[n][j]);
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);

  file_data[11] = 'X';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(n = 0; n <= N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.18lf\t", 0.5 * (X[n][j] + X[n][j+1]));
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);

  file_data[11] = 'A';
  file_data[12] = 'L';
  file_data[13] = 'L';
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  for(j = 0; j < m; ++j)
  {
	for(n = 0; n <= N; ++n)
      fprintf(fp_write, "%.18lf\t", RHO[n][j]);    
	for(n = 0; n <= N; ++n)
      fprintf(fp_write, "%.18lf\t", U[n][j]);
  	for(n = 0; n <= N; ++n)
      fprintf(fp_write, "%.18lf\t", P[n][j]);
 	for(n = 0; n <= N; ++n)
      fprintf(fp_write, "%.18lf\t", Ene[n][j]);
        for(n = 0; n <= N; ++n)
      fprintf(fp_write, "%.18lf\t", 0.5 * (X[n][j] + X[n][j+1]));
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

  fprintf(fp_write, "%s initialized with %d grids.\n\n", scheme, m);

  fprintf(fp_write, "Configurated:\n");
  fprintf(fp_write, "q   = %g\n", config[0]);
  fprintf(fp_write, "tau = %g\n", config[1]);
  fprintf(fp_write, "h   = %g\n", config[2]);
  fprintf(fp_write, "eps = %g\n", config[3]);
  fprintf(fp_write, "tim = %d\n\n", (int)config[4]);

  fprintf(fp_write, "%d time steps computed.\n", N);
  fprintf(fp_write, "CPU time for each step:");
  for(n = 1; n <= N; ++n)
  {
    fprintf(fp_write, "%.18lf  ", cpu_time[n]);
    sum[n] = sum[n-1] + cpu_time[n];
  }
  fprintf(fp_write, "\nTotal CPU time at each step:");
  for(n = 1; n <= N; ++n)
    fprintf(fp_write, "%.18lf  ", sum[n]);

  fclose(fp_write);
}
