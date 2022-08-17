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
 * @brief This function reads the initial data file to generate the initial data.
 * @param[in]  fp: The pointer to the input file.
 * @param[out]  U: The pointer to the data array.
 * @param[in] num: The number of the numbers in the input file. 
 * @return It returns 0 if successfully read the file,
 *         while returns the index of the wrong entry.
 */
int _1D_file_read(FILE * fp, double * U, const int num)
{
  int idx = 0, j = 0; // j is a frequently used index for spatial variables.
  char number[100]; // A string that stores a number.
  char ch;
  // int sign = 1;
  char *endptr;
    
  while((ch = getc(fp)) != EOF)
  {
    if(isspace(ch) && idx)
    {
      number[idx] = '\0';
      idx = 0;

      /* // Before format_string() and str2num() are deprecated.
      sign = format_string(number);
      if(!sign)
	return j+1;
      else if(j == num)
	return j;
      U[j] = sign * str2num(number);
      */
      errno = 0;
      U[j] = strtod(number,&endptr);
      if (errno == ERANGE || *endptr != '\0')
	  return j+1;
      else if(j == num)
	  return j;
      
      ++j;
    }
    else if((ch == 46) || (ch == 45) || (ch == 69) || (ch == 101) || isdigit(ch))
      number[idx++] = ch;
  }

  return 0;
}

/**
 * @brief This function read the configuration data file, and
 *        store the configuration data in the array "config".
 * @details The parameters in the configuration data file refer to 'doc/config.csv'.
 * @param[in]  name:   Name of the test example.
 * @param[in]  add_in: Adress of the initial data folder of the test example.
 */
static void _1D_configurate(const char * name)
{
    FILE * fp_data;
    char add[FILENAME_MAX];
    strcpy(add, add_in);
    strcat(add, "config.txt");
    
  // Open the configuration data file.
  if((fp_data = fopen(add, "r")) == NULL)
  {
    printf("Cannot open configuration data file!\n");
    exit(1);
  }

  // Read the configuration data file.
  if(config_read(fp_data) == 0)
      {
	  fclose(fp_data);
	  exit(2);
      }
  fclose(fp_data);

  printf("%s configurated:\n", name);
  // Check the configuration data.
  config_check();

    printf("delta_x    = %g\n", config[10]);
    printf("delta_t    = %g\n", config[16]);
    printf("bondary    = %d\n", (int)config[17]);
}


/** 
  * @brief      This function reads the initial data file of velocity/pressure/density.
  * @details    The function initialize the extern pointer RHO0/U0/P0 pointing to the
  *             position of a block of memory consisting (m+1) variables* of type double.
  *             The value of first of these variables is m.
  *             The following m variables are the initial value.
  * @param[in]  name: Name of the test example.
  * @param[in]  add_in:  Adress of the initial data folder of the test example.
  */
void _1D_initialize(const char * name)
{
  /* 
   * Read the configuration data.
   * The detail could be seen in the definition of array config
   * referring to file '_1D_f_io.c'.
   */
  _1D_configurate(name);

    char add_in[FILENAME_MAX];
    // Get the address of initial data files.
    example_io(name, add_in, 1);
  
    char addRHO[FILENAME_MAX], addU[FILENAME_MAX], addP[FILENAME_MAX];
    // The address of the velocity/pressure/density file to read in.
    FILE * fp_U, * fp_P, * fp_rho; // The pointer to the above data files.
    int num_U = 0, num_P = 0, num_rho = 0;  // The number of the numbers in the above data files.
    char ch;
    int file_read_state;
    
    // Open the initial data files.
    strcpy(addRHO, add_in);
    strcat(addRHO, "RHO.txt");
    if((fp_rho = fopen(addRHO, "r")) == NULL)
	{
	    strcpy(addRHO, add_in);
	    strcat(addRHO, "RHO.dat");
	}
    if((fp_rho = fopen(addRHO, "r")) == NULL)
	{
	    printf("Cannot open initial data file RHO!\n");
	    exit(1);
	}
    num_rho = flu_var_count(fp_rho, addRHO);

    strcpy(addU, add_in);
    strcat(addU, "U.txt");
    if((fp_U = fopen(addU, "r")) == NULL)
	{
	    strcpy(addU, add_in);
	    strcat(addU, "U.dat");
	}
    if((fp_U = fopen(addU, "r")) == NULL)
	{
	    printf("Cannot open initial data file U!\n");
	    exit(1);
	}
    num_U = flu_var_count(fp_U, addU);
    if(num_U != num_rho)
	{
	    printf("U:Unequal! num_U=%d  num_rho=%d\n", num_rho, num_U);
	    exit(2);
	}

    strcpy(addP, add_in);
    strcat(addP, "P.txt");
    if((fp_P = fopen(addP, "r")) == NULL)
	{
	    strcpy(addP, add_in);
	    strcat(addP, "P.dat");
	}
    if((fp_P = fopen(addP, "r")) == NULL)
	{
	    printf("Cannot open initial data file P!\n");
	    exit(1);
	}
    num_P = flu_var_count(fp_P, addP);
    if(num_P != num_rho)
	{
	    printf("P:Unequal! num_rho=%d  num_P=%d\n", num_rho, num_P);
	    exit(2);
	}

    // Initializes the reading of data.
  RHO0 = (double *)malloc((num_rho + 1) * sizeof(double));
  RHO0[0] = (double)num_rho;
  file_read_state = _1D_file_read(fp_rho, RHO0+1, num_rho);
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
  file_read_state = _1D_file_read(fp_U, U0+1, num_U);
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
  file_read_state = _1D_file_read(fp_P, P0+1, num_P);
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

  printf("%s data initialized, m=%d\n", name, num_rho);
}


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
