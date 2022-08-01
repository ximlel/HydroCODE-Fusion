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

#include "../include/file_io.h"


/*
 * To realize cross-platform programming.
 * MKDIR:  Create a subdirectory.
 * ACCESS: Determine access permissions for files or folders.
 */
#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#include <io.h>
#define MKDIR(a)  _mkdir((a))
#define ACCESS(a) _access((a),0)
#elif __linux__
#include <sys/stat.h>
#include <unistd.h>
#define MKDIR(a)  mkdir((a),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#define ACCESS(a) access((a),0)
#endif


/**
 * @brief      This function counts how many numbers are there in the initial data file. 
 * @param[in]  fp:  The pointer to the input file.
 * @param[in]  add: The address of the input file.
 * @return     The number of the numbers in the initial data file.
 *    @retval  -1: If the given number of column is not coincided with that in the data file.
 */
static int flu_var_count(FILE * fp, const char * add)
{
    int num = 0;  // Data number.
    int flg = 0; /* We read characters one by one from the data file.
		  * "flg" helps us to count.
		  * -# 1: when read a number-using character (0, 1, 2, ..., e, E, minus sign and dot).
		  * -# 0: when read a non-number-using character. 
		  */
    int ch;

    while((ch = getc(fp)) != EOF) // Count the data number.
	{
	    if (ch == 45 || ch == 46 || ch == 69 || ch == 101 || isdigit(ch))
		flg = 1;
	    else if (!isspace(ch))
		{
		    fprintf(stderr, "Input contains illegal character(ASCII=%d) in the file '%s'!\n", ch, add);
		    return 0;
		}
	    else if (flg) // Read in the space.
		{
		    num++;
		    flg = 0;
		}
	}
    
    rewind(fp);
    return num;
}

/** @brief This function produces folder path for data input or output.
 *  @param[in]  name:      Name of the test example.
 *  @param[out] add_mkdir: Folder path for data input or output.
 *  @param[in]  i_or_o:    Conversion parameters for data input/output.
 *                         - 0:             data output.
 *                         - else (e.g. 1): data input.
 */
void example_io(const char * name, char * add_mkdir, const int i_or_o)
{	
    if (i_or_o == 0)
	strcpy(add_mkdir, "../data_out/one-dim/");
    else
	strcpy(add_mkdir, "../data_in/one-dim/");
    strcat(add_mkdir, name);

    if (ACCESS(add_mkdir) == -1) // Folder does not exist.
	{
	    if (i_or_o == 0)
		{
		    if(MKDIR(add_mkdir)) // Creat a new folder.
			{
			    fprintf(stderr, "Output directory '%s' construction failed.\n", add_mkdir);
			    exit(1);
			}
		    else
			printf("Output directory '%s' constructed.\n", add_mkdir);				
		}
	    else
		{
		    fprintf(stderr, "Input directory is not exist!\n");
		    exit(1);
		}
	}
    strcat(add_mkdir, "/");
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
void _1D_initialize(const char * name, const char * add_in)
{
    char addRHO[FILENAME_MAX], addU[FILENAME_MAX], addP[FILENAME_MAX];
    // The address of the velocity/pressure/density file to read in.
    FILE * fp_U, * fp_P, * fp_rho; // The pointer to the above data files.
    int num_U = 0, num_P = 0, num_rho = 0;  // The number of the numbers in the above data files.
    char ch;
    int file_read_state;
    
    // Open the initial data files.
    strcpy(addRHO, add_in);
    strcat(addRHO, "RHO.txt");
    if((fp_rho = fopen(addRHO, "r")) == 0)
	{
	    strcpy(addRHO, add_in);
	    strcat(addRHO, "RHO.dat");
	}
    if((fp_rho = fopen(addRHO, "r")) == 0)
	{
	    printf("Cannot open initial data file RHO!\n");
	    exit(1);
	}
    num_rho = flu_var_count(fp_rho, addRHO);

    strcpy(addU, add_in);
    strcat(addU, "U.txt");
    if((fp_U = fopen(addU, "r")) == 0)
	{
	    strcpy(addU, add_in);
	    strcat(addU, "U.dat");
	}
    if((fp_U = fopen(addU, "r")) == 0)
	{
	    printf("Cannot open initial data file U!\n");
	    exit(1);
	}
    num_U = flu_var_count(fp_U, addU);
    if(num_U != num_rho)
	{
	    printf("U:Unequal! num_U=%d  num_rho=%d\n", num_rho, num_U);
	    exit(3);
	}

    strcpy(addP, add_in);
    strcat(addP, "P.txt");
    if((fp_P = fopen(addP, "r")) == 0)
	{
	    strcpy(addP, add_in);
	    strcat(addP, "P.dat");
	}
    if((fp_P = fopen(addP, "r")) == 0)
	{
	    printf("Cannot open initial data file P!\n");
	    exit(1);
	}
    num_P = flu_var_count(fp_P, addP);
    if(num_P != num_rho)
	{
	    printf("P:Unequal! num_rho=%d  num_P=%d\n", num_rho, num_P);
	    exit(3);
	}

    // Initializes the reading of data.
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

  printf("%s data initialized, m=%d\n", name, num_rho);
}


/**
 * @brief This function read the configuration data file, and
 *        store the configuration data in the array "config".
 * @details The parameters in the configuration data file are as follows.
 *          - config[0] is the constant of the perfect gas
 *          - config[1] is the length of the time step
 *          - config[2] is the spatial grid size
 *          - config[3] is the largest value can be seen as zero
 *          - config[4] is the maximal number of time steps
 *          - config[5] is the total time
 *          - config[6] is the CFL number
 * @param[out] config: Array of the configuration data.
 * @param[in]  name:   Name of the test example.
 * @param[in]  add_in: Adress of the initial data folder of the test example.
 */
void _1D_configurate(double * config, const char * name, const char * add_in)
{
    FILE * fp_data;
    char add[FILENAME_MAX];
    strcpy(add, add_in);
    strcat(add, "config.txt");
    
  // Open the configuration data file.
  if((fp_data = fopen(add, "r")) == 0)
  {
    printf("Cannot open configuration data file!\n");
    free(U0);
    free(P0);
    free(RHO0);
    exit(1);
  }

  int n_conf, state;

  // Read the configuration data file.
  if((n_conf = flu_var_count(fp_data, add)) != N_CONF)
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

  // Check the configuration data.
  if(config[0] < (1.0 + config[3]))
  {
    printf("The constant of the perfect gas(%lf) should be larger than 1.0.\n", config[0]);
    free(U0);
    free(P0);
    free(RHO0);
    exit(2);
  }
  if(config[6] > (1.0 - config[3]))
  {
    printf("The CFL number(%lf) should be smaller than 1.0.\n", config[6]);
    free(U0);
    free(P0);
    free(RHO0);
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
  printf("gamma      = %g\n", config[0]);
  printf("delta_t    = %g\n", config[1]);
  printf("delta_x    = %g\n", config[2]);
  printf("eps        = %g\n", config[3]);
  printf("time step  = %d\n", (int)config[4]);
  printf("total time = %g\n", config[5]);
  printf("CFL number = %g\n", config[6]);
}


/**
 * @brief This function write the solution into output files.
 * @note  It is quite simple so there will be no more comments.
 * @param[in] m: The number of spatial points in the output data.
 * @param[in] N: The number of time steps in the output data.
 * @param[in] RHO,U,P,Ene,X[]: Array of the density/velocity/pressure/energy/position data.
 * @param[in] cpu_time: Array of the CPU time recording.
 * @param[in] config:   Array of the configuration data.
 * @param[in] name:     Name of the test example.
 * @param[in] add_out:  Adress of the data output folder of the test example. 
 */
void _1D_file_write(const int m, const int N, 
                    double * RHO[], double * U[], double * P[], double * Ene[], double * X[], 
                    const double * cpu_time, const double * config, const char * name, const char * add_out)
{
  FILE * fp_write;
  char file_data[FILENAME_MAX] = "";
  char str_time[FILENAME_MAX];

  // Records the time when the program is running.
  struct tm * local_time;
  time_t t;
  t=time(NULL);
  local_time=localtime(&t);
  
//===================Write Output Data File=========================

  strcpy(file_data, add_out);
  strcat(file_data, "/RHO");
  //sprintf(str_time, "_%02d%02d%02d%02d%02d%02d", local_time->tm_year-100, local_time->tm_mon+1, local_time->tm_mday, local_time->tm_hour, local_time->tm_min, local_time->tm_sec);
  //strcat(file_data, str_time);
  strcat(file_data, ".dat");
  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  int j = 0, n = 0;
  for(n = 0; n <= N; ++n)
  {
    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.18lf\t", RHO[n][j]);
    fprintf(fp_write, "\n");
  }
  fclose(fp_write);


  strcpy(file_data, add_out);
  strcat(file_data, "/U");
  strcat(file_data, ".dat");
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

  
  strcpy(file_data, add_out);
  strcat(file_data, "/P");
  strcat(file_data, ".dat");
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

  
  strcpy(file_data, add_out);
  strcat(file_data, "/E");
  strcat(file_data, ".dat");
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

  
  strcpy(file_data, add_out);
  strcat(file_data, "/X");
  strcat(file_data, ".dat");
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

//======================Write Log File============================
  strcpy(file_data, add_out);
  strcat(file_data, "/log");
  strcat(file_data, ".dat");

  if((fp_write = fopen(file_data, "w")) == 0)
  {
    printf("Cannot open log output file!\n");
    exit(1);
  }

  fprintf(fp_write, "%s initialized with %d grids.\n\n", name, m);

  fprintf(fp_write, "Configurated:\n");
  fprintf(fp_write, "gamma = %g\n", config[0]);
  fprintf(fp_write, "tau   = %g\n", config[1]);
  fprintf(fp_write, "h     = %g\n", config[2]);
  fprintf(fp_write, "eps   = %g\n", config[3]);
  fprintf(fp_write, "N_t   = %d\n", (int)config[4]);
  fprintf(fp_write, "t_all = %g\n", config[5]);
  fprintf(fp_write, "CFL   = %g\n", config[6]);

  // fprintf(fp_write, "%d time steps computed.\n", N);
  /*
  double* sum = calloc(N + 1, sizeof(double));
  sum[0] = 0.0;
  fprintf(fp_write, "CPU time for each step:");
  for(n = 1; n <= N; ++n)
  {
    fprintf(fp_write, "%.18lf  ", cpu_time[n]);
    sum[n] = sum[n-1] + cpu_time[n];
  }
  fprintf(fp_write, "\nTotal CPU time at each step:");
  for(n = 1; n <= N; ++n)
    fprintf(fp_write, "%.18lf  ", sum[n]);
  free(sum);
  */
  fclose(fp_write);
}
