/**
 * @file  config_handle.c
 * @brief This is a set of functions which control the read-in of configuration data.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <errno.h>
#include <ctype.h>
#include <limits.h>

#include "../include/var_struc.h"

/*
 * To realize cross-platform programming.
 * ACCESS: Determine access permissions for files or folders.
 */
#ifdef _WIN32
#include <io.h>
/*
 * m=0: Test for existence.
 * m=2: Test for write permission.
 * m=4: Test for read permission.
 */
#define ACCESS(a,m) _access((a),(m))
#elif __linux__
#include <unistd.h>
#define ACCESS(a,m) access((a),(m))
#endif


/**
 * @brief This function check whether the configuration data is reasonable and set the default.
 *        All configuration data is in tha array 'config[]'.
 */
static void config_check(void)
{
    const int dim = (int)config[0];
    printf("  dimension\t= %d\n", dim);

    // Maximum number of time steps
    if(isfinite(config[1]) && config[1] >= 0.0)
	{
	    config[5] = isfinite(config[5]) ? config[5] : (double)INT_MAX;
	    printf("  total time\t= %g\n", config[1]);
	}
    else if(!isfinite(config[5]))
	{
	    fprintf(stderr, "The total time or the maximum number of time steps must be setted properly!\n");
	    exit(2);
	}
    else
	{
	    config[1] = INFINITY;
	    if(isfinite(config[16]))
		{
		    printf("  total time\t= %g * %d = %g\n", config[16], (int)config[5], config[16]*(int)config[5]);
		    printf("  delta_t\t= %g\n", config[16]);
		}
	}
    printf("  time step\t= %d\n", (int)config[5]);
	    
    if(isinf(config[4]))
	config[4] = EPS;
    double eps = config[4];
    if(eps < 0.0 || eps > 0.01)
	{
	    fprintf(stderr, "eps(%f) should in (0, 0.01)!\n", eps);
	    exit(2);
	}
    printf("  eps\t\t= %g\n", eps);

    if(isinf(config[6]))
	config[6] = 1.4;
    else if(config[6] < 1.0 + eps)
	{
	    fprintf(stderr, "The constant of the perfect gas(%f) should be larger than 1.0!\n", config[6]);
	    exit(2);
	}
    printf("  gamma\t\t= %g\n", config[6]);

    if (isinf(config[7]))
	{
	    switch(dim)
		{
		case 1:
		    config[7] = 0.9;  break;
		case 2:
		    config[7] = 0.45; break;
		}
	}
    else if(config[7] > 1.0 - eps)
	{
	    fprintf(stderr, "The CFL number(%f) should be smaller than 1.0.\n", config[7]);
	    exit(2);
	}
    printf("  CFL number\t= %g\n", config[7]);

    if(isinf(config[41]))
	config[41] = 1.9;
    else if(config[41] < -eps || config[41] > 2.0)
	{
	    fprintf(stderr, "The parameter in minmod limiter(%f) should in [0, 2)!\n", config[41]);
	    exit(2);
	}
  
    if(isinf(config[110]))
	config[110] = 0.72;	
    else if(config[110] < eps)
	{
	    fprintf(stderr, "The specific heat at constant volume(%f) should be larger than 0.0!\n", config[110]);
	    exit(2);
	}

    // Specie number
    config[2]   = isfinite(config[2])   ? config[2]   : (double)1;	
    // Coordinate framework (EUL/LAG/ALE)
    config[8]   = isfinite(config[8])   ? config[8]   : (double)0;
    // Reconstruction (prim_var/cons_var)
    config[31]  = isfinite(config[31])  ? config[31]  : (double)0;
    // Output initial data
    config[32]  = isfinite(config[32])  ? config[32]  : (double)true;
    // Dimensional splitting
    config[33]  = isfinite(config[33])  ? config[33]  : (double)false;
    // Slope limiter (Venkatakrishnan/Barth-Jespersen)
    config[40]  = isfinite(config[40])  ? config[40]  : (double)0;
    // Parameter α in minmod limiter
    config[41]  = isfinite(config[41])  ? config[41]  : 1.9;
    // Runge-Kutta time discretization
    config[53]  = isfinite(config[53])  ? config[53]  : (double)false;
    // Conservative variable (U_gamma) ργ
    config[60]  = isfinite(config[60])  ? config[60]  : (double)false;
    // v_fix: Shear velocity
    config[61]  = isfinite(config[61])  ? config[61]  : (double)0;
    // Offset of the upper and downside periodic boundary
    config[70]  = isfinite(config[70])  ? config[70]  : (double)0;
    // offset_x: Grid offset in x direction
    config[210] = isfinite(config[210]) ? config[210] : 0.0;
    // offset_y: Grid offset in y direction
    config[211] = isfinite(config[211]) ? config[211] : 0.0;
    // offset_z: Grid offset in z direction
    config[212] = isfinite(config[212]) ? config[212] : 0.0;
}

/**
 * @brief This function read the configuration data file, and
 *        store the configuration data in the array 'config[]'.
 * @param[in] fp: The pointer to the configuration data file.
 * @return    Configuration data file read status.
 *    @retval 1: Success to read in configuration data file. 
 *    @retval 0: Failure to read in configuration data file. 
 */
static int config_read(FILE * fp)
{	
	char one_line[200]; // String to store one line.
	char *endptr;
	double tmp;
	int i, line_num = 1; // Index of config[*], line number.

	while (fgets(one_line, sizeof(one_line), fp) != NULL)
		{
			// A line that doesn't begin with digits is a comment.
			i = strtoul(one_line, &endptr, 10);
			for ( ; isspace(*endptr); endptr++) ;

			// If the value of config[i] doesn't exit, it is 0 by default.
			if (0 < i && i < N_CONF)
				{
					errno = 0;
					tmp = strtod(endptr, NULL);
					if(errno == ERANGE)
					    {
						fprintf(stderr, "Value range error of %d-th configuration in line %d of configuration file!\n", i, line_num);
						return 1;
					    }
					else if(isinf(config[i]))
					    printf("%3d-th configuration: %g\n", i, config[i] = tmp);
					else if(fabs(config[i] - tmp) > EPS)
					    printf("%3d-th configuration is repeatedly assigned with %g and %g(abandon)!\n", i, config[i], tmp);
				}
			else if (i != 0 || (*endptr != '#' && *endptr != '\0'))
				fprintf(stderr, "Warning: unknown row occurrs in line %d of configuration file!\n", line_num);
			line_num++;
		}
	if (ferror(fp))
		{
			fprintf(stderr, "Read error occurrs in configuration file!\n");
			return 0;
		}
	return 1;
}


/**
 * @brief This function controls configuration data reading and validation.
 * @details The parameters in the configuration data file refer to 'doc/config.csv'.
 * @param[in] add_in: Adress of the initial data folder of the test example.
 */
void configurate(const char * add_in)
{
  FILE * fp_data;
  char add[FILENAME_MAX+40];
  strcpy(add, add_in);
  strcat(add, "config.txt");
  // Open the configuration data file.
  if((fp_data = fopen(add, "r")) == NULL)
      {
	  strcpy(add, add_in);
	  strcat(add, "config.dat");
      }
  if((fp_data = fopen(add, "r")) == NULL)
      {
	  printf("Cannot open configuration data file!\n");
	  perror(add_in);
	  exit(1);
      }

  // Read the configuration data file.
  if(config_read(fp_data) == 0)
      {
	  fclose(fp_data);
	  exit(2);
      }
  fclose(fp_data);

  printf("Configurated:\n");
  // Check the configuration data.
  config_check();
}


/**
 * @brief This function write configuration data and program record into the file 'log.dat'.
 * @details The parameters in the log file refer to 'doc/config.csv'.
 * @remark  The code of recording CPU time is commented out.
 * @param[in] add_out:  Address of the output data folder of the test example.
 * @param[in] cpu_time: Array of the CPU time recording.
 * @param[in] name:     Name of the test example.
 */
void config_write(const char * add_out, const double * cpu_time, const char * name)
{
    char file_data[FILENAME_MAX+40];
	const int dim = (int)config[0];
    FILE * fp_write;

//======================Write Log File============================
  strcpy(file_data, add_out);
  strcat(file_data, "/log");
  strcat(file_data, ".dat");
  if((fp_write = fopen(file_data, "w")) == NULL)
  {
    printf("Cannot open log output file!\n");
    exit(1);
  }

  fprintf(fp_write, "%s is initialized with %d grids.\n\n", name, (int)config[3]);
  fprintf(fp_write, "Configurated:\n");
  fprintf(fp_write, "dim\t\t= %d\n", dim);
  if(isfinite(config[1]))
      fprintf(fp_write, "t_all\t= %d\n", (int)config[1]);
  else if(isfinite(config[16]))
      fprintf(fp_write, "tau\t\t= %g\n", config[16]);
  fprintf(fp_write, "eps\t\t= %g\n", config[4]);
  fprintf(fp_write, "gamma\t= %g\n", config[6]);
  fprintf(fp_write, "CFL\t\t= %g\n", config[7]);
  fprintf(fp_write, "h\t\t= %g\n", config[10]);
  fprintf(fp_write, "bond\t= %d\n", (int)config[17]);
  if(dim == 2)
      {
	  fprintf(fp_write, "h_y\t\t= %g\n", config[11]);
	  fprintf(fp_write, "bond_y\t= %d\n", (int)config[18]);
      }
  fprintf(fp_write, "\nA total of %d time steps are computed.\n", (int)config[5]);
  /*
  double * sum = (double*)calloc(N, sizeof(double));
  sum[0] = 0.0;
  fprintf(fp_write, "CPU time for each step:");
  for(k = 1; k < N; ++k)
  {
    fprintf(fp_write, "%.18f  ", cpu_time[k]);
    sum[k] = sum[k-1] + cpu_time[k];
  }
  fprintf(fp_write, "\nTotal CPU time at each step:");
  for(k = 1; k < N; ++k)
    fprintf(fp_write, "%.18f  ", sum[k]);
  free(sum);
  sum = NULL;
  */
  fclose(fp_write);
}
