/**
 * @file  file_1D_in.c
 * @brief This is a set of functions which control the read-in of one-dimensional data.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"


/**
 * @brief Count out and read in 1-D data of the initial fluid variable 'sfv'.
 */
#define STR_FLU_INI(sfv)						\
    do {								\
	strcpy(add, add_in);						\
	strcat(add, #sfv ".txt");					\
	if((fp = fopen(add, "r")) == NULL)				\
	    {								\
		strcpy(add, add_in);					\
		strcat(add, #sfv ".dat");				\
	    }								\
	if((fp = fopen(add, "r")) == NULL)				\
	    {								\
		printf("Cannot open initial data file: %s!\n", #sfv);	\
		exit(1);						\
	    }								\
	num_cell = flu_var_count(fp, add);				\
	if (num_cell < 1)						\
	    {								\
		printf("Error in counting fluid variables in initial data file: %s!\n", #sfv); \
		fclose(fp);						\
		exit(2);						\
	    }								\
	if(isinf(config[3]))						\
	    config[3] = (double)num_cell;				\
	else if(num_cell != (int)config[3])				\
	    {								\
		printf("Input unequal! num_%s=%d, num_cell=%d.\n", #sfv, num_cell, (int)config[3]); \
		exit(2);						\
	    }								\
	FV0.sfv = (double*)malloc(num_cell * sizeof(double));			\
	if(FV0.sfv == NULL)						\
	    {								\
		printf("NOT enough memory! %s\n", #sfv);		\
		exit(5);						\
	    }								\
	if(flu_var_read(fp, FV0.sfv, num_cell))				\
	    {								\
		fclose(fp);						\
		exit(2);						\
	    }								\
	fclose(fp);							\
    } while(0)

/** 
  * @brief      This function reads the 1-D initial data file of density/velocity/pressure
  *             and performs some other initialization function procedures.
  * @details    The function initialize the extern pointer FV0.RHO/U/P pointing to the
  *             position of a block of memory consisting m variables* of type double.
  *             These m variables are the initial value and the value of m is stored in config[3].
  * @param[in]  name:      Name of the test example.
  * @param[out] N_plot:    Pointer to the number of time steps for plotting.
  * @param[out] time_plot: Pointer to the array of the plotting time recording.
  * @return  \b FV0:  Structure of initial fluid variable data array pointer.
  * @note This function contains the function procedures 'time_plot_read()' and 'configurate()'.
  */
struct flu_var initialize_1D(const char * name, int * N_plot, double * time_plot[])
{
    struct flu_var FV0; // Structure of initial data array pointer.

    char add_in[FILENAME_MAX+40]; 
    // Get the address of the initial data folder of the test example.
    example_io(name, add_in, 1);
    
    /* 
     * Read the configuration data.
     * The detail could be seen in the definition of array config
     * referring to file 'doc/config.csv'.
     */
    configurate(add_in);
    printf("  delta_x\t= %g\n", config[10]);
    printf("  bondary\t= %d\n", (int)config[17]);

    time_plot_read(add_in, N_plot, time_plot);

    char add[FILENAME_MAX+40]; // The address of the velocity/pressure/density file to read in.
    FILE * fp;    // The pointer to the above data files.
    int num_cell; // The number of the numbers in the above data files.

    // Open the initial data files and initializes the reading of data.
    STR_FLU_INI(RHO);
    STR_FLU_INI(U);
    STR_FLU_INI(P);

    printf("'%s' data initialized, grid cell number = %d.\n", add_in, num_cell);
    return FV0;
}
