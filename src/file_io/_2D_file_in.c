/**
 * @file  _2D_file_in.c
 * @brief This is a set of functions which control the read-in of two-dimensional data.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"


/**
 * @brief Count out and read in 2-D data of the initial fluid variable 'sfv'.
 */
#define STR_FLU_INI(sfv)						\
    do {								\
    strcpy(add, add_in);						\
    strcat(add, #sfv ".txt");						\
    if((fp = fopen(add, "r")) == NULL)					\
	{								\
	    strcpy(add, add_in);					\
	    strcat(add, #sfv ".dat");					\
	}								\
    if((fp = fopen(add, "r")) == NULL)					\
	{								\
	    printf("Cannot open initial data file: %s!\n", #sfv);	\
	    exit(1);							\
	}								\
    line = flu_var_count_line(fp, add, &column);			\
    num_cell = line * column;						\
    if (num_cell < 1)							\
	{								\
	    printf("Error in counting fluid variables in initial data file: %s!\n", #sfv); \
	    fclose(fp);							\
	    exit(2);							\
	}								\
    if(isinf(config[3]))						\
	config[3] = (double)num_cell;					\
    if(isinf(config[13]))						\
	config[13] = (double)column;					\
    if(isinf(config[14]))						\
	config[14] = (double)line;					\
    else if(num_cell != (int)config[3] || column != (int)config[13] || line != (int)config[14])	\
	{								\
	    printf("Input unequal! num_%s=%d, num_cell=%d;", #sfv, num_cell, (int)config[3]); \
	    printf(" column=%d, n_x=%d;", column, (int)config[13]);	\
	    printf(" line=%d, n_y=%d.\n", line, (int)config[14]);	\
	    exit(2);							\
	}								\
    FV0.sfv = malloc((num_cell + 2) * sizeof(double));			\
    if(FV0.sfv == NULL)							\
	{								\
	    printf("NOT enough memory! %s\n", #sfv);			\
	    exit(5);							\
	}								\
    FV0.sfv[0] = (double)line;						\
    FV0.sfv[1] = (double)column;					\
    if(flu_var_read(fp, FV0.sfv + 2, num_cell))				\
	{								\
	    fclose(fp);							\
	    exit(2);							\
	}								\
    fclose(fp);								\
    } while(0)

/** 
  * @brief      This function reads the 2-D initial data file of velocity/pressure/density.
  * @details    The function initialize the extern pointer FV0.RHO/U/V/P pointing to the
  *             position of a block of memory consisting (line*column+2) variables* of type double.
  *             The value of first of these variables is (line) number;
  *             The value of second of these variables is (column) number;
  *             The following (line*column) variables are the initial value.
  * @param[in]  name: Name of the test example.
  * @return  \b FV0:  Structure of initial data array pointer.
  */
struct flu_var _2D_initialize(const char * name)
{
    struct flu_var FV0;

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
    printf("  delta_y\t= %g\n", config[11]);
    printf("  bondary_x\t= %d\n", (int)config[17]);
    printf("  bondary_y\t= %d\n", (int)config[18]);

    char add[FILENAME_MAX+40]; // The address of the velocity/pressure/density file to read in.
    FILE * fp; // The pointer to the above data files.
    int num_cell, line, column;  // The number of the numbers in the above data files.

     // Open the initial data files and initializes the reading of data.
    STR_FLU_INI(RHO);
    STR_FLU_INI(U);
    STR_FLU_INI(V);
    STR_FLU_INI(P);

    printf("%s data initialized, line = %d, column = %d.\n", name, line, column);
    return FV0;
}


