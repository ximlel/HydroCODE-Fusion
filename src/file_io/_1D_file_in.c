/**
 * @file  _1D_file_in.c
 * @brief This is a set of functions which control the read-in of one-dimensional data.
 */

#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"


/**
 * @brief This function reads the 1D initial data file to generate the initial data.
 * @param[in]  fp: The pointer to the input file.
 * @param[out]  U: The pointer to the data array of fluid variables.
 * @param[in] num: The number of the numbers in the input file. 
 * @return It returns 0 if successfully read the file,
 *         while returns the index of the wrong entry.
 */
static int _1D_flu_var_read(FILE * fp, double * U, const int num)
{
  int idx = 0, j = 0; // j is a frequently used index for spatial variables.
  char number[100]; // A string that stores a number.
  char ch, *endptr;
  // int sign = 1;
    
  while((ch = getc(fp)) != EOF)
  {
    if(isspace(ch) && idx)
    {
      number[idx] = '\0';
      idx = 0;
      // format_string() and str2num() in 'str_num_common.c' are deprecated.
      /*
      sign = format_string(number);
      if(!sign)
	return j+1;
      else if(j == num)
	return j;
      U[j] = sign * str2num(number);
      */
      errno = 0;
      U[j] = strtod(number, &endptr);
      if (errno == ERANGE || *endptr != '\0')
	  {
	      printf("The %dth entry in the initial data file is not a double-precision floats.\n", j+1);
	      return j+1;
	  }
      else if(j == num)
	  {
	      printf("Error on the initial data file reading!\n");
	      return j;
	  }
      ++j;
    }
    else if((ch == 46) || (ch == 45) || (ch == 69) || (ch == 101) || isdigit(ch))
      number[idx++] = ch;
  }
  return 0;
}


/** 
  * @brief      This function reads the 1D initial data file of velocity/pressure/density.
  * @details    The function initialize the extern pointer FV0.RHO/U/P pointing to the
  *             position of a block of memory consisting (m+1) variables* of type double.
  *             The value of first of these variables is m.
  *             The following m variables are the initial value.
  * @param[in]  name: Name of the test example.
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
	if(isinf(config[3]))						\
	    config[3] = (double)num_cell;				\
	else if(num_cell != (int)config[3])				\
	    {								\
		printf("Input unequal! num_%s=%d,", #sfv, num_cell);	\
		printf(" num_cell=%d.\n", num_cell, (int)config[3]);	\
		exit(2);						\
	    }								\
	FV0.sfv = malloc((num_cell + 1) * sizeof(double));		\
	FV0.sfv[0] = (double)num_cell;					\
	if(_1D_flu_var_read(fp, FV0.sfv+1, num_cell))			\
	    {								\
		fclose(fp);						\
		exit(2);						\
	    }								\
	fclose(fp);							\
    } while(0)

void _1D_initialize(const char * name)
{
    char add_in[FILENAME_MAX+40]; 
    // Get the address of the initial data folder of the test example.
    example_io(name, add_in, 1);
    
    /* 
     * Read the configuration data.
     * The detail could be seen in the definition of array config
     * referring to file 'doc/config.csv'.
     */
    _1D_configurate(add_in);
    printf("  delta_x\t= %g\n", config[10]);
    printf("  delta_t\t= %g\n", config[16]);
    printf("  bondary\t= %d\n", (int)config[17]);
  
    char add[FILENAME_MAX+40]; // The address of the velocity/pressure/density file to read in.
    FILE * fp; // The pointer to the above data files.
    int num_cell;  // The number of the numbers in the above data files.
    
    // Open the initial data files and initializes the reading of data.
    STR_FLU_INI(RHO);
    STR_FLU_INI(U);
    STR_FLU_INI(P);

    printf("%s data initialized, m=%d.\n", name, num_cell);
}
