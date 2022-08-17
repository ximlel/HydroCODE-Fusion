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
static int _1D_flu_var_read(FILE * fp, double * U, const int num)
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
    printf("delta_x    = %g\n", config[10]);
    printf("delta_t    = %g\n", config[16]);
    printf("bondary    = %d\n", (int)config[17]);

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
  file_read_state = _1D_flu_var_read(fp_rho, RHO0+1, num_rho);
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
  file_read_state = _1D_flu_var_read(fp_U, U0+1, num_U);
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
  file_read_state = _1D_flu_var_read(fp_P, P0+1, num_P);
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
