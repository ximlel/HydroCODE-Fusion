/**
 * @file  io_control.c
 * @brief This is a set of common functions which control the input/output data.
 */

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <ctype.h>

#include "../include/var_struc.h"
#include "../include/tools.h"

/*
 * To realize cross-platform programming.
 * ACCESS: Determine access permissions for files or folders.
 *       - mode=0: Test for existence.
 *       - mode=2: Test for write permission.
 *       - mode=4: Test for read permission.
 */
#ifdef _WIN32
#include <io.h>
#define ACCESS(path,mode) _access((path),(mode))
#elif __linux__
#include <unistd.h>
#define ACCESS(path,mode) access((path),(mode))
#endif


/** @brief This function produces folder path for data input or output.
 *  @param[in]  example:   Name of the test example/numerical results.
 *  @param[out] add_mkdir: Folder path for data input or output.
 *  @param[in]  i_or_o:    Conversion parameters for data input/output.
 *    @arg 0:              data output.
 *    @arg 1(non-0 value): data input.
 */
void example_io(const char *example, char *add_mkdir, const int i_or_o)
{
	const int dim   = (int)config[0];
	const int el    = (int)config[8];
	const int order = (int)config[9];

	static int output_const = 0;
	char str_tmp[11], str_order[11];
	switch (dim)
	    {
	    case 1 :
		strcpy(str_tmp, "one-dim/");   break;
	    case 2 :
		strcpy(str_tmp, "two-dim/");   break;
	    case 3 :
		strcpy(str_tmp, "three-dim/"); break;
	    default :
		fprintf(stderr, "Strange computational dimension!\n");
		exit(2);
	    }
	if (i_or_o == 0) // Output
		{
			strcpy(add_mkdir, "../../data_out/");
			strcat(add_mkdir, str_tmp);
			switch (el)
			    {
			    case 0 :
				strcpy(str_tmp, "EUL_"); break;
			    case 1 :
				strcpy(str_tmp, "LAG_"); break;
			    case 2 :
				strcpy(str_tmp, "ALE_"); break;
			    default :
				fprintf(stderr, "Strange description method of fluid motion!\n");
				exit(2);
			    }
			strcat(add_mkdir, str_tmp);
			sprintf(str_order, "%d_order/", order);
			strcat(add_mkdir, str_order);
		}
	else // Input
		{
			strcpy(add_mkdir, "../../data_in/");
			strcat(add_mkdir, str_tmp);
		}
	strcat(add_mkdir, example);

	if (i_or_o == 0)
	    {
		if(CreateDir(add_mkdir) == 1)
		    {
#ifdef _WIN32
			fprintf(stderr, "Output directory '%s' construction failed!\n", add_mkdir);
#elif __linux__
			fprintf(stderr, "\x1b[47;34mOutput directory '%s' construction failed!\x1b[0m\n", add_mkdir);
#endif
			exit(1);
		    }
		else if(output_const == 0)
		    {
#ifdef _WIN32
			printf("Output directory '%s' has been constructed.\n", add_mkdir);				
#elif __linux__
			printf("\x1b[47;34mOutput directory '%s' has been constructed.\x1b[0m\n", add_mkdir);
#endif
			output_const = 1;
		    }
	    }
	else if (ACCESS(add_mkdir,4) == -1)
	    {
		fprintf(stderr, "Input directory '%s' is nonexistent or unreadable!\n", add_mkdir);
		exit(1);
	    }

	strcat(add_mkdir, "/");
}


/**
 * @brief      This function counts how many numbers are there in the initial data file. 
 * @param[in]  fp:  The pointer to the input file.
 * @param[in]  add: The address of the input file.
 * @return \b  num: The number of the numbers in the initial data file.
 */
int flu_var_count(FILE * fp, const char * add)
{
    int num = 0;  // Data number.
    /* We read characters one by one from the data file.
     * "flg" helps us to count.
     * -# 1: when read a number-using character (0, 1, 2, ..., e, E, minus sign and dot).
     * -# 0: when read a non-number-using character. 
     */
    int flg = 0; 
    int ch;

    while((ch = getc(fp)) != EOF) // Count the data number.
	{
	    if (ch == 45 || ch == 46 || ch == 69 || ch == 101 || isdigit(ch))
		flg = 1;
	    else if (!isspace(ch))
		{
		    fprintf(stderr, "Input contains illegal character(ASCII=%d, flag=%d) in the file '%s'!\n", ch, flg, add);
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


/**
 * @brief      This function counts the line and column number of the numbers are there in the initial data file. 
 * @param[in]  fp:  The pointer to the input file.
 * @param[in]  add: The address of the input file.
 * @param[out] n_x: The colume number of the numbers in the initial data file.
 * @return \b line: The line number of the numbers in the initial data file.
 */
int flu_var_count_line(FILE * fp, const char * add, int * n_x)
{
  int line = 0, column = 0;
  /* We read characters one by one from the data file.
   * "flg" helps us to count.
   * -# 1: when read a number-using character (0, 1, 2, ..., e, E, minus sign and dot).
   * -# 0: when read a non-number-using character. 
   */
  int flag = 0;
  int ch;

  do { // Count the data line number.
      ch = getc(fp);
      if(ch == '\n' || ch == EOF)
	  {
	      if(flag)
		  ++column;
	      flag = 0;
	      if(column)
		  {
		      if(!line)
			  *n_x = column;
		      else if(column != *n_x)
			  {
			      printf("Error in input data file '%s', line=%d, column=%d, n_x=%d\n", add, line, column, *n_x);
			      return 0;
			  }
		      ++line;
		      column = 0;
		  }
	  }      
      else if(ch == 45 || ch == 46 || ch == 69 || ch == 101 || isdigit(ch))
	  flag = 1;
      else if (!isspace(ch))
	  {
	      printf("Input contains illigal character(ASCII=%d, flag=%d) in the file '%s', line=%d!\n", ch, flag, add, line);
	      return 0;
	  }
      else if(flag)
	  {
	      ++column;
	      flag = 0;
	  }
  } while(ch != EOF);

  rewind(fp);
  return line;
}


/**
 * @brief This function reads the initial data file to generate the initial data.
 * @param[in]  fp: The pointer to the input file.
 * @param[out]  U: The pointer to the data array of fluid variables.
 * @param[in] num: The number of the numbers in the input file. 
 * @return  It returns 0 if successfully read the file,
 *          while returns the index of the wrong entry.
 */
int flu_var_read(FILE * fp, double * U, const int num)
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
 * @brief Compare function of double for sort function 'qsort()'.
 */
static int compare_double(const void *a,const void *b)
{
    double ret = *(double*)a - *(double*)b;
    if(ret > 0.0)
	return 1;
    else if(ret < 0.0)
	return -1;
    else
	return 0;
}


/**
 * @brief This function reads the time data file for plotting 'time_plot.dat' and 
 *        initialize tha array 'time_plot[]'.
 * @param[in]  add_in:    Adress of the initial data folder of the test example.
 * @param[in]  N_max:     The maximum number of data dimension storing fluid variables in memory.
 * @param[out] N_plot:    Pointer to the number of time steps for plotting.
 * @param[out] time_plot: Pointer to the array of the plotting time recording.
 * @return  It returns the proper number of data dimension storing fluid variables in memory.
 */
int time_plot_read(const char * add_in, const int N_max, int * N_plot, double * time_plot[])
{
    _Bool r = true; // r: Whether to read data file successfully.
    FILE * fp;
    char add[FILENAME_MAX+40];
    strcpy(add, add_in);
    strcat(add, "time_plot.txt");
    // Open the time data file for plotting.
    if((fp = fopen(add, "r")) == NULL)
	{
	    strcpy(add, add_in);
	    strcat(add, "time_plot.dat");
	}
    if((fp = fopen(add, "r")) == NULL)
	{
	    printf("No time data file for plotting! Only the initial data and final result will be plotted.\n");
	    *N_plot = 2;
	    r = false;
	}
    else
	*N_plot = flu_var_count(fp, add) + 2;
    if (*N_plot < 2)
	{
	    printf("Error in counting time data file for plotting!\n");
	    fclose(fp);
	    exit(2);
	}
    *time_plot = (double*)malloc((*N_plot)*sizeof(double));
    if(*time_plot == NULL)
	{
	    printf("NOT enough memory! time_plot[]\n");
	    exit(5);
	}
    (*time_plot)[0] = 0.0;
    (*time_plot)[*N_plot - 1] = config[1];
    if(r)
	{
	    if(flu_var_read(fp, *time_plot + 1, *N_plot - 2))
		{
		    fclose(fp);
		    exit(2);
		}
	    printf("Load time data file for plotting! Plot time step is %d.\n", *N_plot - 2);
	    qsort(*time_plot, *N_plot-1, sizeof(double), compare_double);
	    fclose(fp);
	}
    return N_max<(*N_plot) ? N_max : (*N_plot);
}
