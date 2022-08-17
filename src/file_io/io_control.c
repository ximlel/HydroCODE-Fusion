#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <stdbool.h>

#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/file_io.h"

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

/** @brief This function produces folder path for data input or output.
 *  @param[in]  example:   Name of the test example.
 *  @param[out] add_mkdir: Folder path for data input or output.
 *  @param[in]  i_or_o:    Conversion parameters for data input/output.
 *                         - 0:             data output.
 *                         - else (e.g. 1): data input.
 */
void example_io(const char *example, char *add_mkdir, const int i_or_o)
{
	const int dim = (int)config[0];
	const int el = (int)config[8];
	const int order = (int)config[9];

	char *str_tmp, str_order[10];
	switch (dim)
	    {
	    case 1 :
		str_tmp = "one-dim/";   break;
	    case 2 :
		str_tmp = "two-dim/";   break;
	    case 3 :
		str_tmp = "three-dim/"; break;
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
				str_tmp = "EUL_"; break;
			    case 1 :
				str_tmp = "LAG_"; break;
			    case 2 :
				str_tmp = "ALE_"; break;
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
		if(CreateDir(add_mkdir))
		    {
			fprintf(stderr, "Output directory '%s' construction failed.\n", add_mkdir);
			exit(1);
		    }
		else
		    printf("Output directory '%s' constructed.\n", add_mkdir);				
	    }
	else if (ACCESS(add_mkdir,4) == -1)
	    {
		fprintf(stderr, "Input directory is unreadable!\n");
		exit(1);
	    }

	strcat(add_mkdir, "/");
}


/**
 * @brief      This function counts how many numbers are there in the initial data file. 
 * @param[in]  fp:  The pointer to the input file.
 * @param[in]  add: The address of the input file.
 * @return     The number of the numbers in the initial data file.
 *    @retval  -1: If the given number of column is not coincided with that in the data file.
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
