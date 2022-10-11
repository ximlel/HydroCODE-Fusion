/**
 * @file  sys_pro.c
 * @brief There are some system processing programs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * To realize cross-platform programming.
 * MKDIR:  Create a subdirectory.
 * ACCESS: Determine access permissions for files or folders.
 *       - mode=0: Test for existence.
 *       - mode=2: Test for write permission.
 *       - mode=4: Test for read permission.
 */
#ifdef _WIN32
#include <io.h>
#include <direct.h>
#define ACCESS(path,mode) _access((path),(mode))
#define MKDIR(path)       _mkdir((path))
#elif __linux__
#include <unistd.h>
#include <sys/stat.h>
#define ACCESS(path,mode) access((path),(mode))
#define MKDIR(path)       mkdir((path), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif


/**
 * @brief This function print a progress bar on one line of standard output.
 * @param[in]  pro: Numerator of percent that the process has completed.
 * @param[in]  step: Number of time steps.
 */
void DispPro(const double pro, const int step)
{
	int j;
	static double pro_print = 0.0;
	const double dpro_print = 0.1;
	if (pro >= pro_print)
	    {
		for (j = 0; j < 77; j++)
		    putchar('\b'); // Clears the current line to display the latest progress bar status.
		for (j = 0; j < lround(pro/2); j++)
#ifdef _WIN32
		    putchar('+');  // Print the part of the progress bar that has been completed, denoted by '+'.
#elif __linux__
                printf("\x1b[45m \x1b[0m");
#endif
		for (j = 1; j <= 50-lround(pro/2); j++)
#ifdef _WIN32
		    putchar('-');  // Print how much is left on the progress bar.
#elif __linux__
                printf("\x1b[47m \x1b[0m");
#endif
		fprintf(stdout, "  %6.2f%%   STEP=%-8d", pro, step);  
		fflush(stdout);
		pro_print += dpro_print;
	    }
}

/**
 * @brief This is a function that recursively creates folders.
 * @param[in] pPath: Pointer to the folder Path.
 * @return Folder Creation Status.
 *    @retval -1: The path folder already exists and is readable.
 *    @retval  0: Readable path folders are created recursively.
 *    @retval  1: The path folder is not created properly.
 */
int CreateDir(const char * pPath)
{
	if(0 == ACCESS(pPath,2))
		return -1;

	const char* pCur = pPath;
	char tmpPath[FILENAME_MAX+40];
	memset(tmpPath,0,sizeof(tmpPath));
    
	int pos = 0;
	while(*pCur++!='\0')
		{
			tmpPath[pos++] = *(pCur-1);

			if(*pCur=='/' || *pCur=='\0')
				{
					if(0!=ACCESS(tmpPath,0) && strlen(tmpPath)>0)
						{
							MKDIR(tmpPath);
						}
				}
		}
	if(0 == ACCESS(pPath,2))
		return 0;
	else
		return 1;
}


void init_mem(double * p[], const int n, int ** cell_pt)
{
	for(int k = 0; k < n; ++k)
		{
			p[k] = (double *)calloc(cell_pt[k][0], sizeof(double));
			if(p[k] == NULL)
				{
					printf("Initialize memory fail! DOUBLE data at grid cell points.\n");
					for(int j = 0; j < k; j++)
						{
							free(p[k]);
							p[k] = NULL;
						}
					exit(5);
				}
		}
}


void init_mem_int(int * p[], const int n, int ** cell_pt)
{
	for(int k = 0; k < n; ++k)
		{
			p[k] = (int *)malloc(cell_pt[k][0] * sizeof(int));
			if(p[k] == NULL)
				{
					printf("Initialize memory fail! INT data at grid cell points.\n");
					for(int j = 0; j < k; j++)
						{
							free(p[k]);
							p[k] = NULL;
						}
					exit(5);
				}
		}
}
