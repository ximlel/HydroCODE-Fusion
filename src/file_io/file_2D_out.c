/**
 * @file  file_2D_out.c
 * @brief This is a set of functions which control the readout of two-dimensional data.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"


/**
 * @brief Print out fluid variable 'v' with array data element 'v_print'.
 */
#define PRINT_NC(v, v_print)						\
    do {								\
    strcpy(file_data, add_out);						\
    strcat(file_data, #v);						\
    strcat(file_data, ".dat");						\
    if((fp_write = fopen(file_data, "w")) == NULL)			\
	{								\
	    printf("Cannot open solution output file: %s!\n", #v);	\
	    exit(1);							\
	}								\
    for(k = 0; k < N; ++k)						\
	{								\
	    for(i = 0; i < n_y; ++i)					\
		{							\
		    for(j = 0; j < n_x; ++j)				\
			fprintf(fp_write, "%.10g\t", (v_print));	\
		    fprintf(fp_write, "\n");				\
		}							\
	    fprintf(fp_write, "\n\n");					\
	}								\
    fclose(fp_write);							\
    } while (0)

/**
 * @brief This function write the 2-D solution into output '.dat' files.
 * @note  It is quite simple so there will be no more comments.
 * @param[in] n_x: The number of x-spatial points in the output data.
 * @param[in] n_y: The number of y-spatial points in the output data.
 * @param[in] N:   The number of time steps in the output data.
 * @param[in] CV:  Structure of variable data in computational grid cells.
 * @param[in] X:   Array of the x-coordinate data.
 * @param[in] Y:   Array of the y-coordinate data.
 * @param[in] cpu_time:  Array of the CPU time recording.
 * @param[in] problem:   Name of the numerical results for the test problem.
 * @param[in] time_plot: Array of the plotting time recording.
 */
void file_2D_write(const int n_x, const int n_y, const int N, const struct cell_var_stru CV[],
		    double ** X, double ** Y, const double * cpu_time, const char * problem, const double time_plot[])
{
    char add_out[FILENAME_MAX+40];
    // Get the address of the output data folder of the test example.
    example_io(problem, add_out, 0);
    
    char file_data[FILENAME_MAX+40];
    FILE * fp_write;

//===================Write Solution File=========================

    int k, i, j;
    PRINT_NC(RHO, CV[k].RHO[j][i]);
    PRINT_NC(U,   CV[k].U[j][i]);
    PRINT_NC(V,   CV[k].V[j][i]);
    PRINT_NC(P,   CV[k].P[j][i]);
    PRINT_NC(E,   CV[k].E[j][i]);
    PRINT_NC(X, 0.25*(X[j][i] + X[j][i+1] + X[j+1][i] + X[j+1][i+1]));
    PRINT_NC(Y, 0.25*(Y[j][i] + Y[j][i+1] + Y[j+1][i] + Y[j+1][i+1]));
    
    strcpy(file_data, add_out);
    strcat(file_data, "time_plot.dat");
    if((fp_write = fopen(file_data, "w")) == NULL)
	{
	    printf("Cannot open solution output file: time_plot!\n");
	    exit(1);
	}
    for(k = 0; k < N; ++k)
	fprintf(fp_write, "%.10g\n", time_plot[k]);
    fclose(fp_write);

    config_write(add_out, cpu_time, problem);
}


/**
 * @brief This function write the 2-D solution into Tecplot output files with point data.
 * @param[in] n_x: The number of x-spatial points in the output data.
 * @param[in] n_y: The number of y-spatial points in the output data.
 * @param[in] N:   The number of time steps in the output data.
 * @param[in] CV:  Structure of variable data in computational grid cells.
 * @param[in] X:   Array of the x-coordinate data.
 * @param[in] Y:   Array of the y-coordinate data.
 * @param[in] cpu_time:  Array of the CPU time recording.
 * @param[in] problem:   Name of the numerical results.
 * @param[in] time_plot: Array of the plotting time recording.
 */
void file_2D_write_POINT_TEC(const int n_x, const int n_y, const int N, const struct cell_var_stru CV[],
			double ** X, double ** Y, const double * cpu_time, const char * problem, const double time_plot[])
{
    double const eps = config[4];
    char add_out[FILENAME_MAX+40];
    // Get the address of the output data folder of the test example.
    example_io(problem, add_out, 0);
    
    char file_data[FILENAME_MAX+40];
    FILE * fp;
    int k, i, j;
    char str_tmp[40];

    //===================Write solution File=========================
    strcpy(file_data, add_out);
    sprintf(str_tmp, "FLU_VAR_%.8g.tec", time_plot[N-1] + eps);
    strcat(file_data, str_tmp);
    if ((fp = fopen(file_data, "w")) == NULL)
	{
	    fprintf(stderr, "Cannot open solution output TECPLOT file of '%s'!\n", problem);
	    exit(1);
	}

    fprintf(fp, "TITLE = \"FE-Volume Point Data\"\n");
    fprintf(fp, "VARIABLES = \"X\", \"Y\"");
    fprintf(fp, ", \"P\", \"RHO\", \"U\", \"V\", \"E\"");
    fprintf(fp, "\n");

    for(k = 0; k < N; ++k)
	{
	    // if (k == N-1)
		// continue;
	    fprintf(fp, "ZONE I=%d, J=%d, SOLUTIONTIME=%.10g, DATAPACKING=POINT\n", n_x, n_y, time_plot[k]);
	    for(i = 0; i < n_y; ++i)
		for(j = 0; j < n_x; ++j)
		    {			    
			fprintf(fp, "%.10g\t", 0.25*(X[j][i] + X[j][i+1] + X[j+1][i] + X[j+1][i+1]));
			fprintf(fp, "%.10g\t", 0.25*(Y[j][i] + Y[j][i+1] + Y[j+1][i] + Y[j+1][i+1]));
			fprintf(fp, "%.10g\t", CV[k].P[j][i]);
			fprintf(fp, "%.10g\t", CV[k].RHO[j][i]);
			fprintf(fp, "%.10g\t", CV[k].U[j][i]);
			fprintf(fp, "%.10g\t", CV[k].V[j][i]);
			fprintf(fp, "%.10g\t", CV[k].E[j][i]);
			fprintf(fp, "\n");
		    }
	    fprintf(fp, "\n");
	}
    fclose(fp);
}
