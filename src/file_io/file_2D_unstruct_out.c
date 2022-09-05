/**
 * @file  file_2D_unstruct_out.c
 * @brief This is a set of functions which control the readout of two-dimensional data on unstructured grid.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/file_io.h"


/**
 * @brief Print out 'v'(X/Y) coordinates of spatial points.
 */
#define PRINT_NP(v)					\
    do {						\
	if (mv.v == NULL)				\
	    for(k = 0; k < mv.num_pt; k++)		\
		fprintf(fp, "%.10g\n", 0.0);		\
	else						\
	    for(k = 0; k < mv.num_pt; k++)		\
		fprintf(fp, "%.10g\n", mv.v[k]);	\
	fprintf(fp,"\n");				\
    } while (0)

/**
 * @brief Print out values of fluid variable 'v' in computational grid cells.
 */
#define PRINT_NC(v)					\
    do {						\
	if (FV.v == NULL)				\
	    for(k = 0; k < num_cell; k++)		\
		fprintf(fp, "%.10g\n", 0.0);		\
	else						\
	    for(k = 0; k < num_cell; k++)		\
		fprintf(fp, "%.10g\n", FV.v[k]);	\
	fprintf(fp,"\n");				\
	} while (0)

/**
 * @brief This function write the 2-D solution into Tecplot output files with block data.
 * @param[in] FV: Structure of fluid variable data array in computational grid.
 * @param[in] mv: Structure of meshing variable data.
 * @param[in] problem:   Name of the numerical results for the test problem.
 * @param[in] time_plot: Array of the plotting time recording.
 */
void file_write_2D_BLOCK_TEC(const struct flu_var FV, const struct mesh_var mv, const char * problem, const double time)
{
    const double eps = config[4];
    const int num_cell = (int)config[3];

    int k, num_data;
    int cell_type = 0;
    for (k = 0; k < num_cell; k++)
	cell_type = MAX(mv.cell_pt[0][0], cell_type);
  
    char file_data[FILENAME_MAX];	
    example_io(problem, file_data, 0);

    FILE * fp;
    char str_tmp[40];

	//===================Write solution File=========================
	sprintf(str_tmp, "/FLU_VAR_%.8g.tec", time + eps);
	strcat(file_data, str_tmp);
	if ((fp = fopen(file_data, "w")) == NULL)
		{
			fprintf(stderr, "Cannot open solution output Tecplot file!\n");
			exit(1);
		}
  
	fprintf(fp, "TITLE = \"FE-Volume Brick Data\"\n");
	fprintf(fp, "VARIABLES = \"X\", \"Y\"");
	num_data  = 2;
	fprintf(fp, ", \"P\", \"RHO\", \"U\", \"V\"");
	num_data += 4;
#ifdef MULTIFLUID_BASICS
	fprintf(fp, ", \"Z_a\"");
	num_data += 1;
#ifdef MULTIPHASE_BASICS
	fprintf(fp, ", \"P_b\", \"RHO_b\", \"U_b\", \"V_b\"");
	num_data += 4;
#else
	fprintf(fp, ", \"PHI\", \"gamma\"");
	num_data += 2;
#endif
#endif
	fprintf(fp, "\n");
					
	fprintf(fp, "ZONE T=\"Fluid Region\", SOLUTIONTIME=%.8g \n", time + eps);
	fprintf(fp, "NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, ", mv.num_pt, num_cell);
	
	if (cell_type == 3)
		fprintf(fp, "ZONETYPE=FETRIANGLE\n");
	else if (cell_type == 4)
		fprintf(fp, "ZONETYPE=FEQUADRILATERAL\n");
	else
		{
			printf("NON ZONETYPE!");
			fclose(fp);
			remove(file_data);
			exit(2);
		}

	fprintf(fp, "VARLOCATION=([%d-%d]=CELLCENTERED)\n", 3, num_data);

	PRINT_NP(X);
	PRINT_NP(Y);
	PRINT_NC(P);
	PRINT_NC(RHO);
	PRINT_NC(U);
	PRINT_NC(V);
#ifdef MULTIFLUID_BASICS
	PRINT_NC(Z_a);
#ifdef MULTIPHASE_BASICS
	PRINT_NC(P_b);
	PRINT_NC(RHO_b);
	PRINT_NC(U_b);
	PRINT_NC(V_b);
#else
	PRINT_NC(PHI);
	PRINT_NC(gamma);
#endif
#endif
	
	for(k = 0; k < num_cell; k++)
		{
			for(int i = 1; i <= cell_type; i++)
				{
					if (i <= mv.cell_pt[k][0])
						fprintf(fp, "\t%d", mv.cell_pt[k][i]+1);
					else
						fprintf(fp, "\t%d", mv.cell_pt[k][mv.cell_pt[k][0]]+1);
				}
			fprintf(fp, "\n");
		}

	fclose(fp);
}


/**
 * @brief Print out values of scalar fluid variable 'v' in computational grid cells.
 */
#define PRINT_SCA(v)				\
    do {					\
	fprintf(fp, "SCALARS " #v " double\n");	\
	fprintf(fp, "LOOKUP_TABLE default\n");	\
	for(k = 0; k < num_cell; k++)		\
	    fprintf(fp, "\t%.10g", FV.v[k]);	\
	fprintf(fp, "\n");			\
    } while (0)

/**
 * @brief This function write the 2-D solution into VTK 3D output files with brick data.
 * @param[in] FV: Structure of fluid variable data array in computational grid.
 * @param[in] mv: Structure of meshing variable data.
 * @param[in] problem:   Name of the numerical results for the test problem.
 * @param[in] time_plot: Array of the plotting time recording.
 */
void file_write_3D_VTK(const struct flu_var FV, const struct mesh_var mv, const char * problem, const double time)
{
    const double eps = config[4];
	const int num_cell = (int)config[3];

	char file_data[FILENAME_MAX];	
	example_io(problem, file_data, 0);

	FILE * fp;
    char str_tmp[40];
	int k;

	//===================Write solution File=========================
	sprintf(str_tmp, "/FLU_VAR_%.8g.vtk", time + eps);
	strcat(file_data, str_tmp);
	if ((fp = fopen(file_data, "w")) == NULL)
		{
			fprintf(stderr, "Cannot open solution output file!\n");
			exit(1);
		}

	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "FE-Volume Brick Data\n");
	fprintf(fp, "ASCII\n\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d double\n", mv.num_pt);
	for(k = 0; k < mv.num_pt; k++)
		{
			fprintf(fp, "\t%.10g", mv.X[k]);
			fprintf(fp, "\t%.10g", mv.Y[k]);
			fprintf(fp, "\t%.10g\n", 0.0);
		}
    fprintf(fp, "\n");

	int size = 0;
	for(k = 0; k < num_cell; k++)
		size = size + mv.cell_pt[k][0]+1;	
    fprintf(fp, "CELLS %d %d\n", num_cell, size);
	for(k = 0; k < num_cell; k++)
		{
			for(int i = 0; i <= mv.cell_pt[k][0]; i++)
				fprintf(fp, "\t%d", mv.cell_pt[k][i]);
			fprintf(fp, "\n");
		}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n",num_cell);
	for(k = 0; k < num_cell; k++)
		fprintf(fp, "\t7\n");
	fprintf(fp, "\n");

	fprintf(fp, "CELL_DATA %d\n",num_cell);
	PRINT_SCA(P);
	PRINT_SCA(RHO);
#ifdef MULTIFLUID_BASICS
	PRINT_SCA(Z_a);
#ifdef MULTIPHASE_BASICS
	PRINT_SCA(P_b);
	PRINT_SCA(RHO_b);
	PRINT_SCA(U_b);
	PRINT_SCA(V_b);
#else
	PRINT_SCA(PHI);
	PRINT_SCA(gamma);
#endif
#endif

	fprintf(fp, "VECTORS velocity double\n");
	for(k = 0; k < num_cell; k++)
		fprintf(fp, "\t%.10g %.10g %.10g", FV.U[k], FV.V[k], 0.0);
	fprintf(fp, "\n");	

	fclose(fp);
}
