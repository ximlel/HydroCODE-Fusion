#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"


#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define PRINT_NP(v)									\
	do {											\
		if (mv.v == NULL)							\
			for(k = 0; k < mv.num_pt; k++)			\
				fprintf(fp, "%.15g\n", 0.0);		\
		else										\
			for(k = 0; k < mv.num_pt; k++)			\
				fprintf(fp, "%.15g\n", mv.v[k]);	\
		fprintf(fp,"\n");							\
	} while (0)

#define PRINT_NC(v)									\
	do {											\
		if (FV.v == NULL)							\
			for(k = 0; k < num_cell; k++)			\
				fprintf(fp, "%.15g\n", 0.0);		\
		else										\
			for(k = 0; k < num_cell; k++)			\
				fprintf(fp, "%.15g\n", FV.v[k]);	\
		fprintf(fp,"\n");							\
	} while (0)

void file_write_TEC(const struct flu_var FV, const struct mesh_var mv, const char * problem, const double time, const int dim_plot)
{
	int k;
	
	const int num_cell = (int)config[3];
	const int dim = (int)config[0];

	int cell_type = 0;
	for (k = 0; k < num_cell; k++)
		cell_type= MAX(mv.cell_pt[0][0], cell_type);

	int num_data;
	
	FILE * fp;
  
	char file_data[FILENAME_MAX];	
	example_io(problem, file_data, 0);	

	char tmp[50];
	sprintf(tmp, "/FLU_VAR_%.15g.tec", time);
	strcat(file_data, tmp);

	//===================Write solution File=========================
	
	if ((fp = fopen(file_data, "w")) == NULL)
		{
			fprintf(stderr, "Cannot open solution output file!\n");
			exit(1);
		}
  
	fprintf(fp, "TITLE = \"FE-Volume Brick Data\"\n");
	fprintf(fp, "VARIABLES = \"X\"");
	if (dim_plot > 1)
		fprintf(fp, ", \"Y\"");
	if (dim_plot > 2)
		fprintf(fp, ", \"Z\"");
	fprintf(fp, ", \"P\", \"RHO\", \"U\"");
	if (dim_plot > 1)
		fprintf(fp, ", \"V\"");
	if (dim_plot > 2)
		fprintf(fp, ", \"W\"");
	if ((int)config[2] == 2)		 							
		fprintf(fp, ", \"PHI\"");
	fprintf(fp, "\n");
					
	fprintf(fp, "ZONE T=\"Fluid Region\", SOLUTIONTIME=%.15g \n", time);
	fprintf(fp, "NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, ", mv.num_pt, num_cell);
	
	if (cell_type < 2)
		{
			printf("NON ZONETYPE!");
			fclose(fp);
			remove(file_data);
			exit(2);
		}
	else if (cell_type <= 2 && dim == 1)
		fprintf(fp, "ZONETYPE=FELINESEG\n");
	else if (cell_type <= 3 && dim == 2)
		fprintf(fp, "ZONETYPE=FETRIANGLE\n");
	else if (cell_type <= 4 && dim == 2)
		fprintf(fp, "ZONETYPE=FEQUADRILATERAL\n");
	else if (cell_type <= 4 && dim == 3)
		fprintf(fp, "ZONETYPE=FETETRAHEDRON\n");
	else if (cell_type <= 8 && dim == 3)
		fprintf(fp, "ZONETYPE=FEBRICK\n");
	else
		{
			printf("NON ZONETYPE!");
			fclose(fp);
			remove(file_data);
			exit(2);
		}
	
	num_data = dim_plot + 2 + dim_plot;
	if ((int)config[2] == 2)
		num_data++;
	fprintf(fp, "VARLOCATION=([%d-%d]=CELLCENTERED)\n", dim_plot + 1, num_data);

	PRINT_NP(X);
	if (dim_plot > 1)
		PRINT_NP(Y);
	if (dim_plot > 2)
		PRINT_NP(Z);
	PRINT_NC(P);
	PRINT_NC(RHO);
	PRINT_NC(U);
	if (dim_plot > 1)
		PRINT_NC(V);
	if (dim_plot > 2)
		PRINT_NC(W);
	if ((int)config[2] == 2)
		PRINT_NC(PHI);									
	
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


#define PRINT_SCA(v)							\
	do {										\
		fprintf(fp, "SCALARS " #v " double\n");	\
		fprintf(fp, "LOOKUP_TABLE default\n");	\
		for(k = 0; k < num_cell; k++)			\
			fprintf(fp, "\t%.15g", FV.v[k]);	\
		fprintf(fp, "\n");						\
	} while (0)

void file_write_VTK_3D(const struct flu_var FV, const struct mesh_var mv, const char * problem)
{
	int k;
	
	const int num_cell = (int)config[3];
	
	FILE * fp;
  
	char file_data[FILENAME_MAX];	
	example_io(problem, file_data, 0);	

	strcat(file_data, "/FLU_VAR.vtk");

	//===================Write solution File=========================
	
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
			fprintf(fp, "\t%.15g", mv.X[k]);
			fprintf(fp, "\t%.15g", mv.Y[k]);
			if (mv.Z == NULL)
				fprintf(fp, "\t%.15g\n", 0.0);
			else
				fprintf(fp, "\t%.15g\n", mv.Z[k]);
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
	if (mv.cell_type == NULL)
		for(k = 0; k < num_cell; k++)
			fprintf(fp, "\t7\n");
	else
		for(k = 0; k < num_cell; k++)
			fprintf(fp, "\t7\n");
//			fprintf(fp, "\t%d\n", mv.cell_type[k]);	
	fprintf(fp, "\n");

	fprintf(fp, "CELL_DATA %d\n",num_cell);
	PRINT_SCA(P);
	PRINT_SCA(RHO);
	if ((int)config[2] == 2)
		PRINT_SCA(PHI);

	fprintf(fp, "VECTORS velocity double\n");
	
			
	for(k = 0; k < num_cell; k++)
		{			
			if (FV.W == NULL)
				fprintf(fp, "\t%.15g %.15g %.15g", FV.U[k], FV.V[k], 0.0);
			else
				fprintf(fp, "\t%.15g %.15g %.15g", FV.U[k], FV.V[k], FV.W[k]);	
		}
	fprintf(fp, "\n");	

	fclose(fp);
}
