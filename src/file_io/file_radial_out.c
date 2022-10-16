#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"


void file_radial_write_TEC(const struct flu_var FV, const double * R, const char * problem, const double time)
{
    double const eps    =      config[4];
    int    const Ncell  = (int)config[3];  // Number of computing cells in r direction
    int    const Tcell  = (int)config[14]; // Number of computing cells in \theta direction
    double const dtheta =      config[11];

    char file_data[FILENAME_MAX];
    example_io(problem, file_data, 0);

    FILE * out;
    int i, j;
    char str_tmp[40];

    //===================Write solution File=========================
    sprintf(str_tmp, "FLU_VAR_%.8g.tec", time + eps);
    strcat(file_data, str_tmp);
    if ((out = fopen(file_data, "w")) == NULL)
	{
	    fprintf(stderr, "Cannot open solution output Tecplot file!\n");
	    exit(1);
	}

    fprintf(out, "TITLE = \"Planar Plot of Radially Symmetric Data\"\n");
    fprintf(out, "VARIABLES = \"X\", \"Y\"");
    fprintf(out, ", \"RHO\", \"U\", \"P\"");
#ifdef  MULTIFLUID_BASICS
#ifndef MULTIPHASE_BASICS
    fprintf(out, ", \"gamma\"");
#endif
#endif
    fprintf(out, "\n");

    fprintf(out, "ZONE I=%d, J=%d, F=POINT, SOLUTIONTIME=%.8g\n", Tcell+1, Ncell+1, time);
    for(i=0; i<=Ncell; i++)
	for(j=0; j<=Tcell; j++)
	    {
		fprintf(out,"%.10g\t",R[i]*cos(j*dtheta));
		fprintf(out,"%.10g\t",R[i]*sin(j*dtheta));
		fprintf(out,"%.10g\t",FV.RHO[i]);
		fprintf(out,"%.10g\t",FV.U[i]);
		fprintf(out,"%.10g\t",FV.P[i]);
#ifdef  MULTIFLUID_BASICS
#ifndef MULTIPHASE_BASICS
		fprintf(out,"%.10g\t",FV.gamma[i]);
#endif
#endif
		fprintf(out,"\n");
	    }
    fclose(out);
}
