#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"


void file_spher_write_TEC(const struct flu_var FV, const struct spher_mesh_var smv, const char * problem, const double time)
{
    double const eps    =      config[4];
    int    const Ncell  = (int)config[3];  // Number of computing cells in r direction
    int    const Tcell  = (int)config[14]; // Number of computing cells in \theta direction
    double const dtheta =      config[11];

    char file_data[FILENAME_MAX];
    example_io(problem, file_data, 0);

    FILE * out;
    char str_tmp[40];

    //===================Write solution File=========================
    sprintf(str_tmp, "/FLU_VAR_%.8g.tec", time + eps);
    strcat(file_data, str_tmp);
    if ((out = fopen(file_data, "w")) == NULL)
	{
	    fprintf(stderr, "Cannot open solution output Tecplot file!\n");
	    exit(1);
	}

    fprintf(out, "TITLE = \"Planar Plot of Radially Symmetric Data\"\n");
    fprintf(out, "VARIABLES = X, Y, D, U, P, Gamma\n");
    fprintf(out, "ZONE I=%d, J=%d, F=POINT, SOLUTIONTIME=%.8g\n",Tcell+1,Ncell+1,time);
    for(int j=0; j<=Ncell; j++)
	for(int i=0; i<=Tcell; i++)
	    {
		fprintf(out,"%.10g",smv.RR[i]*cos(j*dtheta));
		fprintf(out,"%.10g",smv.RR[i]*sin(j*dtheta));
		fprintf(out,"%.10g",FV.RHO[i]);
		fprintf(out,"%.10g",FV.U[i]);
		fprintf(out,"%.10g",FV.P[i]);
#ifdef  MULTIFLUID_BASICS
#ifndef MULTIPHASE_BASICS
		fprintf(out,"%.10g",FV.gamma[i]);
#endif
#endif
		fprintf(out,"\n");
	    }
    fclose(out);
}
