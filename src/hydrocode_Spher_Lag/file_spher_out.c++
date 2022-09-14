#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir((path))
#elif __linux__
#include <sys/stat.h>
#define MKDIR(path)  mkdir((path), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif


static void wrin2s(FILE *out,double **R,double **Z,double **D,double **U,double **P,double **Gamma,double t)//for mesh
{
	fprintf(out,"Solution For Euler Equation\n");
	fprintf(out,"variables=Z,R,D,U,P,Gamma\n");
	fprintf(out,"zone I=%d,J=%d,F=POINT,SOLUTIONTIME=%lf\n",Tcell_plot+1,Ncell+1,t);
	int it,jr;
	for(jr=0;jr<=Ncell;jr++)
			for(it=0;it<=Tcell_plot;it++)
				fprintf(out,"%lf %lf %lf %lf %lf %lf\n",Z[jr][it],R[jr][it],D[jr][it],U[jr][it],P[jr][it],Gamma[jr][it]);
}


static void Write(FILE *out,double *XX,int N)// write data in file
{
	int i;
	for(i=0;i<N;i++)
		{
			fprintf(out,"%lf ",XX[i]);
		}
}


void file_spher_write_TEC(const int n_x, const int n_y, const int N, const struct cell_var_stru CV[],
			double ** X, double ** Y, const double * cpu_time, const char * problem, const double time_plot[])
{
	outs=fopen("../datas_fin.m","w");
	fprintf(outs,"RR=[");
	Write(outs,RR,Ncell);
	fprintf(outs,"];\n");
	fprintf(outs,"DD=[");
	Write(outs,DD,Ncell);
	fprintf(outs,"];\n");
	fprintf(outs,"UU=[");
	Write(outs,UU,Ncell);
	fprintf(outs,"];\n");
	fprintf(outs,"PP=[");
	Write(outs,PP,Ncell);
	fprintf(outs,"];\n");
	fclose(outs);

	char file_data[FILENAME_MAX];
	FILE *out,*outs;
	MKDIR(DATAOUT);

	for(i=0;i<=Ncell;i++)
		for(j=0;j<=Tcell;j++)
			{
				rb[i][j]=Rb[i]/cos(0.5*dtheta)*sin(j*dtheta);
				zb[i][j]=Rb[i]/cos(0.5*dtheta)*cos(j*dtheta);
			}
	for(i=1;i<=Ncell;i++)
		for(j=0;j<=Tcell;j++)
			{
				DD2[i][j]=0.5*(DD[i-1]+DD[i]);
				UUxi2[i][j]=0.5*(UU[i-1]+UU[i]);
				PP2[i][j]=0.5*(PP[i-1]+PP[i]);
				GammaGamma2[i][j]=0.5*(GammaGamma[i-1]+GammaGamma[i]); //GammaGamma[i];
			}
	for(j=0;j<=Tcell;j++)
		{
			DD2[0][j]=DD[0];
			UUxi2[0][j]=UU[0];
			PP2[0][j]=PP[0];
			GammaGamma2[0][j]=GammaGamma[0];
		}
	sprintf(file_data, "%s/FLU_VAR_%.5g.dat", "../data_out", time);
	out=fopen(file_data,"w");
	wrin2s(out,rb,zb,DD2,UUxi2,PP2,GammaGamma2,time);
	fclose(out);
}
