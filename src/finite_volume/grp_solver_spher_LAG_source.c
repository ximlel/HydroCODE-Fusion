#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "../include_cii/mem.h"
#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/file_io.h"
#include "../include/meshing.h"
#include "../include/riemann_solver.h"
#include "../include/inter_process.h"


// M=2,
// M=1 planar; M=2 cylindrical ; M=3 spherical
void grp_solver_spher_LAG_source(struct flu_var * FV, struct cell_var_stru * CV, struct radial_mesh_var * smv, double * R[],
				 const int M, const char * problem, double * cpu_time, const int N_plot , double time_plot[])
{
    int i, k=0;
    
    clock_t tic, toc;
    double cpu_time_sum = 0.0;

    //parameters
    double const Timeout = config[1];      // Output time
    double const eps     = config[4];
    int    const N       = (int)config[5]; // the maximum number of time steps
    double const CFL     = config[7];      // CFL condition
    double const dtheta  = config[11];     //initial d_angle
    int    const Ncell   = (int)config[3]; // Number of computing cells in r direction
    int    const Md      = Ncell+2;        // max vector dimension
    double       dt      = config[16];     // the length of the time step

    double U_T,V_T;
    double Rb_NStep,Lb_NStep;
    //double Rb_side[Md],Lb_side[Md],Rbh_side[Md],Lbh_side[Md],Sh[Md];

    double wave_speed[2], dire[4], mid[4];

    double Smax_dr;
    double time_c=0.0;
    int nt = 1;

    struct i_f_var ifv_L = {0}, ifv_R = {0};

    //initial value
    double *DD = CV->RHO[0]; // D:Density;U,V:Velocity;P:Pressure
    double *UU = CV->U[0];
    double *PP = CV->P[0];
    double *EE = CV->E[0];
#ifdef  MULTIFLUID_BASICS
    double *GammaGamma = CV->gamma[0]; // Ratio of special heats
#else
    double *GammaGamma = (double*)ALLOC(Md*sizeof(double)); // Ratio of special heats
    for(i=0; i<Md; i++)//center cell is cell 0
	GammaGamma[i] = config[6];
#endif

    // the slopes of variable values
    double *TmV = (double*)CALLOC(Md, sizeof(double));
    double *DmU = (double*)CALLOC(Md, sizeof(double));
    double *DmD = (double*)CALLOC(Md, sizeof(double));
    double *DmP = (double*)CALLOC(Md, sizeof(double));
    // CV.d_v   = TmV;
    // CV.d_u   = DmU;
    // CV.d_rho = DmD;
    // CV.d_p   = DmP;

    //GRP variables
    double *VLmin = (double*)ALLOC(Md*sizeof(double));
    double *Umin  = (double*)ALLOC(Md*sizeof(double));
    double *Pmin  = (double*)ALLOC(Md*sizeof(double));
    double *DLmin = (double*)ALLOC(Md*sizeof(double));
    double *DRmin = (double*)ALLOC(Md*sizeof(double));
    double *U_t   = (double*)ALLOC(Md*sizeof(double));
    double *P_t   = (double*)ALLOC(Md*sizeof(double));
    double *DL_t  = (double*)ALLOC(Md*sizeof(double));
    double *DR_t  = (double*)ALLOC(Md*sizeof(double));

    double *Rb   = smv->Rb;  //radius and length of outer cell boundary
    double *Lb   = smv->Lb;
    double *RR   = smv->RR;  //centroidal radius and variable in cells
    double *DdrL = smv->DdrL;//distance from boundary to center in a cell
    double *DdrR = smv->DdrR;
    double *Ddr  = smv->Ddr;
    double *Rbh  = smv->Rbh; //h: half time step
    double *Lbh  = smv->Lbh;
    double *dRc  = smv->dRc; //(derivative)centers distance
    double *vol  = smv->vol;

    //flux, conservative variable and wave speed
    double *F_e  = (double*)ALLOC(Md*sizeof(double));
    double *F_u  = (double*)ALLOC(Md*sizeof(double));
    double *F_u2 = (double*)ALLOC(Md*sizeof(double));
    double *mass = (double*)ALLOC(Md*sizeof(double));
    for(i=0; i<=Ncell; i++)//center cell is cell 0
	mass[i] = DD[i]*vol[i];

    for(k=1; k<=N; k++)
	{
	    tic = clock();
	    if (time_c > time_plot[nt] && nt < (N_plot-1))
		{
#ifndef NOTECPLOT
			file_radial_write_TEC(*FV, *smv, problem, time_plot[nt]);
#endif
			nt++;
		}

	    Smax_dr=0.0;

	    GammaGamma[0] = GammaGamma[1];
	    UU[0] = 0.0;
	    DD[0] = DD[1];
	    PP[0] = PP[1];
	    EE[0] = PP[0]/(GammaGamma[0] - 1.0)/DD[0];

	    GammaGamma[Ncell+1] = GammaGamma[Ncell];
	    UU[Ncell+1] = UU[Ncell];
	    DD[Ncell+1] = DD[Ncell];
	    PP[Ncell+1] = PP[Ncell];
	    EE[Ncell+1] = EE[Ncell];

	    VIP_limiter_radial(Ncell, (_Bool)(k-1), DmU, TmV, UU, smv);
	    minmod_limiter_radial(Ncell, (_Bool)(k-1), DmD, DD, smv);
	    minmod_limiter_radial(Ncell, (_Bool)(k-1), DmP, PP, smv);
/*
	    TmV[Ncell+1]=TmV[Ncell];
	    DmU[Ncell+1]=DmU[Ncell];
	    DmD[Ncell+1]=DmD[Ncell];
	    DmP[Ncell+1]=DmP[Ncell];
	    TmV[Ncell+1]=TmV[Ncell];
	    DmU[Ncell+1]=DmU[Ncell];
	    DmD[Ncell+1]=DmD[Ncell];
	    DmP[Ncell+1]=DmP[Ncell];
*/
	    for(i=0; i<=Ncell; i++)
		{
		    ifv_L.gamma = GammaGamma[i];
		    ifv_R.gamma = GammaGamma[i+1];
		    ifv_L.d_rho = DmD[i];
		    ifv_R.d_rho = DmD[i+1];
		    ifv_L.d_p   = DmP[i];
		    ifv_R.d_p   = DmP[i+1];
		    ifv_L.d_u   = DmU[i];
		    ifv_R.d_u   = DmU[i+1];
		    ifv_L.RHO   = DD[i]   + DdrL[i]  *ifv_L.d_rho;
		    ifv_R.RHO   = DD[i+1] - DdrR[i+1]*ifv_R.d_rho;
		    ifv_L.P     = PP[i]   + DdrL[i]  *ifv_L.d_p;
		    ifv_R.P     = PP[i+1] - DdrR[i+1]*ifv_R.d_p;
		    ifv_L.U     = UU[i]   + DdrL[i]  *ifv_L.d_u;
		    ifv_R.U     = UU[i+1] - DdrR[i+1]*ifv_R.d_u;
		    if(ifv_L.P < eps || ifv_R.P < eps || ifv_L.RHO < eps || ifv_R.RHO < eps)
			{
			    printf("<0.0 error on [%d, %d] (t_n, x) - Reconstruction\n", k, i);
			    goto return_NULL;
			}
		    if(!isfinite(ifv_L.d_p)|| !isfinite(ifv_R.d_p)|| !isfinite(ifv_L.d_u)|| !isfinite(ifv_R.d_u)|| !isfinite(ifv_L.d_rho)|| !isfinite(ifv_R.d_rho))
			{
			    printf("NAN or INFinite error on [%d, %d] (t_n, x) - Slope\n", k, i); 
			    goto return_NULL;
			}

		    GRPsolverSLag(wave_speed, dire, mid, &ifv_L, &ifv_R, Rb[i+1], M, eps, eps);

		    if(mid[2] < eps || mid[0] < eps || mid[3] < eps)
			{
			    printf("<0.0 error on [%d, %d] (t_n, x) - STAR\n", k, i);
			    time_c = Timeout;
			}
		    if(!isfinite(mid[1])|| !isfinite(mid[2])|| !isfinite(mid[0])|| !isfinite(mid[3]))
			{
			    printf("NAN or INFinite error on [%d, %d] (t_n, x) - STAR\n", k, i); 
			    time_c = Timeout;
			}
		    if(!isfinite(dire[1])|| !isfinite(dire[2])|| !isfinite(dire[0])|| !isfinite(dire[3]))
			{
			    printf("NAN or INFinite error on [%d, %d] (t_n, x) - DIRE\n", k, i); 
			    time_c = Timeout;
			}

		    Smax_dr=Smax_dr>fabs(wave_speed[0])/Ddr[i] ? Smax_dr:fabs(wave_speed[0])/Ddr[i];
		    Smax_dr=Smax_dr>fabs(wave_speed[1])/Ddr[i+1]?Smax_dr:fabs(wave_speed[1])/Ddr[i+1];

		    Umin[i+1]  = mid[1];
		    Pmin[i+1]  = mid[2];
		    DLmin[i+1] = mid[0];
		    DRmin[i+1] = mid[3];
		    U_t[i+1]   = dire[1];
		    P_t[i+1]   = dire[2];
		    DL_t[i+1]  = dire[0];
		    DR_t[i+1]  = dire[3];
		}

	    dt=CFL/Smax_dr;
	    if(time_c<Timeout&&(time_c+dt)>Timeout)
		dt=Timeout-time_c;//compute for time step

	    for(i=1;i<=Ncell;i++)
		{
		    ifv_L.gamma = GammaGamma[i];
		    ifv_L.RHO   = DD[i]+(0.5*(Rb[i]+Rb[i+1])-RR[i])*ifv_L.d_rho;
		    ifv_L.P     = PP[i]+(0.5*(Rb[i]+Rb[i+1])-RR[i])*ifv_L.d_p;
		    U_T         = UU[i]+(0.5*(Rb[i]+Rb[i+1])-RR[i])*DmU[i];
		    V_T         = 0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta)*TmV[i];
		    ifv_L.U     = -U_T*sin(0.5*dtheta)+V_T*cos(0.5*dtheta);
		    ifv_R       =  ifv_L;
		    ifv_R.U     = -ifv_L.U;
		    ifv_L.d_p   = DmP[i]*cos(0.5*dtheta);
		    ifv_L.d_u   = DmU[i]+TmV[i];

		    AcousticSLagTangent(dire, mid, &ifv_L, &ifv_R, 0.5*(Rb[i]+Rb[i+1])/cos(0.5*dtheta), M, eps);
		    F_u2[i]  = mid[2]+dt*dire[2];
		    VLmin[i] = (ifv_L.U*cos(0.5*dtheta)+ifv_L.V*sin(0.5*dtheta)+dt*dire[1])*sin(0.5*dtheta);
		}

	    for(i=0; i<=Ncell; i++)
		{
		    Umin[i+1] += 0.5*dt*U_t[i+1];
		    Pmin[i+1] += 0.5*dt*P_t[i+1];

		    F_u[i+1] = Pmin[i+1];
		    F_e[i+1] = Pmin[i+1]*Umin[i+1];

		    Rb_NStep=Rb[i+1]+Umin[i+1]*dt;
		    Lb_NStep=2.0*Rb_NStep*tan(0.5*dtheta);
		    Lbh[i+1]=0.5*(Lb[i+1]+Lb_NStep);
		    Rbh[i+1]=(Rb[i+1]*(2.*Lb[i+1]+Lb_NStep)+Rb_NStep*(Lb[i+1]+2.*Lb_NStep))/(3.*(Lb[i+1]+Lb_NStep));
		    Rb[i+1]=Rb_NStep;
		    Lb[i+1]=Lb_NStep;
		    RR[i]=Rb[i+1]-(2.*Lb[i]+Lb[i+1])/(3.*(Lb[i]+Lb[i+1]))*(Rb[i+1]-Rb[i]);
		    /*
		      Rb_side[i]=0.25*(Lb[i]+Lb[i+1])/sin(0.5*dtheta);
		      Lb_side[i]=0.5*(Lb[i+1]-Lb[i])/sin(0.5*dtheta);
		      Rbh_side[i]=0.25*(Lbh[i]+Lbh[i+1])/sin(0.5*dtheta);
		      Lbh_side[i]=0.5*(Lbh[i+1]-Lbh[i])/sin(0.5*dtheta);
		      Sh[i]=Rbh[i+1]*Lbh[i+1]-Rbh[i]*Lbh[i]-sin(dtheta)*Rbh_side[i]*Lbh_side[i];
		    */

		    Umin[i+1]  += 0.5*dt*U_t[i+1];
		    Pmin[i+1]  += 0.5*dt*P_t[i+1];
		    DLmin[i+1] +=     dt*DL_t[i+1];
		    DRmin[i+1] +=     dt*DR_t[i+1];
		}

	    radial_mesh_update(smv);

	    for(i=1;i<=Ncell;i++)//m=2
		{
		    DD[i]=mass[i]/vol[i];
		    UU[i]=UU[i] - dt/mass[i]*((F_u[i+1]-F_u2[i])*Rbh[i+1]*Lbh[i+1]-(F_u[i]-F_u2[i])*Rbh[i]*Lbh[i]);
		    EE[i]=EE[i] - dt/mass[i]*( F_e[i+1]         *Rbh[i+1]*Lbh[i+1]- F_e[i]         *Rbh[i]*Lbh[i]);
		    PP[i]=(EE[i]-0.5*UU[i]*UU[i])*(GammaGamma[i]-1.)*DD[i];
		    if(isnan(PP[i])||isnan(UU[i])||isnan(DD[i]))
			{
			    printf("variable is nan,error!\n");
			    time_c = Timeout;
			}
		    else if (PP[i]<0 || UU[i])
			{
			    printf("p<0,error!\n");
			    time_c = Timeout;
			}
		    DmD[i]=(DD[i]-DD[i-1])/dRc[i];
		    DmU[i]=(UU[i]-UU[i-1])/dRc[i];
		    DmP[i]=(PP[i]-PP[i-1])/dRc[i];
		}
	    UU[0]=0.0;
	    EE[0]=EE[0]-dt/mass[0]*(F_e[1]*Rbh[1]*Lbh[1]);

	    time_c=time_c+dt;
	    if(isfinite(Timeout))
		DispPro(time_c*100.0/Timeout, k);
	    else
		DispPro(k*100.0/N, k);
	    if(time_c > (Timeout - eps) || isinf(time_c))
		break;

	    for(i = 0; i < Md; ++i)
		{
		    CV->RHO[nt-1][i] = CV->RHO[nt][i];
		    CV->U[nt-1][i]   =   CV->U[nt][i];
		    CV->E[nt-1][i]   =   CV->E[nt][i];  
		    CV->P[nt-1][i]   =   CV->P[nt][i];
		}

	    toc = clock();
	    cpu_time[nt] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
	    cpu_time_sum += cpu_time[nt];
	}

    printf("\nTime is up at time step %d.\n", k);
    printf("The cost of CPU time for 1D-GRP Eulerian scheme for this problem is %g seconds.\n", cpu_time_sum);

 return_NULL:
    config[5] = (double)k;
  if(fabs(time_plot[1]) < eps || isinf(time_plot[1]))
      {
	  if(isfinite(time_c))
	      {
		  time_plot[N_plot-2] = time_c - dt;
		  time_plot[N_plot-1] = time_c;
	      }
	  else
	      {
		  time_plot[N_plot-2] = N*dt - dt;
		  time_plot[N_plot-1] = N*dt;
	      }
      }

    DD = NULL;
    UU = NULL;
    PP = NULL;
    EE = NULL;
    GammaGamma = NULL;
    FREE(DmD);
    FREE(DmU);
    FREE(TmV);
    FREE(DmP);
/*
    CV.d_rho = NULL;
    CV.d_u   = NULL;
    CV.d_v   = NULL;
    CV.d_p   = NULL;
*/
    FREE(Umin);
    FREE(VLmin);
    FREE(Pmin);
    FREE(DLmin);
    FREE(DRmin);
    FREE(U_t);
    FREE(P_t);
    FREE(DL_t);
    FREE(DR_t);
    FREE(mass);

    free(F_u);
    free(F_e);
}
