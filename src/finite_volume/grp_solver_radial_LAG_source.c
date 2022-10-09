#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#include "../include_cii/mem.h"
#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/file_io.h"
#include "../include/meshing.h"
#include "../include/riemann_solver.h"
#include "../include/inter_process.h"


// M=1 planar; M=2 cylindrical √; M=3 spherical
void grp_solver_radial_LAG_source(struct cell_var_stru CV, struct radial_mesh_var * smv, double * R[],
				 const int M, const char * problem, double * cpu_time, int * N_plot , double time_plot[])
{
    int i, k=0;
    
    clock_t tic, toc;
    double cpu_time_sum = 0.0;

    //parameters
    double const Timeout = config[1];       // Output time
    double const eps     = config[4];
    int    const N       = (int)config[5];  // the maximum number of time steps
    double const CFL     = config[7];       // CFL condition
    double const dtheta  = config[11];      //initial d_angle
    int    const Ncell   = (int)config[3];  // Number of computing cells in r direction
    int    const Md      = Ncell+2;         // max vector dimension
    double       dt      = config[16];      // the length of the time step

    double U_T,V_T;
    double Rb_NStep,Lb_NStep;
    //double Rb_side[Md],Lb_side[Md],Rbh_side[Md],Lbh_side[Md],Sh[Md];

    double wave_speed[2], dire[4], mid[4];

    double Smax_dr;
    double time_c = 0.0;
    _Bool stop_t = false;
    int nt = 0;

    struct i_f_var ifv_L = {0}, ifv_R = {0};
    struct flu_var FV;

    // initial value
    double *DD = CV.RHO[0]; // D:Density;U,V:Velocity;P:Pressure
    double *UU = CV.U[0];
    double *PP = CV.P[0];
    double *EE = CV.E[0];
#ifdef MULTIFLUID_BASICS
    double *GammaGamma = CV.gamma[0]; // Ratio of special heats
    GammaGamma[0] = GammaGamma[1];
    GammaGamma[Ncell+1] = GammaGamma[Ncell];
#else
    double *GammaGamma = (double*)ALLOC(Md*sizeof(double)); // Ratio of special heats
    for(i = 0; i < Md; i++) //center cell is cell 0
	GammaGamma[i] = config[6];
#endif
    UU[0] = 0.0;
    DD[0] = DD[1];
    PP[0] = PP[1];
    EE[0] = PP[0]/(GammaGamma[0] - 1.0)/DD[0];

    // the slopes of variable values
    double *TmV = (double*)CALLOC(Md, sizeof(double));
    double *DmU = (double*)CALLOC(Md, sizeof(double));
    double *DmD = (double*)CALLOC(Md, sizeof(double));
    double *DmP = (double*)CALLOC(Md, sizeof(double));

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
    double *dRc  = smv->dRc; //(derivative)centers distance
    double *vol  = smv->vol;
    double *Rbh  = (double*)ALLOC(Md*sizeof(double)); //h: half time step
    double *Lbh  = (double*)ALLOC(Md*sizeof(double));

    //flux, conservative variable and wave speed
    double *F_e  = (double*)ALLOC(Md*sizeof(double));
    double *F_u  = (double*)ALLOC(Md*sizeof(double));
    double *F_u2 = (double*)ALLOC(Md*sizeof(double));
    double *mass = (double*)ALLOC(Md*sizeof(double));
    for(i = 1; i <= Ncell; i++)//center cell is cell 0
	mass[i] = DD[i] * vol[i];

    for(k = 1; k <= N; k++)
	{
	    tic = clock();

	    if (time_c >= time_plot[nt] && nt < (*N_plot-1))
		{
#ifndef NOTECPLOT
		    FV.RHO   = CV.RHO[nt];
		    FV.U     = CV.U[nt];
		    FV.P     = CV.P[nt];
		    FV.gamma = CV.gamma[0];
		    file_radial_write_TEC(FV, *smv, problem, time_plot[nt]);
#endif
		    for(i = 0; i < Md; ++i)
			{
			    CV.RHO[nt+1][i] = CV.RHO[nt][i];
			    CV.U[nt+1][i]   =   CV.U[nt][i];
			    CV.E[nt+1][i]   =   CV.E[nt][i];  
			    CV.P[nt+1][i]   =   CV.P[nt][i];
			}
		    DD = CV.RHO[nt+1];
		    UU = CV.U[nt+1];
		    PP = CV.P[nt+1];
		    EE = CV.E[nt+1];
		    nt++;
		}

	    Smax_dr = 0.0; // S_max/dr = 0.0

	    VIP_limiter_radial   (Ncell, (_Bool)(k-1), DmU, TmV, UU, smv);
	    minmod_limiter_radial(Ncell, (_Bool)(k-1), DmD,      DD, smv);
	    minmod_limiter_radial(Ncell, (_Bool)(k-1), DmP,      PP, smv);

	    UU[Ncell+1]  = UU[Ncell];
	    DD[Ncell+1]  = DD[Ncell];
	    PP[Ncell+1]  = PP[Ncell];
	    EE[Ncell+1]  = EE[Ncell];
	    // TmV[0]       = 0.0;
	    DmU[0]       = 0.0;
	    DmD[0]       = 0.0;
	    DmP[0]       = 0.0;
	    // TmV[Ncell+1] = TmV[Ncell];
	    DmU[Ncell+1] = 0.0;
	    DmD[Ncell+1] = 0.0;
	    DmP[Ncell+1] = 0.0;

	    for(i = 0; i <= Ncell; i++)
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
		    if(ifvar_check(&ifv_L, &ifv_R, 1))
			{
			    printf(" on [%d, %d] (t_n, x).\n", k, i);
			    goto return_NULL;
			}

		    GRPsolverRLag(wave_speed, dire, mid, &ifv_L, &ifv_R, Rb[i+1], M, eps, eps);

		    if(star_dire_check(mid, dire, 1))
			{
			    printf(" on [%d, %d] (t_n, x).\n", k, i);
			    stop_t = true;
			}

		    Smax_dr = Smax_dr > fabs(wave_speed[0])/Ddr[i]   ? Smax_dr : fabs(wave_speed[0])/Ddr[i];
		    Smax_dr = Smax_dr > fabs(wave_speed[1])/Ddr[i+1] ? Smax_dr : fabs(wave_speed[1])/Ddr[i+1];

		    Umin[i+1]  = mid[1];
		    Pmin[i+1]  = mid[2];
		    DLmin[i+1] = mid[0];
		    DRmin[i+1] = mid[3];
		    U_t[i+1]   = dire[1];
		    P_t[i+1]   = dire[2];
		    DL_t[i+1]  = dire[0];
		    DR_t[i+1]  = dire[3];
		}

	    if(isfinite(Timeout) || !isfinite(config[16]) || config[16] <= 0.0) //compute for time step
		{
		    dt = CFL/Smax_dr;
		    if(dt < eps)
			{
			    printf("\nThe length of the time step is so small on [%d, %g, %g] (t_n, time_c, dt)\n", k, time_c, dt);
			    stop_t = true;
			}
		    else if((time_c + dt) > (Timeout - eps))
			dt = Timeout - time_c;
		    else if(!isfinite(dt))
			{
			    printf("NAN or INFinite error on [%d, %g, %g] (t_n, time_c, dt) - CFL\n", k, time_c, dt); 
			    goto return_NULL;
			}
		}

	    for(i = 1; i <= Ncell; i++)
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

		    AcousticRLagTangent(dire, mid, &ifv_L, &ifv_R, 0.5*(Rb[i]+Rb[i+1])/cos(0.5*dtheta), M, eps);
		    F_u2[i]  = mid[2]+dt*dire[2];
		    VLmin[i] = (ifv_L.U*cos(0.5*dtheta)+ifv_L.V*sin(0.5*dtheta)+dt*dire[1])*sin(0.5*dtheta);
		}

	    for(i = 0; i <= Ncell; i++)
		{
		    Umin[i+1] += 0.5 * dt * U_t[i+1];
		    Pmin[i+1] += 0.5 * dt * P_t[i+1];

		    F_u[i+1] = Pmin[i+1];
		    F_e[i+1] = Pmin[i+1] * Umin[i+1];

		    Rb_NStep = Rb[i+1] + Umin[i+1]*dt;
		    Lb_NStep = 2.0 * Rb_NStep * tan(0.5*dtheta);
		    Lbh[i+1] = 0.5 * (Lb[i+1] + Lb_NStep);
		    Rbh[i+1] = (Rb[i+1]*(2.*Lb[i+1]+Lb_NStep)+Rb_NStep*(Lb[i+1]+2.*Lb_NStep))/(3.*(Lb[i+1]+Lb_NStep));
		    Rb[i+1]  = Rb_NStep;
		    Lb[i+1]  = Lb_NStep;
		    RR[i]    = Rb[i+1]-(2.*Lb[i]+Lb[i+1])/(3.*(Lb[i]+Lb[i+1]))*(Rb[i+1]-Rb[i]);
		    /*
		      Rb_side[i] =0.25*(Lb[i]+Lb[i+1])  /sin(0.5*dtheta);
		      Lb_side[i] =0.5 *(Lb[i+1]-Lb[i])  /sin(0.5*dtheta);
		      Rbh_side[i]=0.25*(Lbh[i]+Lbh[i+1])/sin(0.5*dtheta);
		      Lbh_side[i]=0.5 *(Lbh[i+1]-Lbh[i])/sin(0.5*dtheta);
		      Sh[i]=Rbh[i+1]*Lbh[i+1]-Rbh[i]*Lbh[i]-sin(dtheta)*Rbh_side[i]*Lbh_side[i];
		    */

		    Umin[i+1]  += dt *  U_t[i+1] * 0.5;
		    Pmin[i+1]  += dt *  P_t[i+1] * 0.5;
		    DLmin[i+1] += dt * DL_t[i+1];
		    DRmin[i+1] += dt * DR_t[i+1];
		}

	    radial_mesh_update(smv);

	    DD[0] = mass[0]/vol[0];
	    UU[0] = 0.0;
	    EE[0] = EE[0] - dt/mass[0]*(F_e[1]*Rbh[1]*Lbh[1]);
	    for(i = 1; i <= Ncell; i++) //m=2
		{
		    DD[i] = mass[i]/vol[i];
		    UU[i] = UU[i] - dt/mass[i]*((F_u[i+1]-F_u2[i])*Rbh[i+1]*Lbh[i+1]-(F_u[i]-F_u2[i])*Rbh[i]*Lbh[i]);
		    EE[i] = EE[i] - dt/mass[i]*( F_e[i+1]         *Rbh[i+1]*Lbh[i+1]- F_e[i]         *Rbh[i]*Lbh[i]);
		}
	    for(i = 0; i <= Ncell; i++) //m=2
		{
		    PP[i] = (EE[i] - 0.5*UU[i]*UU[i]) * (GammaGamma[i]-1.0) * DD[i];
		    if (PP[i] < eps)
			{
			    printf("p<0.0 error on [%d, %d] (t_n, x) - Update\n", k, i);
			    stop_t = true;
			}
		}

	    for(i = 1; i <= Ncell; i++)
		{
		    DmU[i]=(Umin[i+1] -Umin[i]) /Ddr[i];
		    DmP[i]=(Pmin[i+1] -Pmin[i]) /Ddr[i];
		    DmD[i]=(DLmin[i+1]-DRmin[i])/Ddr[i];
		}
	    DmU[0]      =(Umin[1]-UU[0]) /dRc[0];
	    DmP[0]      =(Pmin[1]-PP[0]) /dRc[0];
	    DmD[0]      =(DLmin[1]-DD[0])/dRc[0];
	    DmU[Ncell+1]=(Umin[Ncell+1] -UU[Ncell])/dRc[Ncell];
	    DmP[Ncell+1]=(Pmin[Ncell+1] -PP[Ncell])/dRc[Ncell];
	    DmD[Ncell+1]=(DLmin[Ncell+1]-DD[Ncell])/dRc[Ncell];

	    time_c=time_c+dt;
	    if(isfinite(Timeout))
		DispPro(time_c*100.0/Timeout, k);
	    else
		DispPro(k*100.0/N, k);
	    if(stop_t || time_c > (Timeout - eps) || !isfinite(time_c))
		break;

	    toc = clock();
	    cpu_time_sum += ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;
	    cpu_time[nt]  = cpu_time_sum;
	}

    printf("\nTime is up at time step %d.\n", k);
    printf("The cost of CPU time for 1D-GRP Eulerian scheme for this problem is %g seconds.\n", cpu_time_sum);

 return_NULL:
    config[5] = (double)k;
    *N_plot = nt+1;
    if(isfinite(time_c))
	time_plot[nt] = time_c;
    else if(isfinite(Timeout))
	time_plot[nt] = Timeout;
    else if(isfinite(dt))
	time_plot[nt] = k*dt;

    DD = NULL;
    UU = NULL;
    PP = NULL;
    EE = NULL;
    GammaGamma = NULL;
    FREE(DmD);
    FREE(DmU);
    FREE(TmV);
    FREE(DmP);
    FREE(Umin);
    FREE(VLmin);
    FREE(Pmin);
    FREE(DLmin);
    FREE(DRmin);
    FREE(U_t);
    FREE(P_t);
    FREE(DL_t);
    FREE(DR_t);
    free(F_u);
    free(F_e);
    FREE(mass);
}
