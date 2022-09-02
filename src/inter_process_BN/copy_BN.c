#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/finite_volume.h"

//center_var to U
void BN_C2U(struct center_var C, double *U, int i, int j, int x_or_y)
{	
    U[0] = C.ZRHO_gC[i][j];
    U[4] = C.ZRHO_sC[i][j];
    switch(x_or_y) {
    case 0: //x direction
	U[1] = C.RHO_U_gC[i][j];
	U[2] = C.RHO_V_gC[i][j];
	U[5] = C.RHO_U_sC[i][j];
	U[6] = C.RHO_V_sC[i][j];	
	break;
    case 1: //y direction		
	U[1] = C.RHO_V_gC[i][j];
	U[2] = C.RHO_U_gC[i][j];
	U[5] = C.RHO_V_sC[i][j];
	U[6] = C.RHO_U_sC[i][j];	
	break;
    }
    U[3] = C.E_gC[i][j]-0.5*U[2]*U[2]/C.ZRHO_gC[i][j];
    U[7] = C.E_sC[i][j]-0.5*U[6]*U[6]/C.ZRHO_sC[i][j];
}

//U_L + U_R to primitive_var in center_var
void BN_ULR2prim(struct U_var U_L, struct U_var U_R, struct center_var C, int i, int j, int x_or_y)
{	
    C.RHO_gC[i][j] = 0.5*(U_L.rho_g+U_R.rho_g);
    C.RHO_sC[i][j] = 0.5*(U_L.rho_s+U_R.rho_s);
    C.P_gC[i][j]   = 0.5*(U_L.p_g +U_R.p_g);
    C.P_sC[i][j]   = 0.5*(U_L.p_s +U_R.p_s);
    switch(x_or_y) {		   		
    case 0: //x direction
	C.U_gC[i][j]   = 0.5*(U_L.u_g +U_R.u_g);
	C.V_gC[i][j]   = 0.5*(U_L.v_g +U_R.v_g);
	C.U_sC[i][j]   = 0.5*(U_L.u_s +U_R.u_s);
	C.V_sC[i][j]   = 0.5*(U_L.v_s +U_R.v_s);
	break;
    case 1: //y direction	
	C.U_gC[i][j]   = 0.5*(U_L.v_g +U_R.v_g);
	C.V_gC[i][j]   = 0.5*(U_L.u_g +U_R.u_g);
	C.U_sC[i][j]   = 0.5*(U_L.v_s +U_R.v_s);
	C.V_sC[i][j]   = 0.5*(U_L.u_s +U_R.u_s);
	break;
    }
}

//U_L + U_R to conservative_var in center_var
void BN_ULR2cons(struct U_var U_L, struct U_var U_R, struct center_var C, int i, int j, int x_or_y)
{	
    C.ZRHO_gC[i][j]  = 0.5*(U_L.U_rho_g+U_R.U_rho_g);
    C.ZRHO_sC[i][j]  = 0.5*(U_L.U_rho_s+U_R.U_rho_s);
    C.E_gC[i][j]     = 0.5*(U_L.U_e_g +U_R.U_e_g);
    C.E_sC[i][j]     = 0.5*(U_L.U_e_s +U_R.U_e_s);
    switch(x_or_y) {		   		
    case 0: //x direction
	C.RHO_U_gC[i][j] = 0.5*(U_L.U_u_g +U_R.U_u_g);
	C.RHO_V_gC[i][j] = 0.5*(U_L.U_v_g +U_R.U_v_g);
	C.RHO_U_sC[i][j] = 0.5*(U_L.U_u_s +U_R.U_u_s);
	C.RHO_V_sC[i][j] = 0.5*(U_L.U_v_s +U_R.U_v_s);
	break;
    case 1: //y direction	
	C.RHO_U_gC[i][j] = 0.5*(U_L.U_v_g +U_R.U_v_g);
	C.RHO_V_gC[i][j] = 0.5*(U_L.U_u_g +U_R.U_u_g);
	C.RHO_U_sC[i][j] = 0.5*(U_L.U_v_s +U_R.U_v_s);
	C.RHO_V_sC[i][j] = 0.5*(U_L.U_u_s +U_R.U_u_s);
	break;
    }
}

void RI_LR_ave(struct RI_var *RI, struct RI_var RI_L,struct RI_var RI_R)
{
    RI->eta_g=0.5*(RI_L.eta_g+RI_R.eta_g);
    RI->Q    =0.5*(RI_L.Q    +RI_R.Q);
    RI->P    =0.5*(RI_L.P    +RI_R.P);
    RI->H    =0.5*(RI_L.H    +RI_R.H);
    RI->rho_s=0.5*(RI_L.rho_s+RI_R.rho_s);
    RI->u_s  =0.5*(RI_L.u_s  +RI_R.u_s); 
}

void BN_RI2Cx(struct RI_var RI, struct center_var C, int i, int j)
{
    C.Q_xd[i][j]=RI.Q;
    C.P_xd[i][j]=RI.P;
    C.H_xd[i][j]=RI.H;
    C.eta_g_xd[i][j]=RI.eta_g;
}

void BN_RI2Cy(struct RI_var RI, struct center_var C, int i, int j)
{
    C.Q_yd[i][j]=RI.Q;
    C.P_yd[i][j]=RI.P;
    C.H_yd[i][j]=RI.H;
    C.eta_g_yd[i][j]=RI.eta_g;
}

void GRP_var_init(struct GRP_LR_var *G, struct slope_var SV, struct U_var U, double d, int i, int j, int pm_xy)
{
    G->rho_gx=SV.RHO_gx[i][j];
    G->p_gx  =SV.P_gx[i][j];
    G->rho_sx=SV.RHO_sx[i][j];
    G->p_sx  =SV.P_sx[i][j];	
    G->rho_gy=SV.RHO_gy[i][j];
    G->p_gy  =SV.P_gy[i][j];
    G->rho_sy=SV.RHO_sy[i][j];
    G->p_sy  =SV.P_sy[i][j];
    if (pm_xy < 2) { //x-direction			
	G->u_gx=SV.U_gx[i][j];
	G->v_gx=SV.V_gx[i][j];
	G->u_sx=SV.U_sx[i][j];
	G->v_sx=SV.V_sx[i][j];
	G->u_gy=SV.U_gy[i][j];	
	G->v_gy=SV.V_gy[i][j];
	G->u_sy=SV.U_sy[i][j];
	G->v_sy=SV.V_sy[i][j];
    }
    else { //y-direction
	G->u_gx=SV.V_gx[i][j];
	G->v_gx=SV.U_gx[i][j];
	G->u_sx=SV.V_sx[i][j];
	G->v_sx=SV.U_sx[i][j];
	G->u_gy=SV.V_gy[i][j];	
	G->v_gy=SV.U_gy[i][j];
	G->u_sy=SV.V_sy[i][j];
	G->v_sy=SV.U_sy[i][j];		
    }
    switch(pm_xy) {
    case 0: //x-direction: left side
	G->rho_g =U.rho_g+d/2*G->rho_gx;
	G->p_g =U.p_g+d/2*G->p_gx;
	G->u_g =U.u_g+d/2*G->u_gx;
	G->v_g =U.v_g+d/2*G->v_gx;
	G->rho_s =U.rho_s+d/2*G->rho_sx;
	G->p_s =U.p_s+d/2*G->p_sx;
	G->u_s =U.u_s+d/2*G->u_sx;
	G->v_s =U.v_s+d/2*G->v_sx;
	break;
    case 1: //x-direction: right side
	G->rho_g =U.rho_g-d/2*G->rho_gx;
	G->p_g =U.p_g-d/2*G->p_gx;
	G->u_g =U.u_g-d/2*G->u_gx;
	G->v_g =U.v_g-d/2*G->v_gx;
	G->rho_s =U.rho_s-d/2*G->rho_sx;
	G->p_s =U.p_s-d/2*G->p_sx;
	G->u_s =U.u_s-d/2*G->u_sx;
	G->v_s =U.v_s-d/2*G->v_sx;
	break;
    case 2: //y-direction: left side
	G->rho_g =U.rho_g+d/2*G->rho_gy;
	G->p_g =U.p_g+d/2*G->p_gy;
	G->u_g =U.u_g+d/2*G->u_gy;
	G->v_g =U.v_g+d/2*G->v_gy;
	G->rho_s =U.rho_s+d/2*G->rho_sy;
	G->p_s =U.p_s+d/2*G->p_sy;
	G->u_s =U.u_s+d/2*G->u_sy;
	G->v_s =U.v_s+d/2*G->v_sy;
	break;
    case 3: //y-direction: right side
	G->rho_g =U.rho_g-d/2*G->rho_gy;
	G->p_g =U.p_g-d/2*G->p_gy;
	G->u_g =U.u_g-d/2*G->u_gy;
	G->v_g =U.v_g-d/2*G->v_gy;
	G->rho_s =U.rho_s-d/2*G->rho_sy;
	G->p_s =U.p_s-d/2*G->p_sy;
	G->u_s =U.u_s-d/2*G->u_sy;
	G->v_s =U.v_s-d/2*G->v_sy;
	break;
    }
}
	
void GRP_RI_var_init(struct GRP_LR_var *G, struct slope_var SV, struct center_var C, double d, int i, int j, int pm_xy)
{
    G->Qx=SV.Q_x[i][j];
    G->Px=SV.P_x[i][j];
    G->Hx=SV.H_x[i][j];
    G->eta_gx=SV.eta_g_x[i][j];
    G->Qy=SV.Q_y[i][j];
    G->Py=SV.P_y[i][j];
    G->Hy=SV.H_y[i][j];
    G->eta_gy=SV.eta_g_y[i][j];
    switch(pm_xy) {
    case 0: //x-direction: left side
	G->Q     =C.Q_xd[i][j]+d/2*G->Qx;
	G->P     =C.P_xd[i][j]+d/2*G->Px;
	G->H     =C.H_xd[i][j]+d/2*G->Hx;
	G->eta_g =C.eta_g_xd[i][j]+d/2*G->eta_gx;
	break;
    case 1: //x-direction: right side
	G->Q     =C.Q_xd[i][j]-d/2*G->Qx;
	G->P     =C.P_xd[i][j]-d/2*G->Px;
	G->H     =C.H_xd[i][j]-d/2*G->Hx;
	G->eta_g =C.eta_g_xd[i][j]-d/2*G->eta_gx;			
	break;			
    case 2: //y-direction: left side
	G->Q     =C.Q_yd[i][j]+d/2*G->Qy;
	G->P     =C.P_yd[i][j]+d/2*G->Py;
	G->H     =C.H_yd[i][j]+d/2*G->Hy;
	G->eta_g =C.eta_g_yd[i][j]+d/2*G->eta_gy;
	break;
    case 3: //y-direction: right side
	G->Q     =C.Q_yd[i][j]-d/2*G->Qy;
	G->P     =C.P_yd[i][j]-d/2*G->Py;
	G->H     =C.H_yd[i][j]-d/2*G->Hy;
	G->eta_g =C.eta_g_yd[i][j]-d/2*G->eta_gy;			
	break;
    }	
}

void G_LR_RI2U(struct GRP_LR_var *G, double z_s, int x_or_y)
{	
    struct U_var U;
    struct RI_var RI;
    RI.Q=G->Q;
    RI.P=G->P;
    RI.H=G->H;
    RI.eta_g=G->eta_g;
    RI.z_s=z_s;
    RI.u_s=G->u_s;
    RI.rho_s=G->rho_s;
    RI2U_cal(&U, &RI, z_s, G->rho_g);
    G->rho_g =U.rho_g;
    G->rho_s =U.rho_s;
    G->p_g =U.p_g;
    G->p_s =U.p_s;		
    G->u_g =U.u_g;
    G->u_s =U.u_g;
    /*	
	G->rho_gx =0.0;
	G->p_gx =0.0;  
	G->u_gx =0.0;
	G->v_gx =0.0;
	G->rho_sx =0.0;
	G->p_sx =0.0;
	G->u_sx =0.0;
	G->v_sx =0.0;
	G->rho_gy =0.0;
	G->p_gy =0.0;
	G->u_gy =0.0;
	G->v_gy =0.0;
	G->rho_sy =0.0;
	G->p_sy =0.0;
	G->u_sy =0.0;
	G->v_sy =0.0;
    */
}
/* x方向的边界条件
 */
void boundary_cond_x(struct center_var C, int cond, int l)
{
    const int n_y = (int)config[14]+2, n_x = (int)config[13]+2;
    int i,k;
    double ***p;
    for(i = 0; i < n_y; ++i) {
	C.Z_sC[i][n_x-2]  = C.Z_sC[i][n_x-3];
	for(k=0, p=&C.Z_sC; k<sizeof(struct center_var)/sizeof(double **); k++, p++) {
	    if (cond != 1 && k < 1) {
		if (cond != -1 || l < 2)
		    (*p)[i][n_x-1] = (*p)[i][n_x-2];
		(*p)[i][0]     = (*p)[i][1];		    
	    }
	    else if (cond != 1) {
		if (cond != -1 || l < 2)
		    (*p)[i][n_x-1] = (*p)[i][n_x-2];
		(*p)[i][1]     = (*p)[i][2];		    
		(*p)[i][0]     = (*p)[i][1];	    
	    }
	    else if (k < 1) { // wall condition
		(*p)[i][n_x-1] = (*p)[i][1];		
		(*p)[i][0]     = (*p)[i][n_x-2];
	    }
	    else {
		(*p)[i][n_x-1] = (*p)[i][2];
		(*p)[i][0]     = (*p)[i][n_x-3];
		(*p)[i][1]     = (*p)[i][n_x-2];	
	    }		
	}
    }
}
/* y方向的边界条件
 */
void boundary_cond_y(struct center_var C, int cond, int l)
{
    const int n_y = (int)config[14]+2, n_x = (int)config[13]+2;
    int j,k;
    double ***p;
    for(j = 0; j < n_x; ++j) {
	for(k=0, p=&C.Z_sC; k<sizeof(struct center_var)/sizeof(double **); k++, p++) {
	    if (cond != 1 && k < 1) {			    
		if (cond != -1 || l < 2)
		    (*p)[n_y-1][j] = (*p)[n_y-2][j];
		(*p)[0][j]     = (*p)[1][j];
	    }
	    else if (cond != 1) {
		if (cond != -1 || l < 2)
		    (*p)[n_y-1][j] = (*p)[n_y-2][j];
		(*p)[1][j]     = (*p)[2][j];		    
		(*p)[0][j]     = (*p)[1][j];		    
	    }
	    else if (k < 1) { // wall condition
		(*p)[n_y-1][j] = (*p)[1][j];
		(*p)[0][j]     = (*p)[n_y-2][j];		    
	    }
	    else {
		(*p)[n_y-1][j] = (*p)[2][j];		
		(*p)[0][j]     = (*p)[n_y-3][j];
		(*p)[1][j]     = (*p)[n_y-2][j];   
	    }
	}
    }
}

void boundary_cond_slope_x(struct slope_var SV, int cond, int l)
{
    const int n_y = (int)config[14]+2, n_x = (int)config[13]+2;
    int i,k;
    double ***p;
    for(i = 0; i < n_y; ++i) {
	for(k=0, p=&SV.Z_sx; k<sizeof(struct slope_var)/sizeof(double **); k++, p++) {
	    if (cond != 1 && k < 1) {
		if (cond != -1 || l < 2)
		    (*p)[i][n_x-1] = (*p)[i][n_x-2];
		(*p)[i][0]     = (*p)[i][1];		    
	    }
	    else if (cond != 1) {
		if (cond != -1 || l < 2)
		    (*p)[i][n_x-1] = (*p)[i][n_x-2];
		(*p)[i][1]     = (*p)[i][2];		    
		(*p)[i][0]     = (*p)[i][1];		    
	    }
	    else if (k < 1) { // wall condition
		(*p)[i][n_x-1] = (*p)[i][1];		
		(*p)[i][0]     = (*p)[i][n_x-2];
	    }
	    else {
		(*p)[i][n_x-1] = (*p)[i][2];
		(*p)[i][0]     = (*p)[i][n_x-3];
		(*p)[i][1]     = (*p)[i][n_x-2];				
	    }
	}
    }
}

void boundary_cond_slope_y(struct slope_var SV, int cond, int l)
{
    const int n_y = (int)config[14]+2, n_x = (int)config[13]+2;
    int j,k;
    double ***p;
    for(j = 0; j < n_x; ++j) {		
	for(k=0, p=&SV.Z_sx; k<sizeof(struct slope_var)/sizeof(double **); k++, p++) {
	    if (cond != 1 && k < 2) {			    
		if (cond != -1 || l < 2)
		    (*p)[n_y-1][j] = (*p)[n_y-2][j];
		(*p)[0][j]     = (*p)[1][j];
	    }
	    else if (cond != 1) {
		if (cond != -1 || l < 2)
		    (*p)[n_y-1][j] = (*p)[n_y-2][j];
		(*p)[1][j]     = (*p)[2][j];		    
		(*p)[0][j]     = (*p)[1][j];		    
	    }
	    else if (k < 2) { // wall condition
		(*p)[n_y-1][j] = (*p)[1][j];
		(*p)[0][j]     = (*p)[n_y-2][j];		    
	    }
	    else {
		(*p)[n_y-1][j] = (*p)[2][j];		
		(*p)[0][j]     = (*p)[n_y-3][j];
		(*p)[1][j]     = (*p)[n_y-2][j];		    
	    }
	}
    }
}
