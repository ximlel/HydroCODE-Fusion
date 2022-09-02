#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"
#include "../include/tools.h"


//Newton-Rapshon iteration
void NewtonRapshon(double * x_star, double * err, double fun, double dfun, double eps)
{
    double d;
    if (fabs(fun) <= eps)
	d = 0.0;
    else
	d = -fun/dfun;
    * x_star = * x_star + d;
    * err = fabs(d);
}

static inline double V_norm(double * x)
{
    return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
}

void NewtonRapshon_matrix(double * x_star, double * err, double * fun, double * dfun, double eps)
{
    double d[4]={0.0};
    rinv(dfun,4); //Matrix inv of dfun
    int i,j;
    if (V_norm(fun) > eps) {
	for(i=0; i<4; i++)
	    for(j=0; j<4; j++)
		d[i]-=fun[j]*(*(dfun+i+j*4));
    }
    //V_add(x_star, x_star, d);
    for(i=0; i<4; i++)
	x_star[i] = x_star[i]+d[i];
    * err = V_norm(d);
}

static void NewtonRapshon_matrix2(double * x_star, double * err, double * fun, double * dfun, double eps)
{
    double d[2]={0.0};
    rinv(dfun,2); //Matrix inv of dfun
    int i,j;
    if (V_norm(fun) > eps) {
	for(i=0; i<2; i++)
	    for(j=0; j<2; j++)
		d[i]-=fun[j]*(*(dfun+i+j*2));
    }
    //V_add(x_star, x_star, d);
printf("\ndfun:%lf, %lf, %lf, %lf\n",d[0],d[1],x_star[0],x_star[1]);
    for(i=0; i<2; i++)
	x_star[i] = x_star[i]+d[i];
    * err = V_norm(d);
}

//From Riemann invariants to calculate var U.
void RI2U_cal(struct U_var * U, const struct RI_var * RI, double z_s, const double rho_g_start)
{
    const double eps = config[4];
    const double gama_g = config[106];
    double z_g = 1.0-z_s;
    double rho_s=RI->rho_s;
    double u_s=RI->u_s;
    double Q=RI->Q;
    double P=RI->P;
    double H=RI->H;
    double eta_g=RI->eta_g;
    int it_max = 500, k = 0;
    double err1 = 1e50;
    double rho_g=rho_g_start;
    double fun, dfun;
    while (k<it_max && err1>eps) {					
	fun  = H-0.5*pow(Q/z_g,2)/pow(rho_g,2)-gama_g/(gama_g-1.0)*eta_g*pow(rho_g,gama_g-1.0);
	dfun = pow(Q/z_g,2)/pow(rho_g,3)-gama_g*eta_g*pow(rho_g,gama_g-2.0);
	NewtonRapshon(&rho_g, &err1, fun,dfun,eps);
	k=k+1;
    }
    if(k>=it_max)
        printf("\nRI2P_err:%lf! %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", err1, z_s,rho_s,u_s,Q,P,H,eta_g,rho_g);
    U->p_g = pow(rho_g,gama_g)*eta_g;
    U->u_g = Q/z_g/rho_g+u_s;
    U->p_s = (P-Q*(U->u_g-u_s)-z_g*U->p_g)/z_s;
    U->rho_g = rho_g;
    U->rho_s = rho_s;
    U->u_s = u_s;
}

//From var U to calculate Riemann invariants.
void U2RI_cal(const struct U_var * U, struct RI_var * RI)
{
    const double eps = config[4];
    const double gama_g = config[106];
    double z_s = U->z_s, rho_s = U->rho_s, u_s = U->u_s, p_s = U->p_s, rho_g = U->rho_g, u_g = U->u_g, p_g = U->p_g;
    double z_g = 1.0-z_s;
	
    RI->eta_g=p_g/pow(rho_g,gama_g);
    RI->Q=z_g*rho_g*(u_g-u_s);
    RI->P=z_g*rho_g*pow(u_g-u_s,2)+z_g*p_g+z_s*p_s;
    RI->H=0.5*pow(u_g-u_s,2)+gama_g/(gama_g-1.0)*p_g/rho_g;
    RI->rho_s=rho_s;
    RI->u_s=u_s;
}

//compute primitive var
static void primitive_comp_old(double * U, struct U_var * U_L, struct U_var * U_R, double z_sL, double z_sR, double z_sL_out, double z_sR_out, double area_L, double area_R)
{   
    double z_gL=1-z_sL;
    double z_gR=1-z_sR;
    const double gama_g = config[106], gama_s = config[6];
    double eps = config[4];
    double z_s = area_L*z_sL+area_R*z_sR;
    double z_g = 1.0-z_s;
    double U1=U[0], U2=U[1], U3=U[3], U4=U[4], U5=U[5], U6=U[7];
    double rho_gR = U1/z_g;
    double u_gR  = U2/U1;
    double p_gR  = (U3/z_g - 0.5*rho_gR*pow(u_gR,2))*(gama_g-1.0);
    double rho_s  = U4/z_s;
    double u_s   = U5/U4;
    double p_sR  = (U6/z_s - 0.5*rho_s*pow(u_s,2))*(gama_s-1.0);
    U_L->v_g = U[2]/U1;
    U_R->v_g = U_L->v_g;	
    U_L->v_s = U[6]/U4;
    U_R->v_s = U_L->v_s;
    U_L->z_s = z_sL;
    U_R->z_s = z_sR;
    double fun[4], dfun[4][4], x_star[4];
    int it_max = 5000, k = 0;
    double err2 = 1e50;
    while (k<it_max && err2>eps && fabs(z_sL-z_sR)>eps) {			
	fun[0] = U3-area_L*z_gL*(0.5*((U1-area_R*z_gR*rho_gR)/area_L/z_gL)*pow((U2-area_R*z_gR*rho_gR*u_gR)/(U1-area_R*z_gR*rho_gR),2.0)+(p_gR/pow(rho_gR,gama_g)*pow((U1-area_R*z_gR*rho_gR)/area_L/z_gL,gama_g))/(gama_g-1.0))-area_R*z_gR*(0.5*rho_gR*pow(u_gR,2.0)+p_gR/(gama_g-1.0));
	fun[1] = z_gL*((U1-area_R*z_gR*rho_gR)/area_L/z_gL)*(((U2-area_R*z_gR*rho_gR*u_gR)/(U1-area_R*z_gR*rho_gR))-u_s)-z_gR*rho_gR*(u_gR-u_s);
	fun[2] = z_gL*((U1-area_R*z_gR*rho_gR)/area_L/z_gL)*pow(((U2-area_R*z_gR*rho_gR*u_gR)/(U1-area_R*z_gR*rho_gR))-u_s,2.0)+z_gL*(p_gR/pow(rho_gR,gama_g)*pow((U1-area_R*z_gR*rho_gR)/area_L/z_gL,gama_g))+z_sL*(((U6-0.5*z_s*rho_s*pow(u_s,2.0))*(gama_s-1.0)-area_R*z_sR*p_sR)/area_L/z_sL)-z_gR*rho_gR*pow(u_gR-u_s,2.0)-z_gR*p_gR-z_sR*p_sR;
	fun[3] = 0.5*pow(((U2-area_R*z_gR*rho_gR*u_gR)/(U1-area_R*z_gR*rho_gR))-u_s,2.0)+gama_g/(gama_g-1.0)*(p_gR/pow(rho_gR,gama_g)*pow((U1-area_R*z_gR*rho_gR)/area_L/z_gL,gama_g))/((U1-area_R*z_gR*rho_gR)/area_L/z_gL)-0.5*pow(u_gR-u_s,2.0)-gama_g/(gama_g-1.0)*p_gR/rho_gR;
	dfun[0][0] = area_L*z_gL*((gama_g*p_gR*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/(pow(rho_gR,gama_g + 1.0)*(gama_g - 1.0)) - (area_R*z_gR*pow(U2 - area_R*rho_gR*z_gR*u_gR,2.0))/(2.0*area_L*z_gL*pow(U1 - area_R*rho_gR*z_gR,2.0)) + (area_R*z_gR*u_gR*(U2 - area_R*rho_gR*z_gR*u_gR))/(area_L*z_gL*(U1 - area_R*rho_gR*z_gR)) + (area_R*gama_g*p_gR*z_gR*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g - 1.0))/(area_L*pow(rho_gR,gama_g)*z_gL*(gama_g - 1.0))) - (area_R*z_gR*pow(u_gR,2.0))/2.0;
	dfun[1][0] = - (area_R*z_gR)/(gama_g - 1.0) - (area_L*z_gL*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/(pow(rho_gR,gama_g)*(gama_g - 1));
	dfun[2][0] = (area_R*rho_gR*z_gR*(U2 - area_R*rho_gR*z_gR*u_gR))/(U1 - area_R*rho_gR*z_gR) - area_R*rho_gR*z_gR*u_gR;
	dfun[3][0] = 0.0;
	dfun[0][1] = (area_R*z_gR*(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR)))/area_L - ((U1 - area_R*rho_gR*z_gR)*((area_R*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR) - (area_R*z_gR*(U2 - area_R*rho_gR*z_gR*u_gR))/pow(U1 - area_R*rho_gR*z_gR,2.0)))/area_L - z_gR*(u_gR - u_s);
	dfun[1][1] = 0.0;
	dfun[2][1] = - rho_gR*z_gR - (area_R*rho_gR*z_gR)/area_L;
	dfun[3][1] = 0.0;
	dfun[0][2] = (2*(U1 - area_R*rho_gR*z_gR)*((area_R*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR) - (area_R*z_gR*(U2 - area_R*rho_gR*z_gR*u_gR))/pow(U1 - area_R*rho_gR*z_gR,2.0))*(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR)))/area_L - z_gR*pow(u_gR - u_s,2.0) - (area_R*z_gR*pow(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR),2))/area_L - (gama_g*p_gR*z_gL*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/pow(rho_gR,gama_g + 1.0) - (area_R*gama_g*p_gR*z_gR*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g - 1.0))/(area_L*pow(rho_gR,gama_g));
	dfun[1][2] = (z_gL*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/pow(rho_gR,gama_g) - z_gR;
	dfun[2][2] = (2.0*area_R*rho_gR*z_gR*(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR)))/area_L - rho_gR*z_gR*(2.0*u_gR - 2.0*u_s);
	dfun[3][2] = - z_sR - (area_R*z_sR)/area_L;
	dfun[0][3] = ((area_R*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR) - (area_R*z_gR*(U2 - area_R*rho_gR*z_gR*u_gR))/pow(U1 - area_R*rho_gR*z_gR,2.0))*(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR)) + (gama_g*p_gR)/(pow(rho_gR,2.0)*(gama_g - 1.0)) - (area_L*pow(gama_g,2.0)*p_gR*z_gL*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/(pow(rho_gR,gama_g + 1.0)*(U1 - area_R*rho_gR*z_gR)*(gama_g - 1.0)) - (area_R*pow(gama_g,2.0)*p_gR*z_gR*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g - 1.0))/(pow(rho_gR,gama_g)*(U1 - area_R*rho_gR*z_gR)*(gama_g - 1.0)) + (area_L*area_R*gama_g*p_gR*z_gL*z_gR*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/(pow(rho_gR,gama_g)*pow(U1 - area_R*rho_gR*z_gR,2.0)*(gama_g - 1.0));
	dfun[1][3] = (area_L*gama_g*z_gL*pow((U1 - area_R*rho_gR*z_gR)/(area_L*z_gL),gama_g))/(pow(rho_gR,gama_g)*(U1 - area_R*rho_gR*z_gR)*(gama_g - 1)) - gama_g/(rho_gR*(gama_g - 1.0));
	dfun[2][3] = u_s - u_gR + (area_R*rho_gR*z_gR*(u_s - (U2 - area_R*rho_gR*z_gR*u_gR)/(U1 - area_R*rho_gR*z_gR)))/(U1 - area_R*rho_gR*z_gR);
	dfun[3][3] = 0.0;
	x_star[0] = rho_gR;
	x_star[1] = p_gR;
	x_star[2] = u_gR;
	x_star[3] = p_sR;
	NewtonRapshon_matrix(x_star, &err2, fun, dfun[0], eps);
	rho_gR=fmax(x_star[0],eps);
	rho_gR=fmin(rho_gR,U1/area_R/z_gR-eps);
	p_gR =fmax(x_star[1],eps);
	u_gR =x_star[2];
	p_sR =fmax(x_star[3],eps);
	p_sR =fmin(p_sR,(U6-0.5*z_s*rho_s*pow(u_s,2))*(gama_s-1.0)/area_R/z_sR-eps);
	k=k+1;
    }
    if (k>=it_max)
        printf("\nRIeq_err:%lf! %lf, %lf, %lf, %lf, %lf, %lf\n",err2,z_sL,z_sR,rho_gR,p_gR,u_gR,p_sR);

    U_L->rho_g = (U1-area_R*z_gR*rho_gR)/area_L/z_gL;
    U_R->rho_g = rho_gR;
    U_L->u_g = (U2-area_R*z_gR*rho_gR*u_gR)/(U1-area_R*z_gR*rho_gR);
    U_R->u_g = u_gR;
    U_L->p_s = ((U6-0.5*z_s*rho_s*pow(u_s,2.0))*(gama_s-1)-area_R*z_sR*p_sR)/area_L/z_sL;
    U_L->p_g = p_gR/pow(rho_gR,gama_g)*pow((U1-area_R*z_gR*rho_gR)/area_L/z_gL,gama_g);
    U_R->p_s = p_sR;
    U_R->p_g = p_gR;	
    U_L->rho_s= rho_s;
    U_R->rho_s= rho_s;
    U_L->u_s = u_s;
    U_R->u_s = u_s;
    
    struct RI_var RI_L, RI_R;
    U2RI_cal(U_L, &RI_L);
    RI2U_cal(U_L, &RI_L, z_sL_out, U_L->rho_g);
    U2RI_cal(U_R, &RI_R);
    RI2U_cal(U_R, &RI_R, z_sR_out, U_R->rho_g);
    
    U_L->U_rho_g = (1.0-z_sL_out)*U_L->rho_g;
    U_L->U_u_g  = U_L->U_rho_g*U_L->u_g;
    U_L->U_v_g  = U_L->U_rho_g*U_L->v_g;
    U_L->U_e_g  = U_L->U_rho_g*(U_L->p_g/U_L->rho_g/(gama_g-1.0)+0.5*U_L->u_g*U_L->u_g+0.5*U_L->v_g*U_L->v_g);
    U_L->U_rho_s = z_sL_out*U_L->rho_s;
    U_L->U_u_s  = U_L->U_rho_s*U_L->u_s;
    U_L->U_v_s  = U_L->U_rho_s*U_L->v_s;
    U_L->U_e_s  = U_L->U_rho_s*(U_L->p_s/U_L->rho_s/(gama_s-1.0)+0.5*U_L->u_s*U_L->u_s+0.5*U_L->v_s*U_L->v_s);
    U_R->U_rho_g = (1.0-z_sR_out)*U_R->rho_g;
    U_R->U_u_g  = U_R->U_rho_g*U_R->u_g;
    U_R->U_v_g  = U_R->U_rho_g*U_R->v_g;
    U_R->U_e_g  = U_R->U_rho_g*(U_R->p_g/U_R->rho_g/(gama_g-1.0)+0.5*U_R->u_g*U_R->u_g+0.5*U_R->v_g*U_R->v_g);
    U_R->U_rho_s = z_sR_out*U_R->rho_s;
    U_R->U_u_s  = U_R->U_rho_s*U_R->u_s;
    U_R->U_v_s  = U_R->U_rho_s*U_R->v_s;
    U_R->U_e_s  = U_R->U_rho_s*(U_R->p_s/U_R->rho_s/(gama_s-1.0)+0.5*U_R->u_s*U_R->u_s+0.5*U_R->v_s*U_R->v_s);	
}

static void primitive_comp_bak(double * U, struct U_var * U_L, struct U_var * U_R, double z_sL, double z_sR, double z_sL_out, double z_sR_out, double area_L, double area_R)
{   
    double z_gL=1-z_sL;
    double z_gR=1-z_sR;
    const double gama_g = config[106], gama_s = config[6];
    double eps = config[4];
    double z_s = area_L*z_sL+area_R*z_sR;
    double z_g = 1.0-z_s;
    double U1=U[0], U2=U[1], U3=U[3], U4=U[4], U5=U[5], U6=U[7];
    double rho_gR = U1/z_g;
    double rho_s  = U4/z_s;
    double u_s   = U5/U4;
    U_L->v_g = U[2]/U1;
    U_R->v_g = U_L->v_g;	
    U_L->v_s = U[6]/U4;
    U_R->v_s = U_L->v_s;
    U_L->z_s = z_sL;
    U_R->z_s = z_sR;
    double fun, dfun, x_star;
    int it_max = 5000, k = 0;
    double err2 = 1e50;
    while (k<it_max && err2>eps && fabs(z_sL-z_sR)>eps) {			
	fun = ((U3 + (-0.1e1) * 0.5e0 * area_R * z_gR * rho_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + (-0.1e1) * 0.5e0 * (-area_R * z_gR * rho_gR + U1) * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1)) * (gama_g - 1) + 0.5e0 * area_R * ((((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.2e1 * U3 * area_R * gama_g) * pow(rho_gR, 0.3e1) * pow(z_gR, 0.3e1) + (((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + U2 * U2 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * rho_gR * z_gR + (pow(U1, 0.4e1) + 0.2e1 * pow(U1, 0.3e1) * U4 + U1 * U1 * U4 * U4) * u_s * u_s + ((-0.1e1) * 0.2e1 * U1 * U1 * U2 * U4 + (-0.1e1) * 0.2e1 * pow(U1, 0.3e1) * U2) * u_s + U1 * U1 * U2 * U2) * (gama_g - 0.1000000000e1) / z_gR / rho_gR / (-area_R * z_gR * rho_gR + U1) / gama_g / U1) / area_L / z_gL / pow((-area_R * z_gR * rho_gR + U1) / area_L / z_gL, gama_g) + 0.5e0 * ((((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.2e1 * U3 * area_R * gama_g) * pow(rho_gR, 0.3e1) * pow(z_gR, 0.3e1) + (((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + U2 * U2 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * rho_gR * z_gR + (pow(U1, 0.4e1) + 0.2e1 * pow(U1, 0.3e1) * U4 + U1 * U1 * U4 * U4) * u_s * u_s + ((-0.1e1) * 0.2e1 * U1 * U1 * U2 * U4 + (-0.1e1) * 0.2e1 * pow(U1, 0.3e1) * U2) * u_s + U1 * U1 * U2 * U2) * (gama_g - 0.1000000000e1) * pow(z_gR, -0.2e1) / rho_gR / (-area_R * z_gR * rho_gR + U1) / gama_g / U1 / pow(rho_gR, gama_g);

	dfun = (((-0.1e1) * 0.5e0 * area_R * z_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + 0.10e1 * area_R * ((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s) * (U2 - (U1 + U4) * u_s) / rho_gR + 0.5e0 * area_R * z_gR * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1) + (-0.1e1) * 0.10e1 * ((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s) * (U2 - (U1 + U4) * u_s) * area_L * area_R * z_gR / (-area_R * z_gR * rho_gR + U1)) * (gama_g - 1) + 0.5e0 * area_R * ((0.3e1 * ((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + 0.3e1 * ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.3e1 * 0.2e1 * U3 * area_R * gama_g) * rho_gR * rho_gR * pow(z_gR, 0.3e1) + (0.2e1 * ((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + 0.2e1 * ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + 0.2e1 * (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + 0.2e1 * (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + 0.2e1 * U2 * U2 * area_L * area_L * gama_g + 0.2e1 * (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * z_gR) * (gama_g - 0.1000000000e1) / z_gR / rho_gR / (-area_R * z_gR * rho_gR + U1) / gama_g / U1 + (-0.1e1) * 0.5e0 * area_R * ((((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.2e1 * U3 * area_R * gama_g) * pow(rho_gR, 0.3e1) * pow(z_gR, 0.3e1) + (((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + U2 * U2 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * rho_gR * z_gR + (pow(U1, 0.4e1) + 0.2e1 * pow(U1, 0.3e1) * U4 + U1 * U1 * U4 * U4) * u_s * u_s + ((-0.1e1) * 0.2e1 * U1 * U1 * U2 * U4 + (-0.1e1) * 0.2e1 * pow(U1, 0.3e1) * U2) * u_s + U1 * U1 * U2 * U2) * (gama_g - 0.1000000000e1) / z_gR * pow(rho_gR, -0.2e1) / (-area_R * z_gR * rho_gR + U1) / gama_g / U1 + 0.5e0 * area_R * area_R * ((((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.2e1 * U3 * area_R * gama_g) * pow(rho_gR, 0.3e1) * pow(z_gR, 0.3e1) + (((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + U2 * U2 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * rho_gR * z_gR + (pow(U1, 0.4e1) + 0.2e1 * pow(U1, 0.3e1) * U4 + U1 * U1 * U4 * U4) * u_s * u_s + ((-0.1e1) * 0.2e1 * U1 * U1 * U2 * U4 + (-0.1e1) * 0.2e1 * pow(U1, 0.3e1) * U2) * u_s + U1 * U1 * U2 * U2) * (gama_g - 0.1000000000e1) / rho_gR * pow(-area_R * z_gR * rho_gR + U1, -0.2e1) / gama_g / U1) / area_L / z_gL / pow((-area_R * z_gR * rho_gR + U1) / area_L / z_gL, gama_g) + ((U3 + (-0.1e1) * 0.5e0 * area_R * z_gR * rho_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + (-0.1e1) * 0.5e0 * (-area_R * z_gR * rho_gR + U1) * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1)) * (gama_g - 1) + 0.5e0 * area_R * ((((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.2e1 * U3 * area_R * gama_g) * pow(rho_gR, 0.3e1) * pow(z_gR, 0.3e1) + (((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + U2 * U2 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * rho_gR * z_gR + (pow(U1, 0.4e1) + 0.2e1 * pow(U1, 0.3e1) * U4 + U1 * U1 * U4 * U4) * u_s * u_s + ((-0.1e1) * 0.2e1 * U1 * U1 * U2 * U4 + (-0.1e1) * 0.2e1 * pow(U1, 0.3e1) * U2) * u_s + U1 * U1 * U2 * U2) * (gama_g - 0.1000000000e1) / z_gR / rho_gR / (-area_R * z_gR * rho_gR + U1) / gama_g / U1) * gama_g * area_R * z_gR / area_L / z_gL / pow((-area_R * z_gR * rho_gR + U1) / area_L / z_gL, gama_g) / (-area_R * z_gR * rho_gR + U1) + 0.5e0 * ((0.3e1 * ((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + 0.3e1 * ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.3e1 * 0.2e1 * U3 * area_R * gama_g) * rho_gR * rho_gR * pow(z_gR, 0.3e1) + (0.2e1 * ((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + 0.2e1 * ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + 0.2e1 * (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + 0.2e1 * (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + 0.2e1 * U2 * U2 * area_L * area_L * gama_g + 0.2e1 * (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * z_gR) * (gama_g - 0.1000000000e1) * pow(z_gR, -0.2e1) / rho_gR / (-area_R * z_gR * rho_gR + U1) / gama_g / U1 / pow(rho_gR, gama_g) + (-0.1e1) * 0.5e0 * ((((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.2e1 * U3 * area_R * gama_g) * pow(rho_gR, 0.3e1) * pow(z_gR, 0.3e1) + (((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + U2 * U2 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * rho_gR * z_gR + (pow(U1, 0.4e1) + 0.2e1 * pow(U1, 0.3e1) * U4 + U1 * U1 * U4 * U4) * u_s * u_s + ((-0.1e1) * 0.2e1 * U1 * U1 * U2 * U4 + (-0.1e1) * 0.2e1 * pow(U1, 0.3e1) * U2) * u_s + U1 * U1 * U2 * U2) * (gama_g - 0.1000000000e1) * pow(z_gR, -0.2e1) * pow(rho_gR, -0.2e1) / (-area_R * z_gR * rho_gR + U1) / gama_g / U1 / pow(rho_gR, gama_g) + 0.5e0 * ((((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.2e1 * U3 * area_R * gama_g) * pow(rho_gR, 0.3e1) * pow(z_gR, 0.3e1) + (((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + U2 * U2 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * rho_gR * z_gR + (pow(U1, 0.4e1) + 0.2e1 * pow(U1, 0.3e1) * U4 + U1 * U1 * U4 * U4) * u_s * u_s + ((-0.1e1) * 0.2e1 * U1 * U1 * U2 * U4 + (-0.1e1) * 0.2e1 * pow(U1, 0.3e1) * U2) * u_s + U1 * U1 * U2 * U2) * (gama_g - 0.1000000000e1) * area_R / z_gR / rho_gR * pow(-area_R * z_gR * rho_gR + U1, -0.2e1) / gama_g / U1 / pow(rho_gR, gama_g) + (-0.1e1) * 0.5e0 * ((((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.2e1 * U3 * area_R * gama_g) * pow(rho_gR, 0.3e1) * pow(z_gR, 0.3e1) + (((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + U2 * U2 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * rho_gR * z_gR + (pow(U1, 0.4e1) + 0.2e1 * pow(U1, 0.3e1) * U4 + U1 * U1 * U4 * U4) * u_s * u_s + ((-0.1e1) * 0.2e1 * U1 * U1 * U2 * U4 + (-0.1e1) * 0.2e1 * pow(U1, 0.3e1) * U2) * u_s + U1 * U1 * U2 * U2) * (gama_g - 0.1000000000e1) * pow(z_gR, -0.2e1) * pow(rho_gR, -0.2e1) / (-area_R * z_gR * rho_gR + U1) / U1 / pow(rho_gR, gama_g);

	NewtonRapshon(&rho_gR, &err2, fun, dfun, eps);
	rho_gR=fmax(rho_gR,eps);
	rho_gR=fmin(rho_gR,U1/area_R/z_gR-eps);
	k=k+1;
    }
    if (k>=it_max)
        printf("\nRIeq_err:%lf! %lf, %lf, %lf\n",err2,z_sL,z_sR,rho_gR);
        
    U_R->rho_g = rho_gR;
    U_L->rho_g = (U1-area_R*z_gR*rho_gR)/area_L/z_gL;
    U_R->u_g = (U2-(U1+U4)*u_s)/(z_gR*rho_gR)+u_s;
    U_L->u_g = (U2-(U1+U4)*u_s)/(U1-area_R*z_gR*rho_gR)*area_L+u_s;
    
    U_R->p_g = (-0.1e1) * 0.5e0 * ((((0.2e1 * area_R * area_R * gama_g + (0.2e1 * area_L - 0.1e1) * gama_g * area_R) * U1 + 0.2e1 * U4 * area_L * area_R * gama_g + 0.2e1 * U4 * area_R * area_R * gama_g) * u_s * u_s + ((-0.1e1) * 0.2e1 * U2 * area_L * area_R * gama_g + (-0.1e1) * 0.2e1 * U2 * area_R * area_R * gama_g) * u_s + 0.2e1 * U3 * area_R * gama_g) * pow(rho_gR, 0.3e1) * pow(z_gR, 0.3e1) + (((((-0.1e1) * 0.1e1 * gama_g + 0.1e1) * area_R * area_R + (-0.1e1) * 0.2e1 * area_R * gama_g + pow(area_L - 0.1e1, 0.2e1) * gama_g + (-0.1e1) * 0.1e1 * area_L * area_L) * U1 * U1 + ((0.2e1 * U4 + (-0.1e1) * 0.2e1 * U4 * gama_g) * area_R * area_R + (-0.1e1) * 0.2e1 * U4 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L + 0.2e1 * area_L * area_L) * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * area_L * area_L) * U1 + (-U4 * U4 * gama_g + U4 * U4) * area_R * area_R + U4 * U4 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U4 * U4 * area_L * area_L) * u_s * u_s + ((((-0.1e1) * 0.2e1 * U2 + 0.2e1 * U2 * gama_g) * area_R * area_R + 0.2e1 * U2 * area_R * gama_g + ((-0.1e1) * 0.2e1 * area_L * area_L + 0.2e1 * area_L) * U2 * gama_g + 0.2e1 * U2 * area_L * area_L) * U1 + ((-0.1e1) * 0.2e1 * U2 * U4 + 0.2e1 * U2 * U4 * gama_g) * area_R * area_R + 0.2e1 * U2 * U4 * area_L * area_L + (-0.1e1) * 0.2e1 * U2 * U4 * area_L * area_L * gama_g) * u_s + (-0.1e1) * 0.2e1 * U1 * U3 * gama_g + (-U2 * U2 * gama_g + U2 * U2) * area_R * area_R + U2 * U2 * area_L * area_L * gama_g + (-0.1e1) * 0.1e1 * U2 * U2 * area_L * area_L) * rho_gR * rho_gR * z_gR * z_gR + (((gama_g - 0.2e1) * area_R * pow(U1, 0.3e1) + (0.2e1 * U4 * gama_g + (-0.1e1) * 0.4e1 * U4) * area_R * U1 * U1 + (U4 * U4 * gama_g + (-0.1e1) * 0.2e1 * U4 * U4) * area_R * U1) * u_s * u_s + (((-0.1e1) * 0.2e1 * U2 * gama_g + 0.4e1 * U2) * area_R * U1 * U1 + (0.4e1 * U2 * U4 + (-0.1e1) * 0.2e1 * U2 * U4 * gama_g) * area_R * U1) * u_s + (U2 * U2 * gama_g + (-0.1e1) * 0.2e1 * U2 * U2) * area_R * U1) * rho_gR * z_gR + (pow(U1, 0.4e1) + 0.2e1 * pow(U1, 0.3e1) * U4 + U1 * U1 * U4 * U4) * u_s * u_s + ((-0.1e1) * 0.2e1 * U1 * U1 * U2 * U4 + (-0.1e1) * 0.2e1 * pow(U1, 0.3e1) * U2) * u_s + U1 * U1 * U2 * U2) * (gama_g - 0.1000000000e1) * pow(z_gR, -0.2e1) / rho_gR / (-area_R * z_gR * rho_gR + U1) / gama_g / U1;
    U_L->p_g = ((U3 + (-0.1e1) * 0.5e0 * area_R * z_gR * rho_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + (-0.1e1) * 0.5e0 * (-area_R * z_gR * rho_gR + U1) * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1)) * (gama_g - 1) - area_R * z_gR * U_R->p_g) / area_L / z_gL;
    U_R->p_s = (((z_gL * U_L->rho_g + (-0.1e1) * 0.1e1 * z_gR * rho_gR) * u_s * u_s + ((-0.1e1) * 0.2e1 * U_L->rho_g * z_gL * U_L->u_g + 0.2e1 * rho_gR * z_gR * U_R->u_g) * u_s + (U_L->rho_g * U_L->u_g * U_L->u_g + U_L->p_g) * z_gL + ((-0.1e1) * 0.1e1 * rho_gR * U_R->u_g * U_R->u_g + (-0.1e1) * 0.1e1 * U_R->p_g) * z_gR) * area_L + ((-0.1e1) * 0.5e0 * gama_s * rho_s + 0.5e0 * rho_s) * u_s * u_s + U6 * gama_s + (-0.1e1) * 0.1e1 * U6) / z_sR / (area_L + area_R);
    U_L->p_s = ((U6 + (-0.1e1) * 0.5e0 * rho_s * u_s * u_s) * (gama_s - 1) - area_R * (((z_gL * U_L->rho_g + (-0.1e1) * 0.1e1 * z_gR * rho_gR) * u_s * u_s + ((-0.1e1) * 0.2e1 * U_L->rho_g * z_gL * U_L->u_g + 0.2e1 * rho_gR * z_gR * U_R->u_g) * u_s + (U_L->rho_g * U_L->u_g * U_L->u_g + U_L->p_g) * z_gL + ((-0.1e1) * 0.1e1 * rho_gR * U_R->u_g * U_R->u_g + (-0.1e1) * 0.1e1 * U_R->p_g) * z_gR) * area_L + ((-0.1e1) * 0.5e0 * gama_s * rho_s + 0.5e0 * rho_s) * u_s * u_s + U6 * gama_s + (-0.1e1) * 0.1e1 * U6) / (area_L + area_R)) / area_L / z_sL;

    U_L->rho_s= rho_s;
    U_R->rho_s= rho_s;
    U_L->u_s = u_s;
    U_R->u_s = u_s;
    
    struct RI_var RI_L, RI_R;
    U2RI_cal(U_L, &RI_L);
    RI2U_cal(U_L, &RI_L, z_sL_out, U_L->rho_g);
    U2RI_cal(U_R, &RI_R);
    RI2U_cal(U_R, &RI_R, z_sR_out, U_R->rho_g);
    
    U_L->U_rho_g = (1.0-z_sL_out)*U_L->rho_g;
    U_L->U_u_g  = U_L->U_rho_g*U_L->u_g;
    U_L->U_v_g  = U_L->U_rho_g*U_L->v_g;
    U_L->U_e_g  = U_L->U_rho_g*(U_L->p_g/U_L->rho_g/(gama_g-1.0)+0.5*U_L->u_g*U_L->u_g+0.5*U_L->v_g*U_L->v_g);
    U_L->U_rho_s = z_sL_out*U_L->rho_s;
    U_L->U_u_s  = U_L->U_rho_s*U_L->u_s;
    U_L->U_v_s  = U_L->U_rho_s*U_L->v_s;
    U_L->U_e_s  = U_L->U_rho_s*(U_L->p_s/U_L->rho_s/(gama_s-1.0)+0.5*U_L->u_s*U_L->u_s+0.5*U_L->v_s*U_L->v_s);
    U_R->U_rho_g = (1.0-z_sR_out)*U_R->rho_g;
    U_R->U_u_g  = U_R->U_rho_g*U_R->u_g;
    U_R->U_v_g  = U_R->U_rho_g*U_R->v_g;
    U_R->U_e_g  = U_R->U_rho_g*(U_R->p_g/U_R->rho_g/(gama_g-1.0)+0.5*U_R->u_g*U_R->u_g+0.5*U_R->v_g*U_R->v_g);
    U_R->U_rho_s = z_sR_out*U_R->rho_s;
    U_R->U_u_s  = U_R->U_rho_s*U_R->u_s;
    U_R->U_v_s  = U_R->U_rho_s*U_R->v_s;
    U_R->U_e_s  = U_R->U_rho_s*(U_R->p_s/U_R->rho_s/(gama_s-1.0)+0.5*U_R->u_s*U_R->u_s+0.5*U_R->v_s*U_R->v_s);	
}

void primitive_comp(double * U, struct U_var * U_L, struct U_var * U_R, double z_sL, double z_sR, double z_sL_out, double z_sR_out, double area_L, double area_R)
{   
    double z_gL=1-z_sL;
    double z_gR=1-z_sR;
    const double gama_g = config[106], gama_s = config[6];
    double eps = config[4];
    double z_s = area_L*z_sL+area_R*z_sR;
    double z_g = 1.0-z_s;
    double U1=U[0], U2=U[1], U3=U[3], U4=U[4], U5=U[5], U6=U[7];
    double rho_gR = U1/z_g;
    double p_gR  = (U3/z_g - 0.5*rho_gR*pow(U2/U1,2))*(gama_g-1.0);
    double rho_s  = U4/z_s;
    double u_s   = U5/U4;
    U_L->v_g = U[2]/U1;
    U_R->v_g = U_L->v_g;	
    U_L->v_s = U[6]/U4;
    U_R->v_s = U_L->v_s;
    U_L->z_s = z_sL;
    U_R->z_s = z_sR;
    double fun[2], dfun[2][2], x_star[2];
    int it_max = 10, it_max1= 2, k = 0;
    double err2 = 1e50;
    while (k<it_max && err2>eps && fabs(z_sL-z_sR)>eps) {			
	fun[0] =  0.5e0 * pow(U2 - (U1 + U4) * u_s, 0.2e1) * area_L * area_L * pow(-area_R * z_gR * rho_gR + U1, -0.2e1) + gama_g * ((U3 + (-0.1e1) * 0.5e0 * area_R * z_gR * rho_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + (-0.1e1) * 0.5e0 * (-area_R * z_gR * rho_gR + U1) * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1)) * (gama_g - 1) - area_R * z_gR * p_gR) / (gama_g - 1) / (-area_R * z_gR * rho_gR + U1) + (-0.1e1) * 0.5e0 * pow(U2 - (U1 + U4) * u_s, 0.2e1) * pow(z_gR, -0.2e1) * pow(rho_gR, -0.2e1) - gama_g * p_gR / (gama_g - 1) / rho_gR;
	fun[1] = ((U3 + (-0.1e1) * 0.5e0 * area_R * z_gR * rho_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + (-0.1e1) * 0.5e0 * (-area_R * z_gR * rho_gR + U1) * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1)) * (gama_g - 1) - area_R * z_gR * p_gR) / area_L / z_gL / pow((-area_R * z_gR * rho_gR + U1) / area_L / z_gL, gama_g) - p_gR / pow(rho_gR, gama_g);
	dfun[0][0] = ((-0.1e1) * 0.5e0 * area_R * z_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + 0.10e1 * area_R * ((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s) * (U2 - (U1 + U4) * u_s) / rho_gR + 0.5e0 * area_R * z_gR * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1) + (-0.1e1) * 0.10e1 * ((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s) * (U2 - (U1 + U4) * u_s) * area_L * area_R * z_gR / (-area_R * z_gR * rho_gR + U1)) * (gama_g - 1) / area_L / z_gL / pow((-area_R * z_gR * rho_gR + U1) / area_L / z_gL, gama_g) + ((U3 + (-0.1e1) * 0.5e0 * area_R * z_gR * rho_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + (-0.1e1) * 0.5e0 * (-area_R * z_gR * rho_gR + U1) * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1)) * (gama_g - 1) - area_R * z_gR * p_gR) * gama_g * area_R * z_gR / area_L / z_gL / pow((-area_R * z_gR * rho_gR + U1) / area_L / z_gL, gama_g) / (-area_R * z_gR * rho_gR + U1) + p_gR * gama_g / pow(rho_gR, gama_g) / rho_gR;
	dfun[0][1] = -area_R * z_gR / area_L / z_gL / pow((-area_R * z_gR * rho_gR + U1) / area_L / z_gL, gama_g) - 0.1e1 / pow(rho_gR, gama_g);
	dfun[1][0] = 0.10e1 * pow(U2 - (U1 + U4) * u_s, 0.2e1) * area_L * area_L * area_R * z_gR * pow(-area_R * z_gR * rho_gR + U1, -0.3e1) + gama_g * ((-0.1e1) * 0.5e0 * area_R * z_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + 0.10e1 * area_R * ((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s) * (U2 - (U1 + U4) * u_s) / rho_gR + 0.5e0 * area_R * z_gR * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1) + (-0.1e1) * 0.10e1 * ((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s) * (U2 - (U1 + U4) * u_s) * area_L * area_R * z_gR / (-area_R * z_gR * rho_gR + U1)) / (-area_R * z_gR * rho_gR + U1) + gama_g * ((U3 + (-0.1e1) * 0.5e0 * area_R * z_gR * rho_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + (-0.1e1) * 0.5e0 * (-area_R * z_gR * rho_gR + U1) * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1)) * (gama_g - 0.1e1) - area_R * z_gR * p_gR) * area_R * z_gR / (gama_g - 0.1e1) * pow(-area_R * z_gR * rho_gR + U1, -0.2e1) + 0.10e1 * pow(U2 - (U1 + U4) * u_s, 0.2e1) * pow(z_gR, -0.2e1) * pow(rho_gR, -0.3e1) + p_gR * gama_g / (gama_g - 0.1e1) * pow(rho_gR, -0.2e1);
	dfun[1][1] = -gama_g * area_R * z_gR / (gama_g - 1) / (-area_R * z_gR * rho_gR + U1) - gama_g / (gama_g - 1) / rho_gR;
	x_star[0] = rho_gR;
	x_star[1] = p_gR;
	NewtonRapshon_matrix2(x_star, &err2, fun, dfun[0], eps);
//        printf("\nRIeq_err:%lf! %lf, %lf, %lf, %lf\n",err2,z_sL,z_sR,rho_gR,p_gR);
//       printf("\ndfun:%lf, %lf, %lf, %lf\n",dfun[0][0],dfun[1][0],dfun[0][1],dfun[1][1]);
   	rho_gR=fmax(x_star[0],eps);
	rho_gR=fmin(rho_gR,U1/area_R/z_gR-eps);
	p_gR =fmax(x_star[1],eps);
	k=k+1;
	if(k==it_max1)
	    {rho_gR = U1/z_gL;
	    p_gR  = (U3/z_g - 0.5*rho_gR*pow(U2/U1,2))*(gama_g-1.0);}
	if(k==it_max1*2)
	    {rho_gR = U1/z_gR;
	    p_gR  = (U3/z_g - 0.5*rho_gR*pow(U2/U1,2))*(gama_g-1.0);}
	if(k==it_max1*3)
	    {rho_gR = (U2-(U1+U4)*u_s)/z_gR/(U2/U1-u_s);
	    p_gR  = (U3/z_g - 0.5*rho_gR*pow(U2/U1,2))*(gama_g-1.0);}
    }
/*
    if (k>=it_max)
{        printf("\nRIeq_err:%lf! %lf, %lf, %lf, %lf\n",err2,z_sL,z_sR,rho_gR,p_gR);
        printf("\ndfun:%lf, %lf, %lf, %lf\n",dfun[0][0],dfun[1][0],dfun[0][1],dfun[1][1]);
} */    
    U_R->rho_g = rho_gR;
    U_L->rho_g = (U1-area_R*z_gR*rho_gR)/area_L/z_gL;
    U_R->u_g = (U2-(U1+U4)*u_s)/(z_gR*rho_gR)+u_s;
    U_L->u_g = (U2-(U1+U4)*u_s)/(U1-area_R*z_gR*rho_gR)*area_L+u_s;
    
    U_R->p_g = p_gR;
    U_L->p_g =((U3 + (-0.1e1) * 0.5e0 * area_R * z_gR * rho_gR * pow((U2 - (U1 + U4) * u_s) / z_gR / rho_gR + u_s, 0.2e1) + (-0.1e1) * 0.5e0 * (-area_R * z_gR * rho_gR + U1) * pow((U2 - (U1 + U4) * u_s) * area_L / (-area_R * z_gR * rho_gR + U1) + u_s, 0.2e1)) * (gama_g - 1) - area_R * z_gR * p_gR) / area_L / z_gL;
    U_R->p_s = (((z_gL * U_L->rho_g + (-0.1e1) * 0.1e1 * z_gR * rho_gR) * u_s * u_s + ((-0.1e1) * 0.2e1 * U_L->rho_g * z_gL * U_L->u_g + 0.2e1 * rho_gR * z_gR * U_R->u_g) * u_s + (U_L->rho_g * U_L->u_g * U_L->u_g + U_L->p_g) * z_gL + ((-0.1e1) * 0.1e1 * rho_gR * U_R->u_g * U_R->u_g + (-0.1e1) * 0.1e1 * p_gR) * z_gR) * area_L + ((-0.1e1) * 0.5e0 * gama_s * rho_s + 0.5e0 * rho_s) * u_s * u_s + U6 * gama_s + (-0.1e1) * 0.1e1 * U6) / z_sR / (area_L + area_R);
    U_L->p_s = ((U6 + (-0.1e1) * 0.5e0 * rho_s * u_s * u_s) * (gama_s - 1) - area_R * (((z_gL * U_L->rho_g + (-0.1e1) * 0.1e1 * z_gR * rho_gR) * u_s * u_s + ((-0.1e1) * 0.2e1 * U_L->rho_g * z_gL * U_L->u_g + 0.2e1 * rho_gR * z_gR * U_R->u_g) * u_s + (U_L->rho_g * U_L->u_g * U_L->u_g + U_L->p_g) * z_gL + ((-0.1e1) * 0.1e1 * rho_gR * U_R->u_g * U_R->u_g + (-0.1e1) * 0.1e1 * p_gR) * z_gR) * area_L + ((-0.1e1) * 0.5e0 * gama_s * rho_s + 0.5e0 * rho_s) * u_s * u_s + U6 * gama_s + (-0.1e1) * 0.1e1 * U6) / (area_L + area_R)) / area_L / z_sL;
    
    U_L->rho_s= rho_s;
    U_R->rho_s= rho_s;
    U_L->u_s = u_s;
    U_R->u_s = u_s;
    
    struct RI_var RI_L, RI_R;
    U2RI_cal(U_L, &RI_L);
    RI2U_cal(U_L, &RI_L, z_sL_out, U_L->rho_g);
    U2RI_cal(U_R, &RI_R);
    RI2U_cal(U_R, &RI_R, z_sR_out, U_R->rho_g);
    
    U_L->U_rho_g = (1.0-z_sL_out)*U_L->rho_g;
    U_L->U_u_g  = U_L->U_rho_g*U_L->u_g;
    U_L->U_v_g  = U_L->U_rho_g*U_L->v_g;
    U_L->U_e_g  = U_L->U_rho_g*(U_L->p_g/U_L->rho_g/(gama_g-1.0)+0.5*U_L->u_g*U_L->u_g+0.5*U_L->v_g*U_L->v_g);
    U_L->U_rho_s = z_sL_out*U_L->rho_s;
    U_L->U_u_s  = U_L->U_rho_s*U_L->u_s;
    U_L->U_v_s  = U_L->U_rho_s*U_L->v_s;
    U_L->U_e_s  = U_L->U_rho_s*(U_L->p_s/U_L->rho_s/(gama_s-1.0)+0.5*U_L->u_s*U_L->u_s+0.5*U_L->v_s*U_L->v_s);
    U_R->U_rho_g = (1.0-z_sR_out)*U_R->rho_g;
    U_R->U_u_g  = U_R->U_rho_g*U_R->u_g;
    U_R->U_v_g  = U_R->U_rho_g*U_R->v_g;
    U_R->U_e_g  = U_R->U_rho_g*(U_R->p_g/U_R->rho_g/(gama_g-1.0)+0.5*U_R->u_g*U_R->u_g+0.5*U_R->v_g*U_R->v_g);
    U_R->U_rho_s = z_sR_out*U_R->rho_s;
    U_R->U_u_s  = U_R->U_rho_s*U_R->u_s;
    U_R->U_v_s  = U_R->U_rho_s*U_R->v_s;
    U_R->U_e_s  = U_R->U_rho_s*(U_R->p_s/U_R->rho_s/(gama_s-1.0)+0.5*U_R->u_s*U_R->u_s+0.5*U_R->v_s*U_R->v_s);	
}
