#include <math.h>
#include <stdio.h>

/*
 * EXACT RIEMANN SOLVER FOR THE EULER EQUATIONS
 * Purpose: to solve the Riemann problem exactly, for the time dependent one dimensional Euler equations for an ideal gas
 * Input:
 * Output: gap
 * Programer: E. F. Toro
 * Last revision: February 1st 1999
 * Theory is found in Chapter 4 of Reference [1]
 * [1] Toro, E. F., "Riemann Solvers and Numerical Methods for Fluid Dynamics" Springer-Verlag, Second Edition, 1999
 * This program is part of NUMERICA ——
 * A Library of Source Codes for Teaching, Research and Applications, by E. F. Toro
 * Published by NUMERITEK LTD
 */

double Riemann_solver_exact(double * U_star, double * P_star,double gamma, double U_l, double U_r, double P_l, double P_r, double c_l, double c_r, int * CRW, double eps, double tol, int N)
{
    int n = 0;		
    double gap;//relative pressure change
	
    double P_int,U_int;//P_star,U_star
    double P_int_save;
    double f_R,f_L,df_R,df_L;

    double RHO_r=gamma * P_r/c_r/c_r;
    double RHO_l=gamma * P_l/c_l/c_l;
	
    double g1=(gamma -1.0);
    double g2=(gamma+1.0);
    double g3=2.0*gamma/(gamma-1.0);
    double g4=2.0/(gamma-1.0);
    double g5=2.0/(gamma+1.0);
    double g6=(gamma-1.0)/(gamma+1.0);
    double g7=(gamma-1.0)/2.0;
    double g8=gamma-1.0;

    double A_L=2.0/g2/RHO_l;
    double A_R=2.0/g2/RHO_r;
    double B_L=g6*P_l;
    double B_R=g6*P_r;

  //======set the approximate value of p_star================================
  P_int  = pow( (c_l + c_r -  0.5*g8*(U_r-U_l)) / (c_l/pow(P_l,1/g3)+c_r/pow(P_r,1/g3)) , g3);

  //=======THE NEWTON ITERATION=======
  while((n < N))
  {
	  P_int_save=P_int;

	  if(P_int > P_l)
    {
		f_L=(P_int - P_l)*pow(A_L/(P_int+B_L),0.5);
		df_L=pow(A_L/(P_int+B_L),0.5)-0.5*(P_int - P_l)*pow(A_L,0.5)/pow(P_int+B_L,1.5);
    }
    else
    {
		f_L=2.0*c_l/g8*(pow(P_int/P_l,1.0/g3)-1.0);
		df_L=c_l/gamma/P_l*pow(P_int/P_l,1.0/g3-1.0);
	}
    if(P_int > P_r)
    {
		f_R=(P_int - P_r)*pow(A_R/(P_int+B_R),0.5);
		df_R=pow(A_R/(P_int+B_R),0.5)-0.5*(P_int - P_r)*pow(A_R,0.5)/pow(P_int+B_R,1.5);
	}
    else
    {
		f_R=2.0*c_r/g8*(pow(P_int/P_r,1.0/g3)-1.0);
		df_R=c_r/gamma/P_r*pow(P_int/P_r,1.0/g3-1.0);
	}

	if(gap < tol)  break;
	
    P_int=P_int - (f_L - f_R + U_r - U_l)/(df_L-df_R);

    gap = 0.5*fabs(P_int - P_int_save) / (P_int + P_int_save) ;
	++n;
  }

  //===============================
	
  if(P_int > P_l-eps)
      CRW[0]=0;
  else
      CRW[0]=1;
  if(P_int > P_r+eps)
      CRW[1]=0;
  else
      CRW[1]=1;
  
  U_int = 0.5*(U_l+U_r)+ 0.5 *(f_R-f_L);

  *P_star = P_int;
  *U_star = U_int;
  
  return gap;
}
