/**
 * @file  Riemann_solver_exact_Ben.c
 * @brief This is an exact Riemann solver in Ben-Artzi's book.
 */

#include <math.h>
#include <stdio.h>
#include <stdbool.h>


/**
 * @brief EXACT RIEMANN SOLVER FOR A γ-Law Gas
 * @details The purpose of this function is to solve the Riemann problem exactly,
 *          for the time dependent one dimensional Euler equations for a γ-law gas.
 * @param[out] U_star, P_star: Velocity/Pressure in star region.
 * @param[in]  u_L, p_L, c_L: Initial Velocity/Pressure/sound_speed on left  state.
 * @param[in]  u_R, p_R, c_R: Initial Velocity/Pressure/sound_speed on right state.
 * @param[in]  gamma: Ratio of specific heats.
 * @param[out] CRW: Centred Rarefaction Wave (CRW) Indicator of left and right waves.
 *                  - true: CRW
 *                  - false: Shock wave
 * @param[in]  eps: The largest value can be seen as zero.
 * @param[in]  tol: Condition value of 'gap' at the end of the iteration.
 * @param[in]  N:   Maximum iteration step.
 * @return \b gap: Relative pressure change after the last iteration.
 * @par  Reference
 *       Theory is found in Appendix C of Reference [1]. \n
 *       [1] M. Ben-Artzi & J. Falcovitz, "Generalized Riemann problems in computational fluid dynamics", 
 *           Cambridge University Press, 2003
 */
double Riemann_solver_exact_Ben(double * U_star, double * P_star, const double gamma,
			    const double u_L, const double u_R, const double p_L, const double p_R,
			    const double c_L, const double c_R, _Bool * CRW,
			    const double eps, const double tol, const int N)
{
  double mu, nu;
  double delta_p, u_LR, u_RL;
  double k1, k3, p_INT, p_INT0, u_INT;
  double v_L, v_R, gap;
  double temp1, temp2, temp3;
  int n = 0;

  mu = (gamma-1.0) / (2.0*gamma);
  nu = (gamma+1.0) / (2.0*gamma);

  //=====find out the kinds of the 1-wave and the 3-wave, page 132 in the GRP book
  //find out where (u_LR,p_R) lies on the curve of LEFT state
  if(p_R > p_L) // (u_LR,p_R) lies on the shock branch of I1
  {
    delta_p = p_R - p_L;
    u_LR = sqrt(1.0 + nu*delta_p/p_L);
    u_LR = delta_p * c_L / gamma / p_L / u_LR;
    u_LR = u_L - u_LR;
  }
  else // (u_LR,p_R) lies on the rarefaction branch of I1
  {
    u_LR = pow(p_R/p_L, mu) - 1.0;
    u_LR = 2.0 * c_L * u_LR / (gamma-1.0);
    u_LR = u_L - u_LR;
  }
  //find out where (u_RL,p_L) lies on the curve of RIGHT state
  if(p_L > p_R) // (u_RL, p_L) lies on the shock branch of I3
  {
    delta_p = p_L - p_R;
    u_RL = sqrt(1.0 + nu*delta_p/p_R);
    u_RL = delta_p * c_R / gamma / p_R / u_RL;
    u_RL = u_R + u_RL;
  }
  else // (u_RL, p_L) lies on the rarefaction branch of I3
  {
    u_RL = pow(p_L/p_R, mu) - 1.0;
    u_RL = 2.0 * c_R * u_RL / (gamma-1.0);
    u_RL = u_R + u_RL;
  }
  if(u_LR > u_R+eps)
    CRW[1] = false;
  else
    CRW[1] = true;
  if(u_RL > u_L-eps)
    CRW[0] = true;
  else
    CRW[0] = false;

  //======one step of the Newton ietration to get the intersection point of I1 and I3====
  k1 = -c_L / p_L / gamma;//the (p,u)-tangent slope on I1 at (u_L,p_L), i.e. [du/dp](p_L)
  k3 =  c_R / p_R / gamma;//the (p,u)-tangent slope on I3 at (u_R,p_R), i.e. [du/dp](p_R)
  //the intersect of (u-u_L)=k1*(p-p_L) and (u-u_R)=k3*(p-p_R)
  p_INT = (k1*p_L - k3*p_R - u_L + u_R) / (k1 - k3);
  if(p_INT < 0)
    p_INT = (p_L<p_R)? p_L : p_R;
  p_INT = 0.5*p_INT;

  //=======compute the gap between U^n_R and U^n_L(see Appendix C)=======
  if(p_INT > p_L)
  {
    delta_p = p_INT - p_L;
    v_L = sqrt(1.0 + nu*delta_p/p_L);
    v_L = delta_p * c_L / gamma / p_L / v_L;
    v_L = u_L - v_L;
  }
  else
  {
    v_L = pow(p_INT/p_L, mu) - 1.0;
    v_L = 2.0 * c_L * v_L / (gamma-1.0);
    v_L = u_L - v_L;
  }
  if(p_INT > p_R)
  {
    delta_p = p_INT - p_R;
    v_R = sqrt(1.0 + nu*delta_p/p_R);
    v_R = delta_p * c_R / gamma / p_R / v_R;
    v_R = u_R + v_R;
  }
  else
  {
    v_R = pow(p_INT/p_R, mu) - 1.0;
    v_R = 2.0 * c_R * v_R / (gamma-1.0);
    v_R = u_R + v_R;
  }
  gap = fabs(v_L - v_R);


  //=======THE NEWTON ITERATION=======
  while((gap > tol) && (n != N))
  {
    //the (p,u)-tangent slope on I1 at (v_L,p_INT), i.e. [du/dp](p_INT)
    if(p_INT > p_L)
    {
      delta_p = p_INT - p_L;
      temp1 = 1.0 / sqrt(1.0 + nu*delta_p/p_L);
      temp2 = c_L / gamma / p_L;
      temp3 = 0.5 * temp2 * nu / p_L;
      k1 = temp3*delta_p*pow(temp1,3.0) - temp2*temp1;
    }
    else
    {
      temp2 = c_L / gamma / p_L;
      temp1 = 1.0 / pow(p_INT/p_L, nu);
      k1 = -temp1 * temp2;
    }
    //the (p,u)-tangent slope on I3 at (v_R,p_INT), i.e. [du/dp](p_INT)
    if(p_INT > p_R)
    {
      delta_p = p_INT - p_R;
      temp1 = 1.0 / sqrt(1.0 + nu*delta_p/p_R);
      temp2 = c_R / gamma / p_R;
      temp3 = 0.5 * temp2 * nu / p_R;
      k3 = temp2*temp1 - temp3*delta_p*pow(temp1,3.0);
    }
    else
    {
      temp2 = c_R / gamma / p_R;
      temp1 = 1.0 / pow(p_INT/p_R, nu);
      k3 = temp1 * temp2;
    }

    //the intersect of (u-u_L)=k1*(p-p_L) and (u-u_R)=k3*(p-p_R)
    p_INT0 = p_INT + (v_R - v_L) / (k1 - k3);
    if(p_INT0 < 0.0)
      p_INT = 0.5*p_INT;
    else
      p_INT = p_INT0;

    //------the gap------
    ++n;
    if(p_INT > p_L)
    {
      delta_p = p_INT - p_L;
      v_L = sqrt(1.0 + nu*delta_p/p_L);
      v_L = delta_p * c_L / gamma / p_L / v_L;
      v_L = u_L - v_L;
    }
    else
    {
      v_L = pow(p_INT/p_L, mu) - 1.0;
      v_L = 2.0 * c_L * v_L / (gamma-1.0);
      v_L = u_L - v_L;
    }
    if(p_INT > p_R)
    {
      delta_p = p_INT - p_R;
      v_R = sqrt(1.0 + nu*delta_p/p_R);
      v_R = delta_p * c_R / gamma / p_R / v_R;
      v_R = u_R + v_R;
    }
    else
    {
      v_R = pow(p_INT/p_R, mu) - 1.0;
      v_R = 2.0 * c_R * v_R / (gamma-1.0);
      v_R = u_R + v_R;
    }

    gap = fabs(v_L - v_R);
  }

  u_INT = k1*(v_R-v_L)/(k1-k3)+v_L;

  *P_star = p_INT;
  *U_star = u_INT;

  return gap;
}
