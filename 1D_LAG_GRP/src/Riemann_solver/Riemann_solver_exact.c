#include <math.h>
#include <stdio.h>


/**
 * @brief      EXACT RIEMANN SOLVER FOR THE EULER EQUATIONS
 * @details    
 * @param[in]  U_star, P_star Velocity and Pressure values in the star region
 *             gamma          specific heat ratio
 */


double Riemann_solver_exact(double * U_star, double * P_star, double gamma, double u_L, double u_R, double p_L, double p_R, double c_L, double c_R, int * CRW, double eps, double tol, int N)
{
  //double zeta_l, zeta_r;
  double mu, nu, sigma;
  //double c_L, c_R;
  double delta_p, u_LR, u_RL;
  double k1, k3, p_INT, p_INT0, u_INT;
  double rho_star_l, rho_star_r;
  double v_L, v_R, gap;
  double temp1, temp2, temp3;
  double dbg = 2;
  int n = 0;

  mu = (gamma-1.0) / (2.0*gamma);
  nu = (gamma+1.0) / (2.0*gamma);
  sigma = (gamma - 1.0) / (gamma + 1.0);

  //c_L = sqrt(gamma * p_L / rho_l);
  //c_R = sqrt(gamma * p_R / rho_r);
  //rho_l = 1.0 / rho_l;
  //rho_r = 1.0 / rho_r;

  //zeta_l = pow(p_L, (gamma-1.0) / (2.0*gamma));
  //zeta_r = pow(p_R, (gamma-1.0) / (2.0*gamma));

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
    CRW[1] = 0;
  else
    CRW[1] = 1;
  if(u_RL > u_L-eps)
    CRW[0] = 1;
  else
    CRW[0] = 0;

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
    //dbg = pow(p_INT/p_R, mu);
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
