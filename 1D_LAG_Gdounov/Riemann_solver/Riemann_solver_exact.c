#include <math.h>
#include <stdio.h>



void Riemann_solver_exact(double * U_star, double * P_star, double gamma, double U_l, double U_r, double P_l, double P_r, double c_l, double c_r, int * CRW, double tol, int N)
{
  //double zeta_l, zeta_r;
  double mu, nu, sigma;
  //double c_l, c_r;
  double delta_P, U_LR, U_RL;
  double k1, k3, P_int, U_int;
  double rho_star_l, rho_star_r;
  double U_temp_l, U_temp_r, gap;
  double temp1, temp2, temp3;
  double dbg = 2;
  int n = 0;

  mu = (gamma-1.0) / (2.0*gamma);
  nu = (gamma+1.0) / (2.0*gamma);
  sigma = (gamma - 1.0) / (gamma + 1.0);

  //c_l = sqrt(gamma * P_l / rho_l);
  //c_r = sqrt(gamma * P_r / rho_r);
  //rho_l = 1.0 / rho_l;
  //rho_r = 1.0 / rho_r;

  //zeta_l = pow(P_l, (gamma-1.0) / (2.0*gamma));
  //zeta_r = pow(P_r, (gamma-1.0) / (2.0*gamma));

  //=====find out the kind of 1-wave and 3-wave
  if(P_r > P_l)
  {
    delta_P = P_r - P_l;
    U_LR = sqrt(1.0 + nu*delta_P/P_l);
    U_LR = delta_P * c_l / gamma / P_l / U_LR;
    U_LR = U_l - U_LR;
  }
  else
  {
    U_LR = pow(P_r/P_l, mu) - 1.0;
    U_LR = 2.0 * c_l * U_LR / (gamma-1.0);
    U_LR = U_l - U_LR;
  }
  if(P_l > P_r)
  {
    delta_P = P_l - P_r;
    U_RL = sqrt(1.0 + nu*delta_P/P_r);
    U_RL = delta_P * c_r / gamma / P_r / U_RL;
    U_RL = U_r + U_RL;
  }
  else
  {
    U_RL = pow(P_l/P_r, mu) - 1.0;
    U_RL = 2.0 * c_r * U_RL / (gamma-1.0);
    U_RL = U_r + U_RL;
  }
  if(U_LR > U_r+tol)
    CRW[1] = 0;
  else
    CRW[1] = 1;
  if(U_RL > U_l-tol)
    CRW[0] = 1;
  else
    CRW[0] = 0;

  //======one step of the Newton ietration to get the intersection point of L1 and L3====
  k1 = c_l / P_l / gamma;
  k1 = -k1;
  k3 = c_r / P_r / gamma;
  P_int = (k1*P_l - k3*P_r - U_l + U_r) / (k1 - k3);

  //=======compute the gap between U^n_R and U^n_L(see Appendix C)=======
  if(P_int > P_l)
  {
    delta_P = P_int - P_l;
    U_temp_l = sqrt(1.0 + nu*delta_P/P_l);
    U_temp_l = delta_P * c_l / gamma / P_l / U_temp_l;
    U_temp_l = U_l - U_temp_l;
  }
  else
  {
    U_temp_l = pow(P_int/P_l, mu) - 1.0;
    U_temp_l = 2.0 * c_l * U_temp_l / (gamma-1.0);
    U_temp_l = U_l - U_temp_l;
  }
  if(P_int > P_r)
  {
    delta_P = P_int - P_r;
    U_temp_r = sqrt(1.0 + nu*delta_P/P_r);
    U_temp_r = delta_P * c_r / gamma / P_r / U_temp_r;
    U_temp_r = U_r + U_temp_r;
  }
  else
  {
    dbg = pow(P_int/P_r, mu);
    U_temp_r = pow(P_int/P_r, mu) - 1.0;
    U_temp_r = 2.0 * c_r * U_temp_r / (gamma-1.0);
    U_temp_r = U_r + U_temp_r;
  }
  gap = fabs(U_temp_l - U_temp_r);


  //=======THE NEWTON ITERATION=======
  while((gap > tol) && (n != N))
  {
    //------compute the slope of the tangent of L1 and L3------
    if(P_int > P_l)
    {
      delta_P = P_int - P_l;
      temp1 = 1.0 / sqrt(1.0 + nu*delta_P/P_l);
      temp2 = c_l / gamma / P_l;
      temp3 = 0.5 * temp2 * nu / P_l;
      k1 = temp3*delta_P*pow(temp1,3.0) - temp2*temp1;
    }
    else
    {
      temp2 = c_l / gamma / P_l;
      temp1 = 1.0 / pow(P_int/P_l, nu);
      k1 = temp1 * temp2;
      k1 = -k1;
    }
    if(P_int > P_r)
    {
      delta_P = P_int - P_r;
      temp1 = 1.0 / sqrt(1.0 + nu*delta_P/P_r);
      temp2 = c_r / gamma / P_r;
      temp3 = 0.5 * temp2 * nu / P_r;
      k3 = temp2*temp1 - temp3*delta_P*pow(temp1,3.0);
    }
    else
    {
      temp2 = c_r / gamma / P_r;
      temp1 = 1.0 / pow(P_int/P_r, nu);
      k3 = temp1 * temp2;
    }

    //------compute the intersetion point of L1 and L3------
    P_int = P_int + (U_temp_r - U_temp_l) / (k1 - k3);

    //------the gap------
    ++n;
    if(P_int > P_l)
    {
      delta_P = P_int - P_l;
      U_temp_l = sqrt(1.0 + nu*delta_P/P_l);
      U_temp_l = delta_P * c_l / gamma / P_l / U_temp_l;
      U_temp_l = U_l - U_temp_l;
    }
    else
    {
      U_temp_l = pow(P_int/P_l, mu) - 1.0;
      U_temp_l = 2.0 * c_l * U_temp_l / (gamma-1.0);
      U_temp_l = U_l - U_temp_l;
    }
    if(P_int > P_r)
    {
      delta_P = P_int - P_r;
      U_temp_r = sqrt(1.0 + nu*delta_P/P_r);
      U_temp_r = delta_P * c_r / gamma / P_r / U_temp_r;
      U_temp_r = U_r + U_temp_r;
    }
    else
    {
      U_temp_r = pow(P_int/P_r, mu) - 1.0;
      U_temp_r = 2.0 * c_r * U_temp_r / (gamma-1.0);
      U_temp_r = U_r + U_temp_r;
    }

    gap = fabs(U_temp_l - U_temp_r);
  }


  U_int = k1*(U_temp_r-U_temp_l)/(k1-k3)+U_temp_l;

  *P_star = P_int;
  *U_star = U_int;
}
