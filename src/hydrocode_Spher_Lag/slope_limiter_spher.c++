_Bool slope_limiter_spher(const int m, const int n, const int nt, struct cell_var_stru * CV, struct b_f_var * bfv_L, struct b_f_var * bfv_R,
				 struct b_f_var * bfv_D, struct b_f_var * bfv_U, _Bool find_bound_x, const _Bool Slope, const double t_c)
			/*
			for(i=1;i<Ncell;i++)
				{
					sU=(Umin[i+1] -Umin[i]) /Ddr[i];
					sP=(Pmin[i+1] -Pmin[i]) /Ddr[i];
					sD=(DLmin[i+1]-DRmin[i])/Ddr[i];
					DmD[i]=minmod3(Alpha*(DD[i]-DD[i-1])/dRc[i],sD,Alpha*(DD[i+1]-DD[i])/dRc[i+1]);
					DmU[i]=minmod3(Alpha*(UU[i]-UU[i-1])/dRc[i],sU,Alpha*(UU[i+1]-UU[i])/dRc[i+1]);
					DmP[i]=minmod3(Alpha*(PP[i]-PP[i-1])/dRc[i],sP,Alpha*(PP[i+1]-PP[i])/dRc[i+1]);
				}
			DmD[0]=minmod2((DLmin[1]-DD[0])/dRc[0],DmD[1]);
			DmU[0]=minmod2((Umin[1]-UU[0]) /dRc[0],DmU[1]);
			DmP[0]=minmod2((Pmin[1]-PP[0]) /dRc[0],DmP[1]);
			DmD[Ncell]=minmod2((DLmin[Ncell]-DD[Ncell-1])/dRc[Ncell],DmD[Ncell-1]);
			DmU[Ncell]=minmod2((Umin[Ncell]-UU[Ncell-1]) /dRc[Ncell],DmU[Ncell-1]);
			DmP[Ncell]=minmod2((Pmin[Ncell]-PP[Ncell-1]) /dRc[Ncell],DmP[Ncell-1]);
			*/
			//VIP limiter update
			for(i=1;i<Ncell;i++)
				{
					sU=(Umin[i+1] -Umin[i]) /Ddr[i];
					sP=(Pmin[i+1] -Pmin[i]) /Ddr[i];
					sD=(DLmin[i+1]-DRmin[i])/Ddr[i];
					//sV=0.;
					//sV=VLmin[i]/(0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta));
					sV=UU[i]/(0.5*(Rb[i]+Rb[i+1]));
					Vave[0][0] = UU[i+1];
					Vave[0][1] = 0.;
					Vave[1][0] = UU[i-1];
					Vave[1][1] = 0.;
					Vave[2][0] = UU[i]*cos(dtheta);
					Vave[2][1] = UU[i]*sin(dtheta);
					Vave[3][0] = UU[i]*cos(dtheta);
					Vave[3][1] =-UU[i]*sin(dtheta);
					V0[0] = UU[i];
					V0[1] = 0.;
					Vp1[0] = UU[i]+(0.5*(Rb[i]+Rb[i+1])-RR[i])*sU;
					Vp1[1] = 0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta)*sV;
					Vp2[0] = UU[i]+DdrL[i]*sU;
					Vp2[1] = 0.0;
					Vp3[0] = UU[i]-DdrR[i]*sU;
					Vp3[1] = 0.0;
					VIP_lim = fmin(1.0,     useVIPLimiter(4, Vave, V0, Vp1));
					VIP_lim = fmin(VIP_lim, useVIPLimiter(4, Vave, V0, Vp2));
					VIP_lim = fmin(VIP_lim, useVIPLimiter(4, Vave, V0, Vp3));
					DmU[i]=VIP_lim*sU;
					TmV[i]=VIP_lim*sV;
					if (abs(LIMITER_CONF)==1)
						{
							DmD[i]=minmod3(Alpha*(DD[i]-DD[i-1])/dRc[i],sD,Alpha*(DD[i+1]-DD[i])/dRc[i+1]);
							DmP[i]=minmod3(Alpha*(PP[i]-PP[i-1])/dRc[i],sP,Alpha*(PP[i+1]-PP[i])/dRc[i+1]);
							if(LIMITER_CONF>0)
								DmU[i]=minmod3(Alpha*(UU[i]-UU[i-1])/dRc[i],sU,Alpha*(UU[i+1]-UU[i])/dRc[i+1]);
						}
					else if (abs(LIMITER_CONF)==2)
						{
							DmD[i]=minmod3(Alpha*(DD[i]-DD[i-1])/2./DdrR[i],sD,Alpha*(DD[i+1]-DD[i])/2./DdrL[i]);
							DmP[i]=minmod3(Alpha*(PP[i]-PP[i-1])/2./DdrR[i],sP,Alpha*(PP[i+1]-PP[i])/2./DdrL[i]);
							if(LIMITER_CONF>0)
								DmU[i]=minmod3(Alpha*(UU[i]-UU[i-1])/2./DdrR[i],sU,Alpha*(UU[i+1]-UU[i])/2./DdrL[i]);
						}
					if (LIMITER_CONF>0)
						TmV[i]=minmod2(Alpha*(UU[i]*sin(dtheta))/2./(0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta)),sV);
				}
			// i = 0
			sU=(Umin[1]-UU[0]) /dRc[0];
			sP=(Pmin[1]-PP[0]) /dRc[0];
			sD=(DLmin[1]-DD[0])/dRc[0];
			//sV=0.;
			//sV=VLmin[1]/(0.5*Rb[1]*tan(0.5*dtheta));
			sV=UU[0]/Rb[1];
			Vave[0][0] = UU[1];
			Vave[0][1] = 0.;
			Vave[1][0] = UU[0]*cos(dtheta);
			Vave[1][1] = UU[0]*sin(dtheta);
			Vave[2][0] = UU[0]*cos(dtheta);
			Vave[2][1] =-UU[0]*sin(dtheta);
			V0[0] = UU[0];
			V0[1] = 0.;
			Vp1[0] = UU[0]+(0.5*Rb[1]-RR[0])*sU;
			Vp1[1] = 0.5*Rb[1]*tan(0.5*dtheta)*sV;
			Vp2[0] = UU[0]+DdrL[0]*sU;
			Vp2[1] = 0.;
			VIP_lim = fmin(1.0,     useVIPLimiter(3, Vave, V0, Vp1));
			VIP_lim = fmin(VIP_lim, useVIPLimiter(3, Vave, V0, Vp2));
			DmU[0]=VIP_lim*sU;
			TmV[0]=VIP_lim*sV;
			DmD[0]=minmod2(sD,DmD[1]);
			DmP[0]=minmod2(sP,DmP[1]);
			if (LIMITER_CONF>0)
				{
					DmU[0]=minmod2(sU,DmU[1]);
					TmV[0]=minmod2(sV,TmV[1]);
				}
			// i = Ncell
			sU=(Umin[Ncell]-UU[Ncell-1]) /dRc[Ncell];
			sP=(Pmin[Ncell]-PP[Ncell-1]) /dRc[Ncell];
			sD=(DLmin[Ncell]-DD[Ncell-1])/dRc[Ncell];
			DmD[Ncell]=minmod2(sD,DmD[Ncell-1]);
			DmP[Ncell]=minmod2(sP,DmP[Ncell-1]);
			DmU[Ncell]=minmod2(sU,DmU[Ncell-1]); 
					}//end k
}
