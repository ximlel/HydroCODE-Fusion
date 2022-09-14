
Status AcousticSLagTangent(double &DtP,double &DtU,double DU,double DP,double D,double U,double C_star,double r)
//GRP solver for tangential case. Lagrangian version(moving mesh) for cylindrical case
{
	DtP= -D*C_star*C_star*(DU+(M-1)*U/r);
	DtU= -DP/D;
	return OK;
}

Status AcousticSLag(double &DtDL,double &DtDR,double &DtU,double &DtP,double DUL,double DUR,double DPL,double DPR,
					double D,double U,double C_star,double r)
//GRP solver for acoustic case. Lagrangian version(moving mesh) for cylindrical case
{
	DtP= 0.5*C_star*(DPR-DPL)-D*C_star*C_star*(0.5*(DUL+DUR)+(M-1)*U/r);
	DtU=-0.5*(DPR+DPL)/D+0.5*C_star*(DUR-DUL);
	DtDL=DtP/C_star/C_star;
	DtDR=DtDL;
	return OK;
}

Status GRPsolverSLag(double &DtDL,double &DtDR,double &DtU,double &DtP,double UM,double PM,
					 double DL,double DR,double UL,double UR,double PL,double PR,
					 double DDL,double DDR,double DUL,double DUR,double DPL,double DPR,
					 double TDSL,double TDSR,double DpsiL,double DphiR,double r,double GammaL,double GammaR)
//GRP solver for non acoustic case. Lagrangian version(moving mesh) cylindrical case
{
	double aL,bL,dL,aR,bR,dR,CL,CR,C_starL,C_starR,D,U,P,theta,phic,phi1,phi2,phi3,sigmaL,sigmaR,musL,musR;
	musL=(GammaL-1.)/(GammaL+1.);
	musR=(GammaR-1.)/(GammaR+1.);
	CL=sqrt(GammaL*PL/DL);
	CR=sqrt(GammaR*PR/DR);
	if(PM<=PL)//Left rarefaction wave
		{	// middle left state
			D=DL*pow(PM/PL,1./GammaL);
			U=UM;
			P=PM;
			C_starL=sqrt(GammaL*P/D);
			theta=C_starL/CL;
			aL=1.;
			bL=1./(D*C_starL);
			if(fabs(GammaL-5./3.)<EPS)
				phic=-2.*(3.*C_starL*log(theta)+(UL+2.*CL/(GammaL-1.))*(1.-theta));
			else if(fabs(GammaL-3.)<EPS)
				phic=CL-C_starL+(UL+2.*CL/(GammaL-1.))*log(theta);
			else
				phic=(musL-1.)*C_starL/(musL*(4.*musL-1.))*(1.-pow(theta,(1.-4.*musL)/(2.*musL)))+(UL+2.*CL/(GammaL-1.))/(2.*musL-1.)*(1.-pow(theta,(1.-2*musL)/(2.*musL)));
			dL=((1.+musL)/(1.+2.*musL)*pow(theta,0.5/musL)+musL/(1.+2.*musL)*pow(theta,(1.+musL)/musL))*TDSL
				-pow(theta,0.5/musL)*CL*(DpsiL+(M-1.)/(2.*r)*UL)+(M-1.)/(2.*r)*C_starL*(phic-U);
		}
	else//Left shock wave
		{   //middle Left state(behind shock)
			D=DL*(PM/PL+(GammaL-1.)/(GammaL+1.))/(PM/PL*(GammaL-1.)/(GammaL+1.)+1.);
			U=UM;
			P=PM;
			C_starL=sqrt(GammaL*P/D);
			phi1=0.5*sqrt((1.-musL)/(DL*(PM+musL*PL)))*(PM+PL*(1.+2.*musL))/(PM+musL*PL);
			phi2=-0.5*sqrt((1.-musL)/(DL*(PM+musL*PL)))*(PM*(2.+musL)+musL*PL)/(PM+musL*PL);
			phi3=-0.5*(PM-PL)/DL*sqrt((1.-musL)/(DL*(PM+musL*PL)));
			sigmaL=(D*U-DL*UL)/(D-DL);
			aL=1.-D*(sigmaL-U)*phi1;
			bL=phi1-(sigmaL-U)/(D*C_starL*C_starL);
			dL=(-(sigmaL-UL)*phi2-1./DL)*DPL+(sigmaL-UL+DL*phi3+CL*CL*DL*phi2)*DUL-(sigmaL-UL)*phi3*DDL+(DL*UL*CL*CL*phi2+DL*UL*phi3+(sigmaL-U)*U)*(M-1.)/r;
		}

	if(PM<=PR)//Right rarefaction wave
		{	//middle right state
			D=DR*pow(PM/PR,1./GammaR);
			U=UM;
			P=PM;
			C_starR=sqrt(GammaR*P/D);
			aR=1.;
			bR=-1./(D*C_starR);
			theta=C_starR/CR;
			if(fabs(GammaR-5./3.)<EPS)
				phic=-2.*(3.*C_starR*log(theta)-(UR-2.*CR/(GammaR-1.))*(1.-theta));
			else if(fabs(GammaR-3.)<EPS)
				phic=CR-C_starR-(UR-2.*CR/(GammaR-1.))*log(theta);
			else
				phic=(musR-1.)*C_starR/(musR*(4.*musR-1.))*(1.-pow(theta,(1.-4.*musR)/(2.*musR)))-(UR-2.*CR/(GammaR-1.))/(2.*musR-1.)*(1.-pow(theta,(1.-2*musR)/(2.*musR)));
			dR=((1.+musR)/(1.+2.*musR)*pow(theta,0.5/musR)+musR/(1.+2.*musR)*pow(theta,(1.+musR)/musR))*TDSR
				+pow(theta,0.5/musR)*CR*(DphiR+(M-1.)/(2.*r)*UR)+(M-1.)/(2.*r)*C_starR*(phic+U);
		}
	else//Right shock wave
		{   //middle Right state (behind shock)
			D=DR*(PM/PR+(GammaR-1.)/(GammaR+1.))/(PM/PR*(GammaR-1.)/(GammaR+1.)+1.);
			U=UM;
			P=PM;
			C_starR=sqrt(GammaR*P/D);
			phi1=0.5*sqrt((1.-musR)/(DR*(PM+musR*PR)))*(PM+PR*(1.+2.*musR))/(PM+musR*PR);
			phi2=-0.5*sqrt((1.-musR)/(DR*(PM+musR*PR)))*(PM*(2.+musR)+musR*PR)/(PM+musR*PR);
			phi3=-0.5*(PM-PR)/DR*sqrt((1.-musR)/(DR*(PM+musR*PR)));
			sigmaR=(D*U-DR*UR)/(D-DR);
			aR=1.+D*(sigmaR-UM)*phi1;
			bR=-phi1-(sigmaR-UM)/(D*C_starR*C_starR);
			dR=((sigmaR-UR)*phi2-1./DR)*DPR+(sigmaR-UR-DR*phi3-CR*CR*DR*phi2)*DUR+(sigmaR-UR)*phi3*DDR-(DR*UR*CR*CR*phi2+DR*UR*phi3-(sigmaR-U)*U)*(M-1.)/r;
		}
	DtU=(dL*bR-dR*bL)/(aL*bR-aR*bL);
	DtP=(dL*aR-dR*aL)/(bL*aR-bR*aL);
	DtDL=DtP/C_starL/C_starL;
	DtDR=DtP/C_starR/C_starR;
	return OK;
}
#endif
