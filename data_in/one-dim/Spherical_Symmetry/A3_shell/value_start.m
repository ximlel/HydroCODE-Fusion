Ncell=3000;  % Number of computing cells in r direction
Diaph1=10.0;
Diaph2=10.2;
Diaph3=10.25;
Domlen=15.0; % Domain length
Ncell1=round(Ncell*Diaph1/Domlen);
Ncell2=round(Ncell*Diaph2/Domlen);
Ncell3=round(Ncell*Diaph3/Domlen);

DL0=0.00129
DR0=19.237
DL1=0.000129
UL0=0.0
UR0=-200.0
UL1=0.0
PL0=1.01325
PR0=1.01325
PL1=0.101325

rho=zeros(Ncell,1);
u  =zeros(Ncell,1);
p  =zeros(Ncell,1);
phi=zeros(Ncell,1);
for i=1:Ncell1
    rho(i)=DL0;
    u(i)  =UL0;
    p(i)  =PL0;
    phi(i)=1.0;
end
for i=(Ncell1+1):Ncell2
    rho(i)=DR0;
    u(i)  =UL0;
    p(i)  =PR0;
    phi(i)=0.0;
end
for i=(Ncell2+1):Ncell3
    rho(i)=DR0;
    u(i)  =UR0;
    p(i)  =PR0;
    phi(i)=0.0;
end
for i=(Ncell3+1):Domlen
    rho(i)=DL1;
    u(i)  =UL1;
    p(i)  =PL1;
    phi(i)=1.0;
end

fid = fopen('RHO.dat','wt');
for j=1:Tcell
fprintf(fid,'%12.10g\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('U.dat','wt');
for j=1:Tcell
fprintf(fid,'%12.10g\t',u);
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('P.dat','wt');
for j=1:Tcell
fprintf(fid,'%12.10g\t',p);
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('PHI.dat','wt');
for j=1:Tcell
fprintf(fid,'%12.10g\t',phi);
fprintf(fid,'\n');
end
fclose(fid);


Timeout=0.23;     % Output time
D_PLOT_T=0.001;   % Output time interval
time_plot=0:D_PLOT_T:Timeout;
fid = fopen('time_plot.dat','wt');
fprintf(fid,'%12.10g\t',time_plot);
fprintf(fid,'\n');
fclose(fid);


eps=1e-9;
step=25000;
GAMMAL=1.4;       % Ratio of special heats
GAMMAR=3.0;
CFL=0.45;         % CFL condition
Tcell=400;        % Number of computing cells in \theta direction
Tcell_plot=100.0; % Output zoom
Epsilon=1.0;      % r_0=Epsilon*dr

fid = fopen('config.dat','wt');
fprintf(fid,'1\t%g\n',Timeout);
fprintf(fid,'2\t2\n');
fprintf(fid,'4\t%g\n',eps);
fprintf(fid,'5\t%i\n',step);
fprintf(fid,'6\t%g\n',GAMMAL);
fprintf(fid,'7\t%g\n',CFL);
fprintf(fid,'10\t%g\n',Domlen/Ncell);
fprintf(fid,'13\t%i\n',Ncell);
fprintf(fid,'14\t%i\n',Tcell);
fprintf(fid,'106\t%g\n',GAMMAR);
fclose(fid);
