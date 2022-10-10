function value_start(Ncell)
if nargin < 1
    Ncell = 3000   % Number of computing cells in r direction
elseif strcmpi(Ncell,{'INPUT'})
    Ncell = input('Please input the Ncell number: (Default 3000 6000) ')
else
    Ncell
    if ischar(Ncell)
        error("Not a string 'INPUT' was entered to represent the input!");
    elseif Ncell ~= fix(Ncell) || Ncell <= 0
        error("Not a positive integer was entered to represent variable 'Ncell'!")
    end
end

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
%UR0=-200.0
UR0=-10
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
for i=(Ncell3+1):Ncell
    rho(i)=DL1;
    u(i)  =UL1;
    p(i)  =PL1;
    phi(i)=1.0;
end

fid = fopen('RHO.dat','wt');
fprintf(fid,'%12.10g\t',rho(2:Ncell));
fprintf(fid,'\n');
fclose(fid);

fid = fopen('U.dat','wt');
fprintf(fid,'%12.10g\t',u(2:Ncell));
fprintf(fid,'\n');
fclose(fid);

fid = fopen('P.dat','wt');
fprintf(fid,'%12.10g\t',p(2:Ncell));
fprintf(fid,'\n');
fclose(fid);

fid = fopen('PHI.dat','wt');
fprintf(fid,'%12.10g\t',phi(2:Ncell));
fprintf(fid,'\n');
fclose(fid);

%Timeout=0.23;     % Output time
Timeout=0.8;
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
fprintf(fid,'3\t%i\n',Ncell-1);
fprintf(fid,'4\t%g\n',eps);
fprintf(fid,'5\t%i\n',step);
fprintf(fid,'6\t%g\n',GAMMAL);
fprintf(fid,'7\t%g\n',CFL);
fprintf(fid,'10\t%g\n',Domlen/Ncell);
fprintf(fid,'11\t%g\n',0.5*pi/Tcell);
fprintf(fid,'13\t%i\n',Ncell-1);
fprintf(fid,'14\t%i\n',Tcell);
fprintf(fid,'17\t-24\n');
fprintf(fid,'20\t%i\n',Tcell);
fprintf(fid,'106\t%g\n',GAMMAR);
fclose(fid);
end
