column=50;
line=50;
d=0.2;

GAMMA=5.0;
gamma=1.4;

rho_u_f   = @(x,y) -GAMMA/2.0/pi.*y.*exp((1-x.^2-y.^2)/2.0).*(1.0-(gamma-1.0)*GAMMA^2/8.0/gamma/pi^2*exp(1-x.^2-y.^2)).^(1.0/(gamma-1.0));
rho_v_f   = @(x,y)  GAMMA/2.0/pi.*x.*exp((1-x.^2-y.^2)/2.0).*(1.0-(gamma-1.0)*GAMMA^2/8.0/gamma/pi^2*exp(1-x.^2-y.^2)).^(1.0/(gamma-1.0));
rho_f = @(x,y) (1.0-(gamma-1.0)*GAMMA^2/8.0/gamma/pi^2*exp(1-x.^2-y.^2)).^(1.0/(gamma-1.0));
E_f   = @(x,y) (1.0-(gamma-1.0)*GAMMA^2/8.0/gamma/pi^2*exp(1-x.^2-y.^2)).^(gamma/(gamma-1.0))/(gamma-1.0)+0.5*(1.0-(gamma-1.0)*GAMMA^2/8.0/gamma/pi^2*exp(1-x.^2-y.^2)).^(1.0/(gamma-1.0)).*((-GAMMA/2.0/pi.*y.*exp((1-x.^2-y.^2)/2.0)).^2+(GAMMA/2.0/pi.*x.*exp((1-x.^2-y.^2)/2.0)).^2);

center = 25;

RHO = zeros(line,column);
U = zeros(line,column);
V = zeros(line,column);
P = zeros(line,column);
for j=1:column
for i=1:line
       RHO(i,j) = quad2d(rho_f,d*(j-center-1),d*(j-center),d*(i-center-1),d*(i-center),'AbsTol',1e-11);
       U(i,j)   = quad2d(rho_u_f,  d*(j-center-1),d*(j-center),d*(i-center-1),d*(i-center),'AbsTol',1e-11)/RHO(i,j);
       V(i,j)   = quad2d(rho_v_f,  d*(j-center-1),d*(j-center),d*(i-center-1),d*(i-center),'AbsTol',1e-11)/RHO(i,j);
       P(i,j)   = (quad2d(E_f,  d*(j-center-1),d*(j-center),d*(i-center-1),d*(i-center),'AbsTol',1e-11)-0.5*RHO(i,j)*(U(i,j)*U(i,j)+V(i,j)*V(i,j)))*(gamma-1.0);   
end
end
U=U+0.5;
%V=V+1.0;


fid = fopen('RHO.dat','wt');
for j=1:line
fprintf(fid,'%16.11g\t',RHO(j,:));
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('U.dat','wt');
for j=1:line
fprintf(fid,'%16.11g\t',U(j,:));
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('V.dat','wt');
for j=1:line
fprintf(fid,'%16.11g\t',V(j,:));
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('P.dat','wt');
for j=1:line
fprintf(fid,'%16.11g\t',P(j,:));
fprintf(fid,'\n');
end
fclose(fid);
