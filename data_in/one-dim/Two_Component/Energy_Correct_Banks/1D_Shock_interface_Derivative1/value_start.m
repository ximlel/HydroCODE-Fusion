line=1;
column=400;
interface=40;
shock=32;

rho_1=1
u_1=0.5
p_1=1
phi_1=0
rho_2=0.0875
u_2=0.5
p_2=1
phi_2=1
phi_3=0

gamma=1.35;
M_s=1.5;
M=abs(u_1/sqrt(gamma*p_1/rho_1)-M_s);
f=1/(2/(gamma+1)/M/M+(gamma-1)/(gamma+1));
g=2*gamma/(gamma+1)*M*M-(gamma-1)/(gamma+1);
rho_3=rho_1*f
u_3=(1-1/f)*(u_1+sqrt(gamma*p_1/rho_1)*M)+u_1/f
p_3=p_1*g

rho=zeros(column,1);
for i=1:shock
    rho(i)=rho_3;
end
for i=(shock+1):interface
    rho(i)=rho_1;
end
for i=(interface+1):column
    rho(i)=rho_2;
end
fid = fopen('RHO.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
for i=1:shock
    u(i)=u_3;
end
for i=(shock+1):interface
    u(i)=u_1;
end
for i=(interface+1):column
    u(i)=u_2;
end
fid = fopen('U.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',u);
fprintf(fid,'\n');
end
fclose(fid);

v=zeros(column,1);
fid = fopen('V.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',v);
fprintf(fid,'\n');
end
fclose(fid);

p=zeros(column,1);
for i=1:shock
    p(i)=p_3;
end
for i=(shock+1):interface
    p(i)=p_1;
end
for i=(interface+1):column
    p(i)=p_2;
end
fid = fopen('P.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',p);
fprintf(fid,'\n');
end
fclose(fid);

phi=zeros(column,1);
for i=1:shock
    phi(i)=phi_3;
end
for i=(shock+1):interface
    phi(i)=phi_1;
end
for i=(interface+1):column
    phi(i)=phi_2;
end
fid = fopen('PHI.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',phi);
fprintf(fid,'\n');
end
fid = fopen('Z_a.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',phi);
fprintf(fid,'\n');
end
fclose(fid);
