line=1;
column=1000;
interface=500;

rho_1=1.0
u_1=0
p_1=25.0
phi_1=0
rho_2=1.0
u_2=0
p_2=0.01
phi_2=1

rho=zeros(column,1);
for i=1:interface
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
for i=1:interface
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
for i=1:interface
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
for i=1:interface
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
fclose(fid);
