% Almost stationary shock
line=1;
column=100;
interface=50;

rho_1= 4
u_1  = -0.3
p_1  = 5/3
rho_2= 1
u_2  = -1.3
p_2  = 1e-6

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
