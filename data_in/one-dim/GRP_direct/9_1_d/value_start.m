% Interacting blast wave problem
line=1;
column=800;
interface1=column/10;
interface2=column/10*9;

rho_1= 1
u_1  = 0
p_1  = 1000
rho_2= 1
u_2  = 0
p_2  = 0.01
rho_3= 1
u_3  = 0
p_3  = 100

rho=zeros(column,1);
for i=1:interface1
    rho(i)=rho_1;
end
for i=(interface1+1):interface2
    rho(i)=rho_2;
end
for i=(interface2+1):column
    rho(i)=rho_3;
end
fid = fopen('RHO.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
for i=1:interface1
    u(i)=u_1;
end
for i=(interface1+1):interface2
    u(i)=u_2;
end
for i=(interface2+1):column
    u(i)=u_3;
end
fid = fopen('U.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',u);
fprintf(fid,'\n');
end
fclose(fid);

p=zeros(column,1);
for i=1:interface1
    p(i)=p_1;
end
for i=(interface1+1):interface2
    p(i)=p_2;
end
for i=(interface2+1):column
    p(i)=p_3;
end
fid = fopen('P.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',p);
fprintf(fid,'\n');
end
fclose(fid);
