% The Noh shock reflection problem
line=1;
column=100;

rho_1= 1
u_1  = -1
p_1  = 1e-6

rho=zeros(column,1);
for i=1:column
    rho(i)=rho_1;
end
fid = fopen('RHO.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
for i=1:column
    u(i)=u_1;
end
fid = fopen('U.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',u);
fprintf(fid,'\n');
end
fclose(fid);

p=zeros(column,1);
for i=1:column
    p(i)=p_1;
end
fid = fopen('P.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',p);
fprintf(fid,'\n');
end
fclose(fid);
