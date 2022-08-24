% Sod's shock tube problem - x-direction
line=100;
column=10;
interface=50;

rho_1= 1
u_1  = 0
v_1  = 0
p_1  = 1
rho_2= 0.125
u_2  = 0
v_2  = 0
p_2  = 0.1

fid = fopen('RHO.dat','wt');
for j=1:interface
rho=rho_1*ones(column,1);
fprintf(fid,'%12.10g\t',rho);
fprintf(fid,'\n');
end
rho=rho_2*ones(column,1);
for j=(interface+1):line
fprintf(fid,'%12.10g\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('U.dat','wt');
u=u_1*ones(column,1);
for j=1:interface
fprintf(fid,'%12.10g\t',u);
fprintf(fid,'\n');
end
u=u_2*ones(column,1);
for j=(interface+1):line
fprintf(fid,'%12.10g\t',u);
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('V.dat','wt');
v=v_1*ones(column,1);
for j=1:interface
fprintf(fid,'%12.10g\t',v);
fprintf(fid,'\n');
end
v=v_2*ones(column,1);
for j=(interface+1):line
fprintf(fid,'%12.10g\t',v);
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('P.dat','wt');
p=p_1*ones(column,1);
for j=1:interface
fprintf(fid,'%12.10g\t',p);
fprintf(fid,'\n');
end
p=p_2*ones(column,1);
for j=(interface+1):line
fprintf(fid,'%12.10g\t',p);
fprintf(fid,'\n');
end
fclose(fid);
