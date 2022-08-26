line=400;
column=400;
line_discon=200;
column_discon=200;

rho_1=1.5
u_1=0
v_1=0
p_1=1.5

rho_2=0.5323
u_2=1.206
v_2=0
p_2=0.3

rho_3=0.138
u_3=1.206
v_3=1.206
p_3=0.029

rho_4=0.5323
u_4=0
v_4=1.206
p_4=0.3


rho=zeros(column,1);
fid = fopen('RHO.dat','wt');
for i=1:column_discon
    rho(i)=rho_3;
end
for i=(column_discon+1):column
    rho(i)=rho_4;
end
for j=1:line_discon
fprintf(fid,'%g\t',rho);
fprintf(fid,'\n');
end
for i=1:column_discon
    rho(i)=rho_2;
end
for i=(column_discon+1):column
    rho(i)=rho_1;
end
for j=(line_discon+1):line
fprintf(fid,'%g\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
fid = fopen('U.dat','wt');
for i=1:column_discon
    u(i)=u_3;
end
for i=(column_discon+1):column
    u(i)=u_4;
end
for j=1:line_discon
fprintf(fid,'%g\t',u);
fprintf(fid,'\n');
end
for i=1:column_discon
    u(i)=u_2;
end
for i=(column_discon+1):column
    u(i)=u_1;
end
for j=(line_discon+1):line
fprintf(fid,'%g\t',u);
fprintf(fid,'\n');
end
fclose(fid);

v=zeros(column,1);
fid = fopen('V.dat','wt');
for i=1:column_discon
    v(i)=v_3;
end
for i=(column_discon+1):column
    v(i)=v_4;
end
for j=1:line_discon
fprintf(fid,'%g\t',v);
fprintf(fid,'\n');
end
for i=1:column_discon
    v(i)=v_2;
end
for i=(column_discon+1):column
    v(i)=v_1;
end
for j=(line_discon+1):line
fprintf(fid,'%g\t',v);
fprintf(fid,'\n');
end
fclose(fid);

p=zeros(column,1);
fid = fopen('P.dat','wt');
for i=1:column_discon
    p(i)=p_3;
end
for i=(column_discon+1):column
    p(i)=p_4;
end
for j=1:line_discon
fprintf(fid,'%g\t',p);
fprintf(fid,'\n');
end
for i=1:column_discon
    p(i)=p_2;
end
for i=(column_discon+1):column
    p(i)=p_1;
end
for j=(line_discon+1):line
fprintf(fid,'%g\t',p);
fprintf(fid,'\n');
end
fclose(fid);
