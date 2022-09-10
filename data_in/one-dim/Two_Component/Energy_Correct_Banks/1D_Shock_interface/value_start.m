function value_start(column, line)
if nargin < 1
    column=1600
    line  =1
elseif strcmpi(column,{'INPUT'})
    column = input('Please input the column number: (Default 1600 640) ')
    line   = input('Please input the line   number: (Default 1 3) ')
else
    column
    if nargin < 2
        line=1
    else
        line
    end
    if isstring(column)
        error("Not a string 'INPUT' was entered to represent the input!");
    elseif column ~= fix(column) || column <= 0
        error("Not a positive integer was entered to represent variable 'column'!")
    elseif line   ~= fix(line)   || line   <= 0
        error("Not a positive integer was entered to represent variable 'line'!")
    end
end

interface=column/2;
shock=column/10;

rho_1=1
u_1=0
p_1=1
phi_1=0
rho_2=1.9
u_2=0
p_2=1
phi_2=1
rho_3=2.7647
u_3=1.4833
p_3=4.4468
phi_3=0

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
fprintf(fid,'%12.10f\t',rho);
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
fprintf(fid,'%12.10f\t',u);
fprintf(fid,'\n');
end
fclose(fid);

v=zeros(column,1);
fid = fopen('V.dat','wt');
for j=1:line
fprintf(fid,'%12.10f\t',v);
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
fprintf(fid,'%12.10f\t',p);
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
fprintf(fid,'%12.10f\t',phi);
fprintf(fid,'\n');
end
fid = fopen('Z_a.dat','wt');
for j=1:line
fprintf(fid,'%12.10f\t',phi);
fprintf(fid,'\n');
end
fclose(fid);


eps=1e-9;
t_all=0.25;
step=50000;
gamma_a = 5.0;
gamma_b = 1.35;
Cv_a = 1.5;
Cv_b = 2.4;
L_x = 1.0/column;
L_y = L_x;

fid = fopen('config.dat','wt');
fprintf(fid,'1\t%g\n',t_all);
fprintf(fid,'2\t2\n');
fprintf(fid,'4\t%g\n',eps);
fprintf(fid,'5\t%i\n',step);
fprintf(fid,'6\t%g\n',gamma_a);
fprintf(fid,'10\t%g\n',L_x);
fprintf(fid,'11\t%g\n',L_y);
fprintf(fid,'13\t%g\n',column);
fprintf(fid,'14\t%g\n',line);
fprintf(fid,'17\t-4\n');
fprintf(fid,'18\t-4\n');
fprintf(fid,'61\t0\n');
fprintf(fid,'70\t0\n');
fprintf(fid,'106\t%g\n',gamma_b);
fprintf(fid,'110\t%g\n',Cv_a);
fprintf(fid,'111\t%g\n',Cv_b);
fclose(fid);
