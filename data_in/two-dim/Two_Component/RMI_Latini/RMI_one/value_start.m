function value_start(column)
if nargin < 1
    column = 128
elseif strcmpi(column,{'INPUT'})
    column = input('Please input the column number: (Default 128 256 512) ')
else
    column
    if ischar(column)
        error("Not a string 'INPUT' was entered to represent the input!");
    elseif column ~= fix(column) || column <= 0
        error("Not a positive integer was entered to represent variable 'column'!")
    end
end

line=1683*column/128;
lambda=0.593333;
d_x=lambda/column
d_y=7.8/line
shock=round(0.1/d_y);


phi_1=0.0;
rho_2=1.351;
p_2=9.56;
v_2=0.0;
phi_2=0.0;
rho_3=5.494;
p_3=p_2;
v_3=v_2;
phi_3=1.0;

gamma=1.276;
M=1.21;
f=1/(2/(gamma+1)/M/M+(gamma-1)/(gamma+1));
g=2*gamma/(gamma+1)*M*M-(gamma-1)/(gamma+1);
rho_1=rho_2*f;
p_1=p_2*g;
v_1=(1-1/f)*(v_2+sqrt(gamma*p_2/rho_2)*M)+v_2/f;


CC = zeros(line,column);
PHI = zeros(line,column);
RHO = zeros(line,column);
U = zeros(line,column);
V = zeros(line,column);
P = zeros(line,column);
for j=1:column
for i=1:shock
       CC(i,j)  = phi_1;
       RHO(i,j) = rho_1;
       PHI(i,j) = phi_1;
       U(i,j)   = 0.0;
       V(i,j)   = v_1;
       P(i,j)   = p_1;      
end
end
for j=1:column
for i=(shock+1):(round(0.37/d_y)+1)
    CC(i,j) = index((j-0.5)*d_x,(i-0.5)*d_y);
end
end
for j=1:column
for i=(round(0.37/d_y)+2):line
    CC(i,j) = phi_3;
end
end
for j=1:column
for i=(shock+1):line
       RHO(i,j) = CC(i,j)*rho_3+(1.0-CC(i,j))*rho_2;
       PHI(i,j) = CC(i,j)*rho_3/RHO(i,j);
       U(i,j)   = 0.0;
       V(i,j)   = v_2;
       P(i,j)   = p_2;      
end
end


fid = fopen('PHI.dat','wt');
for j=1:line
fprintf(fid,'%10.8g\t',PHI(j,:));
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('Z_a.dat','wt');
for j=1:line
fprintf(fid,'%10.8g\t',CC(j,:));
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('RHO.dat','wt');
for j=1:line
fprintf(fid,'%10.8g\t',RHO(j,:));
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('U.dat','wt');
for j=1:line
fprintf(fid,'%10.8g\t',U(j,:));
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('V.dat','wt');
for j=1:line
fprintf(fid,'%10.8g\t',V(j,:));
fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('P.dat','wt');
for j=1:line
fprintf(fid,'%10.8g\t',P(j,:));
fprintf(fid,'\n');
end
fclose(fid);

eps=1e-9;
t_all=30;
step=5000000;
CFL=0.45;
gamma_a=1.093;

fid = fopen('config.dat','wt');
fprintf(fid,'1\t%g\n',t_all);
fprintf(fid,'2\t2\n');
fprintf(fid,'4\t%g\n',eps);
fprintf(fid,'5\t%i\n',step);
fprintf(fid,'6\t%g\n',gamma_a);
fprintf(fid,'7\t%g\n',CFL);
fprintf(fid,'10\t%g\n',d_x);
fprintf(fid,'11\t%g\n',d_y);
fprintf(fid,'13\t%g\n',column);
fprintf(fid,'14\t%g\n',line);
fprintf(fid,'17\t-7\n');
fprintf(fid,'18\t-4\n');
fprintf(fid,'106\t%g\n',gamma);
fclose(fid);
end
