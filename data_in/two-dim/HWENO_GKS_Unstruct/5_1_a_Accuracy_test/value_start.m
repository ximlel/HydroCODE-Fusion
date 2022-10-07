function value_start(column)
if nargin < 1
    column = 20
elseif strcmpi(column,{'INPUT'})
    column = input('Please input the column number: (Default 20 40 80 160 320 640 1280) ')
else
    column
    if ischar(column)
        error("Not a string 'INPUT' was entered to represent the input!");
    elseif column ~= fix(column) || column <= 0
        error("Not a positive integer was entered to represent variable 'column'!")
    end
end

line=3;
d_x=2.0/column;

u0=1
p0=1
phi0=1

rho=zeros(column,1);
f=@(x) 1+0.2*sin(pi*x);

for i=1:column
    rho(i)=quad(f,(i-1)*d_x,i*d_x,1e-18)/d_x;
end
fid = fopen('RHO.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
for i=1:column
    u(i)=u0;
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
for i=1:column
    p(i)=p0;
end
fid = fopen('P.dat','wt');
for j=1:line
fprintf(fid,'%12.10g\t',p);
fprintf(fid,'\n');
end
fclose(fid);


eps=1e-9;
t_all=10;
step=500000;
gamma=1.4;
CFL=0.45;

fid = fopen('config.dat','wt');
fprintf(fid,'1\t%g\n',t_all);
fprintf(fid,'4\t%g\n',eps);
fprintf(fid,'5\t%i\n',step);
fprintf(fid,'6\t%g\n',gamma);
fprintf(fid,'7\t%g\n',CFL);
fprintf(fid,'10\t%g\n',d_x);
fprintf(fid,'11\t%g\n',d_x);
fprintf(fid,'13\t%i\n',column);
fprintf(fid,'14\t%i\n',line);
fprintf(fid,'17\t-7\n');
fprintf(fid,'18\t-7\n');
fclose(fid);
end
