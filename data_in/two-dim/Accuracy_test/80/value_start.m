line=3;
column=80
d_x=2.0/column;

u0=1
p0=1
phi0=1

rho=zeros(column,1);
f=@(x) 1+0.2*sin(pi*x);

for i=1:column
    rho(i)=quad(f,(i-1)*d_x,i*d_x,1e-18,0)/d_x;
end
fid = fopen('RHO.dat','wt');
for j=1:line
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
for i=1:column
    u(i)=u0;
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
for i=1:column
    p(i)=p0;
end
fid = fopen('P.dat','wt');
for j=1:line
fprintf(fid,'%12.10f\t',p);
fprintf(fid,'\n');
end
fclose(fid);
