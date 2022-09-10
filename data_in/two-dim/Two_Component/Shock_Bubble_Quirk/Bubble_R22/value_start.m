line=890;
column=2500;
shock=125;
delta_x=0.89/line;
R=0.25/delta_x;
center_x=0.5/delta_x;
center_y=0.89/2/delta_x;
idx = @(y,x) sqrt((x-center_x).^2+(y-center_y).^2)>R;

phi_1=1.0;
rho_2=1.0;
u_2=0.0;
p_2=1.0;
phi_2=1.0;
rho_3=rho_2*0.4*0.72/0.249/0.365;
%rho_3=rho_2*0.287/0.091;
u_3=u_2;
p_3=p_2;
phi_3=0;

gamma=1.4;
M=1.22;
f=1/(2/(gamma+1)/M/M+(gamma-1)/(gamma+1));
g=2*gamma/(gamma+1)*M*M-(gamma-1)/(gamma+1);
rho_1=rho_2*f;
u_1=(1-1/f)*(u_2+sqrt(gamma*p_2/rho_2)*M)+u_2/f;
p_1=p_2*g;


CC = zeros(line,column);
PHI = zeros(line,column);
RHO = zeros(line,column);
U = zeros(line,column);
V = zeros(line,column);
P = zeros(line,column);
for j=1:line
for i=1:shock
       CC(j,i)  = 1.0;
       RHO(j,i) = rho_1;
       PHI(j,i) = 1.0;
       U(j,i)   = u_1;
       V(j,i)   = 0.0;
       P(j,i)   = p_1;      
end
end
for j=1:line
for i=(shock+1):column
    if(i<center_x+0.5)
       far_p(1)  = center_x-i+1;
       near_p(1) = center_x-i; 
    elseif(i>center_x+0.5)
       far_p(1)  = i-center_x;
       near_p(1) = i-center_x-1;
    end
    if(j<center_y+0.5)
       far_p(2)  = center_y-j+1;
       near_p(2) = center_y-j;
    elseif(j>center_y+0.5)
       far_p(2)  = j-center_y;
       near_p(2) = j-center_y-1;
    end
    if(norm(near_p)<R && norm(far_p)>R)
       CC(j,i) = quad2d(idx,j-1,j,i-1,i,'AbsTol',1e-3);
    elseif(norm(near_p)>=R)
       CC(j,i) = 1.0;
    else
       CC(j,i) = 0.0;
    end
end
end
for j=1:line
for i=(shock+1):column
       RHO(j,i) = CC(j,i)*rho_2+(1.0-CC(j,i))*rho_3;
       PHI(j,i) = CC(j,i)*rho_2/RHO(j,i);
       U(j,i)   = u_2;
       V(j,i)   = 0.0;
       P(j,i)   = p_2;      
end
end
RHO=RHO(:,end:-1:1);
PHI=PHI(:,end:-1:1);
U=-1.0*U(:,end:-1:1);
V=V(:,end:-1:1);
P=P(:,end:-1:1);

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
