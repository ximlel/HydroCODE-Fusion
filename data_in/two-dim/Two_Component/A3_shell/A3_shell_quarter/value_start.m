% dimension:g,cm,ms
% gamma_1  = 1.4;
% gamma_23 = 3.0;
% gamma_4  = 1.4;
line = 600;
delta_x = 15/line
R_1 = 10/delta_x;
R_2 = 10.2/delta_x;
R_3 = 10.25/delta_x;
R_4 = 10.5/delta_x;

p   = 1.01325;%p_123
p_4 = 0.101325;
rho_1  = 0.00129;
rho_23 = 19.237;
rho_4  = 0.000129;
% u_1 = 0.0;
% u_2 = 0.0;
% u_3 = 200.0;
% u_4 = 0.0;

CC = zeros(line,line);
PHI = zeros(line,line);
RHO = zeros(line,line);
U = zeros(line,line);
V = zeros(line,line);
P = zeros(line,line);

idx_1 = @(y,x) sqrt(x.^2+y.^2)<R_1;
idx_3 = @(y,x) sqrt(x.^2+y.^2)>R_3;
u_fun = @(y,x) -x./sqrt(x.^2+y.^2)*200.0;
v_fun = @(y,x)  y./sqrt(x.^2+y.^2)*200.0;
u_full_fun = @(y,x) -x./sqrt(x.^2+y.^2)*200.0.*(sqrt(x.^2+y.^2)>R_2&sqrt(x.^2+y.^2)<R_3);
v_full_fun = @(y,x)  y./sqrt(x.^2+y.^2)*200.0.*(sqrt(x.^2+y.^2)>R_2&sqrt(x.^2+y.^2)<R_3);

far_p = zeros(1,2);
near_p = zeros(1,2);

for j=1:line
for i=1:line
       far_p(1) = i;
       far_p(2) = j;
       near_p(1) = i-1;
       near_p(2) = j-1;
       
    if(norm(far_p)<R_1)
       CC(j,i) = 1;
       PHI(j,i) = 1;
       RHO(j,i) = rho_1;
       U(j,i) = 0;
       V(j,i) = 0;
       P(j,i) = p;
    elseif(norm(near_p)>R_1 && norm(far_p)<R_2)
       CC(j,i) = 0;
       PHI(j,i) = 0;
       RHO(j,i) = rho_23;   
       U(j,i) = 0;
       V(j,i) = 0;
       P(j,i) = p;
    elseif(norm(near_p)>R_2 && norm(far_p)<R_3)
       CC(j,i) = 0;
       PHI(j,i) = 0;
       RHO(j,i) = rho_23;
       U(j,i) = quad2d(u_fun,-j,-j+1,i-1,i,'AbsTol',1e-8);
       V(j,i) = quad2d(v_fun,-j,-j+1,i-1,i,'AbsTol',1e-8);
       P(j,i) = p;
    elseif(norm(near_p)>R_3)
       CC(j,i) = 1;
       PHI(j,i) = 1;
       RHO(j,i) = rho_4;
       U(j,i) = 0;
       V(j,i) = 0;
       P(j,i) = p_4;
    elseif(norm(near_p)<=R_1)
       CC(j,i) = quad2d(idx_1,j-1,j,i-1,i,'AbsTol',1e-3);
       RHO(j,i) = CC(j,i)*rho_1+(1.0-CC(j,i))*rho_23;
       PHI(j,i) = CC(j,i)*rho_1/RHO(j,i);
       U(j,i) = 0;
       V(j,i) = 0;
       P(j,i) = p;
    elseif(norm(near_p)<=R_2)
       CC(j,i) = 0;
       PHI(j,i) = 0;
       RHO(j,i) = rho_23;
       U(j,i) = quad2d(u_full_fun,-j,-j+1,i-1,i,'AbsTol',1e-1);
       V(j,i) = quad2d(v_full_fun,-j,-j+1,i-1,i,'AbsTol',1e-1);
       P(j,i) = p;
    else
       CC(j,i) = quad2d(idx_3,j-1,j,i-1,i,'AbsTol',1e-3);
       RHO(j,i) = CC(j,i)*rho_4+(1.0-CC(j,i))*rho_23;
       PHI(j,i) = CC(j,i)*rho_4/RHO(j,i);
       U(j,i) = quad2d(u_full_fun,-j,-j+1,i-1,i,'AbsTol',1e-1)*rho_23/RHO(j,i);
       V(j,i) = quad2d(v_full_fun,-j,-j+1,i-1,i,'AbsTol',1e-1)*rho_23/RHO(j,i);
       P(j,i) = CC(j,i)*p_4+(1.0-CC(j,i))*p;      
    end
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
