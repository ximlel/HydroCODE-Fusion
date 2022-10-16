function value_start(A)
if nargin < 1
    A = 8  %amplitude (number of grid)
elseif strcmpi(A,{'INPUT'})
    A = input('Please input the A number: (Default 50 25 20 10 5 2 0) ')
else
    A
    if ischar(A)
        error("Not a string 'INPUT' was entered to represent the input!");
    elseif A ~= fix(A) || A <= 0
        error("Not a positive integer was entered to represent variable 'A'!")
    end
end

column=200;
line=400;
K=1;  %wavenumber

gamma = 1.35;
Ma=1.5;
L_x=0.5/column;
L_y=1.0/line;

shock=round(line*0.1);
center=round(line*0.5)+A;  %center of the interface seperated two fluid.

u=0.0;
rho_uH=1;
rho_sH=2.7647;
p_uH=1;  %p_L and p_uH are same.
p_sH=4.4468;
v_uH=0;  %v_L and v_uH are same.
v_sH=1.4833;

rho_L=1.9;



index = @(y,x) y<A*cos(x*K/column*2*pi);

CC=zeros(line,column);
for j=1:(center-A)
CC(j,:)=ones(1,column);
end
for l=0:(K-1)
for i=1:round(column/K)
	y_1=A*cos((i-1)*K/column*2*pi);
    y_2=A*cos(i*K/column*2*pi);
for j=(center-A+1):(center+A)
    if(min(y_1,y_2)>(j-center))
       CC(j,i+l*round(column/K))=1; 
    elseif(max(y_1,y_2)<(j-center-1))
       CC(j,i+l*round(column/K))=0;
    else
       CC(j,i+l*round(column/K))=quad2d(index,j-center-1,j-center,i-1,i,'AbsTol',1e-3);
    end
end
end
end


fid = fopen('PHI.dat','wt');
for j=1:line
fprintf(fid,'%14.12g\t',CC(j,:));
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('RHO.dat','wt');
for j=1:shock
fprintf(fid,'%14.12g\t',rho_sH*ones(column,1));
fprintf(fid,'\n');
end
for j=(shock+1):line
fprintf(fid,'%14.12g\t',CC(j,:)*rho_uH+(1-CC(j,:))*rho_L);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('U.dat','wt');
for j=1:line
fprintf(fid,'%14.12g\t',u*ones(column,1));
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('V.dat','wt');
for j=1:shock
fprintf(fid,'%14.12g\t',v_sH*ones(column,1));
fprintf(fid,'\n');
end
for j=(shock+1):line
fprintf(fid,'%14.12g\t',v_uH*ones(column,1));
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('P.dat','wt');
for j=1:shock
fprintf(fid,'%14.12g\t',p_sH*ones(column,1));
fprintf(fid,'\n');
end
for j=(shock+1):line
fprintf(fid,'%14.12g\t',p_uH*ones(column,1));
fprintf(fid,'\n');
end
fclose(fid);
end
