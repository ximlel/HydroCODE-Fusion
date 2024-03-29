Path="u_0/"
load(fullfile(Path,"RHO.dat"));
load(fullfile(Path,"U.dat"));
load(fullfile(Path,"V.dat"));
load(fullfile(Path,"P.dat"));
%U=U(end:-1:1,:);
%V=V(end:-1:1,:);
line = 50;
xa=linspace(0.1,9.9,line);
ya=xa;
[x,y]=meshgrid(xa,ya);

S = P./(RHO.^1.4);
figure(1);
mesh(x,y,RHO);
figure(2);
mesh(x,y,S);
figure(3);
quiver(x,y,U,V);
