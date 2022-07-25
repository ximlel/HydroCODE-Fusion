N_STEP=48;
load RHO.txt
load U.txt
load X.txt
rho_p=RHO(N_STEP,:);
u_p  =U(N_STEP,:);
x_p  =X(N_STEP,:);
figure(1)
plot(x_p,rho_p,'.');
ylim([0,1]);
xlabel('x');
ylabel('\rho')
figure(2)
plot(x_p,u_p,'.');
ylim([0,1]);
xlabel('x');
ylabel('u')
