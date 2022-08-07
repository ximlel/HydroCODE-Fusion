N_STEP=2;
load RHO.dat
load U.dat
load X.dat
rho_p=RHO(N_STEP,:);
u_p  =U(N_STEP,:);
x_p  =X(N_STEP,:);
figure(1)
hold on
plot(x_p,rho_p,'.');
ylim([0,1]);
xlabel('x');
ylabel('\rho')
figure(2)
hold on
plot(x_p,u_p,'.');
ylim([0,1]);
xlabel('x');
ylabel('u')
