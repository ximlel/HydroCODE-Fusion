N_STEP=2;
load RHO.dat
load X.dat
rho_p=RHO(N_STEP,:);
x_p  =X(N_STEP,:)-50;
figure(1)
hold on
plot(x_p,rho_p,'b.');
xlim([-10,400]);
ylim([0,3.5]);
xlabel('x');
ylabel('\rho');
