N_STEP=2;
load RHO.dat
load U.dat
load P.dat
load X.dat
rho_p=RHO(N_STEP,:);
u_p  =U(N_STEP,:);
p_p  =P(N_STEP,:);
x_p  =X(N_STEP,:)-40;
figure(1)
hold on
plot(x_p,rho_p,'b.');
xlim([-10,90]);
ylim([0,3]);
xlabel('x');
ylabel('\rho');
figure(2)
hold on
plot(x_p,u_p,'b.');
xlim([-10,90]);
ylim([0,2.5]);
xlabel('x');
ylabel('u');
figure(3)
hold on
plot(x_p,p_p,'b.');
xlim([-10,90]);
ylim([0,5.5]);
xlabel('x');
ylabel('p');
