N_STEP=2;
load RHO.dat
load U.dat
load E.dat
load X.dat
rho_p=RHO(N_STEP,:);
u_p  =U(N_STEP,:);
e_p  =E(N_STEP,:)-0.5*U(N_STEP,:).^2;
x_p  =X(N_STEP,:);
figure(1)
hold on
plot(x_p,rho_p,'rx');
ylim([0,1]);
xlabel('x');
ylabel('\rho');
figure(2)
hold on
plot(x_p,u_p,'r.');
ylim([0,1]);
xlabel('x');
ylabel('u');
figure(4)
hold on
plot(x_p,e_p,'r.');
ylim([1.7,3.3]);
xlabel('x');
ylabel('e');
