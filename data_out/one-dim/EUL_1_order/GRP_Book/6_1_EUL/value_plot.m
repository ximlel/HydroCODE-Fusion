I = h5info('FLU_VAR.h5');
[TN, seq_T] = sort({I.Groups(:).Name});
N = size(seq_T,2);
rho_p = h5read('FLU_VAR.h5',[TN{2},'/RHO']);
u_p   = h5read('FLU_VAR.h5',[TN{2},'/U']);
e_p   = h5read('FLU_VAR.h5',[TN{2},'/E'])-0.5*u_p.^2;
x_p   = h5read('FLU_VAR.h5',[TN{2},'/X']);
time_plot = h5readatt('FLU_VAR.h5',TN{2},'time_plot');
rho_p=RHO(N_STEP,:);
u_p  =U(N_STEP,:);
x_p  =X(N_STEP,:);
figure(1)
hold on
plot(x_p,rho_p,'r.');
ylim([0,1]);
xlabel('x');
ylabel('\rho');
figure(2)
hold on
plot(x_p,u_p,'r.');
ylim([0,1]);
xlabel('x');
ylabel('u');
