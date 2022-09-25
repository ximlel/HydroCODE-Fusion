I = h5info('FLU_VAR.h5');
[TN, seq_T] = sort({I.Groups(:).Name});
N = size(seq_T,2);
rho_p = h5read('FLU_VAR.h5',[TN{2},'/RHO']);
u_p   = h5read('FLU_VAR.h5',[TN{2},'/U']);
e_p   = h5read('FLU_VAR.h5',[TN{2},'/E'])-0.5*u_p.^2;
x_p   = h5read('FLU_VAR.h5',[TN{2},'/X']);
time_plot = h5readatt('FLU_VAR.h5',TN{2},'time_plot');
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
