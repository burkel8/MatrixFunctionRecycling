%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;


load("eps0.mat");
eps0 = eigs_monitor;
load("eps0001.mat")
eps0001 = eigs_monitor;
load("eps001.mat")
eps001 = eigs_monitor;
load("eps01.mat");
eps01 = eigs_monitor;
load("eps1.mat");
eps1 = eigs_monitor;

num_systems = 10;
xx = 1:num_systems;


semilogy(xx,eps0,'LineWidth',1);
hold on;
semilogy(xx,eps0001,'LineWidth',1);
hold on;
semilogy(xx,eps001,'LineWidth',1);
hold on;
semilogy(xx,eps01,'LineWidth',1);
hold on;
semilogy(xx,eps1,'LineWidth',1);
hold off;
grid on;
title("QCD sign function - Accuracy of $\mathcal{U}$ as an eigenvector approximation", 'Interpreter', 'latex','FontSize',fontsize);
lgd = legend('$\epsilon=0$', '$\epsilon = 0.001$',"$\epsilon = 0.01$", "$\epsilon=0.1$", "$\epsilon = 1$",'interpreter','latex');
set(lgd,'FontSize',fontsize);