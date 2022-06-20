load("eps0.mat");
eps0 = err_recycle_space;
load("eps0001.mat")
eps0001 = err_recycle_space;
load("eps001.mat")
eps001 = err_recycle_space;
load("eps01.mat");
eps01 = err_recycle_space;
load("eps1.mat");
eps1 = err_recycle_space;

num_systems = 10;
xx = 1:num_systems;


semilogy(xx,eps0);
hold on;
semilogy(xx,eps0001);
hold on;
semilogy(xx,eps001);
hold on;
semilogy(xx,eps01);
hold on;
semilogy(xx,eps1);
hold off;
grid on;
title("QCD sign function - Accuracy of $\mathcal{U}$ as an eigenvector approximation", 'Interpreter', 'latex');
legend('$\epsilon=0$', '$\epsilon = 0.001$',"$\epsilon = 0.01$", "$\epsilon=0.1$", "$\epsilon = 1$",'interpreter','latex');
xlabel('problem index $i$', 'interpreter', 'latex');
ylabel('$ \frac{\|\textbf{U} \textbf{D} - \textbf{A} \textbf{U}\|_{fro}}{\|\textbf{U}\|_{fro}} $','interpreter','latex');