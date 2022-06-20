load("err_recyclepace_smallLQCDsignm50k20numquad1000.mat")
load("err_recyclepace_smallLQCDsignm50k20numquad1000matchanging.mat")

num_systems = 15;
xx = 1:num_systems+1;


semilogy(xx,rhschanging_err_recycle_space);
hold on;
semilogy(xx,matchanging_err_recycle_space);
hold off;
grid on;
title("Accuracy of $\mathcal{U}$ as an eigenvector approximation", 'Interpreter', 'latex');
legend('QCD sign - rhs changing only', 'QCD sign - Matrix and rhs changing $\epsilon = 0.01$','interpreter','latex');
xlabel('system index $i$', 'interpreter', 'latex');
ylabel('fdf','interpreter','latex');