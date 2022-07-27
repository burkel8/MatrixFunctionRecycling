%This script runs an experiment which plots the number of matvecs required
%to reach a specified accuracy.

%%First choose the matrix possible options are
%-- A small lattice QCD matrix of size 3072x3072 ("smallLQCD")
%-- A poisson matrix of size N*N x N*N (user specifies N) ("poisson")
%-- A chemical potential matrix of size N*N x N*N (user specifies N) ("chemical_potantial")
which_matrix = "smallLQCD";   

%%Choose the function . Available options are
% -- inverse function ("inverse")
% -- Sign function ("sign")
% -- log function ("log")
% -- square root function ("sqrt")
problem = 'sign';



k = 20;  %recycle space dimension
N = 100;  %Parameter for Poisson and chemical potential matrix (value 
         %does not matter for other matrices)
num_quad = 300;   %number of quadrature points (add as many differnt points to this list)
tol = 0.000001;
%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;
%%%%%%%%%%%%%%    END USER INPUT HERE  %%%%%%%%%%%%%%%%%%%

[A,n] = return_matrix(which_matrix,N);
[f_scalar, f_matrix] = return_function(problem);

% Define rhs
b = rand(n,1);
b = b/norm(b);

%compute exact solution
exact = f_matrix(A,b);


%Run Arnoldi
%Create a deflation subspace
[U,~] = eigs(A,k,'smallestabs'); % 'Tolerance',10e-09);
 C = A*U;

%Set all the below variables to a dummy value.
err_arnoldi_result = 60;
err_rFOM_v1_result = 60;
err_rFOM_v2_result = 60;
err_rFOM_v3_result = 60;

i = 1; %first Arnoldi cycle length

%keep running each method of approximation with a new cycle
%length until we converge
while err_arnoldi_result > tol

[H,V] = arnoldi( A, b , n,i, 1);
e1 = zeros(i,1); e1(1)=1;
arnoldi_approx = norm(b)*V(:,1:i)*f_matrix(H(1:i,1:i),e1);
err_arnoldi_result = norm(exact - arnoldi_approx)
err_arnoldi(i) = err_arnoldi_result;

if err_rFOM_v1_result > tol
[rFOM_v1_approx] = rFOM2_v1(b,V,H,i,k,U,C,num_quad, f_scalar);
err_rFOM_v1_result = norm(exact - rFOM_v1_approx)
err_rFOM_v1(i) = err_rFOM_v1_result;
end

if err_rFOM_v2_result > tol
[rFOM_v2_approx] = rFOM2_v2(b,V,H,i,k,U,C,num_quad, f_scalar);
err_rFOM_v2_result = norm(exact - rFOM_v2_approx)
err_rFOM_v2(i) = err_rFOM_v2_result;
end

if err_rFOM_v3_result > tol
[rFOM_v3_approx] = rFOM2_v3(b,V,H,i,k,U,C,num_quad, f_scalar, f_matrix);
err_rFOM_v3_result = norm(exact - rFOM_v3_approx)
err_rFOM_v3(i) = err_rFOM_v3_result;
end

i = i + 1;
end

% plot results
semilogy(err_arnoldi,'-s','LineWidth',1);
hold on;
semilogy(err_rFOM_v1,'-v','LineWidth',1 ,'MarkerSize', 9);
hold on;
semilogy(err_rFOM_v2,'-o','LineWidth',1);
hold on;
semilogy(err_rFOM_v3,'-s','LineWidth',1);
hold off;
title(' Reduction of MAT-VECs in rFOM$^{2}$  ','interpreter','latex','FontSize',fontsize)
xlabel('number of MAT-VEC operations','interpreter','latex','FontSize',fontsize);
ylabel('$\| f(\textbf{A})\textbf{b} - \tilde{f}_{i} \|_{2}$','interpreter','latex','FontSize',fontsize);
grid on;
lgd = legend('Arnoldi','rFOM$^{2}$ $\tilde{f}_{1}$','rFOM$^{2}$ $\tilde{f}_{2}$', 'rFOM$^{2}$ $\tilde{f}_{3}$','interpreter','latex');
set(lgd,'FontSize',fontsize);



