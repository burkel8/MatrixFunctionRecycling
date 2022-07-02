%% Step 1: Choose parameters for program

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
problem = 'inverse';


%% Parameters of solve
m = 40;  %Arnoldi cycle length
k = 20;  %recycle space dimension
N = 80;  %Parameter for Poisson and chemical potential matrix (value 
         %does not matter for other matrices)
num_quad = [100,200,300,400,500,600];   %number of quadrature points (add as many differnt points to this list)

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;
%%%%%%%%%%%%%%    END USER INPUT HERE  %%%%%%%%%%%%%%%%%%%

[A,n] = return_matrix(which_matrix,N);
[f_scalar, f_matrix] = return_function(problem);

%% Define rhs
b = rand(n,1);
b = b/norm(b);

num_tests = size(num_quad,2);

%define exact
exact = f_matrix(A,b);
e1 = zeros(m,1); e1(1)=1;

%% Run Arnoldi
%Create a deflation subspace
[U,~] = eigs(A,k,'smallestabs'); % 'Tolerance',10e-09);
 C = A*U;



%% Run Arnoldi on a close matrix

err_arnoldi = zeros(1,num_tests);
err_quad_arnoldi = zeros(1,num_tests);
err_rFOM_v1 = zeros(1,num_tests);
err_rFOM_v2 = zeros(1,num_tests);


[H,V] = arnoldi( A, b , n,m, 1);

for i=1:num_tests
fprintf('Running test number '); i
arnoldi_approx = norm(b)*V(:,1:m)*f_matrix(H(1:m,1:m),e1);
err_arnoldi(i) = norm(exact - arnoldi_approx);

quad_arnoldi_Approx = quad_arnoldi(b,V,H,m,num_quad(i),f_scalar);
err_quad_arnoldi(i) = norm(exact - quad_arnoldi_Approx);

[rFOM_v1_approx] = rFOM2_v1(b,V,H,m,k,U,C,num_quad(i), f_scalar);
err_rFOM_v1(i) = norm(exact - rFOM_v1_approx);


[rFOM_v2_approx] = rFOM2_v2(b,V,H,m,k,U,C,num_quad(i), f_scalar);
err_rFOM_v2(i) = norm(exact - rFOM_v2_approx);

end

x=1:1:num_tests;
points = num_quad(x);
semilogy(points,err_arnoldi,'-x','LineWidth',1);
hold on;
semilogy(points,err_quad_arnoldi,'-o','LineWidth',1);
hold on;
semilogy(points,err_rFOM_v1,'-x','LineWidth',1);
hold on;
semilogy(points,err_rFOM_v2,'-s','LineWidth',1);
hold off;
title(' approximation accuracy vs number of quadrature nodes','interpreter','latex','FontSize',fontsize)
xlabel('number of quad nodes','interpreter','latex','FontSize',fontsize);
ylabel('$\| f(\textbf{A})\textbf{b} - \textbf{x}_{m} \|_{2}$','interpreter','latex','FontSize',fontsize);
grid on;
lgd = legend('Arnoldi','Arnoldi quad','rFOM$^{2}$ (v1)','rFOM$^{2}$ (v2)','interpreter','latex');
set(lgd,'FontSize',fontsize);
xticks(points)

