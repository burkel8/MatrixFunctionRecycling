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
problem = 'sqrt';


%% Parameters of solve
m = 30;  %Arnoldi cycle length
k = 20;  %recycle space dimension
N = 30;  %Parameter for Poisson and chemical potential matrix (value 
         %does not matter for other matrices)
num_quad = [5000];   %number of quadrature points (add as many differnt points to this list)

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
err_rFOM = zeros(1,num_tests);
err_rGMRES = zeros(1,num_tests);


[H,V] = arnoldi( A, b , n,m, 1);

for i=1:num_tests
fprintf('Running test number '); i
arnoldi_approx = norm(b)*V(:,1:m)*f_matrix(H(1:m,1:m),e1);
err_arnoldi(i) = norm(exact - arnoldi_approx);

quad_arnoldi_Approx = quad_arnoldi(b,V,H,m,num_quad(i),f_scalar);
err_quad_arnoldi(i) = norm(exact - quad_arnoldi_Approx);


[rFOM_approx] = rFOM2(b,V,H,m,k,U,C,num_quad(i), f_scalar);
err_rFOM(i) = norm(exact - rFOM_approx);

[rGMRES_approx] = rGMRES2(b,V,H,m,k,U,C,num_quad(i), f_scalar);
err_rGMRES(i) = norm(exact - rGMRES_approx);
end

x=1:1:num_tests;
points = num_quad(x);
semilogy(points,err_arnoldi,'-x');
hold on;
semilogy(points,err_quad_arnoldi,'-o');
hold on;
semilogy(points,err_rFOM,'-s');
title(' approximation accuracy vs number of quadrature nodes','interpreter','latex')
xlabel('number of quad nodes','interpreter','latex');
ylabel('$\| f(A)b - x_{m} \|_{fro}$','interpreter','latex');
grid on;
legend('Arnoldi','Arnoldi quad','rFOM');
xticks(points)

