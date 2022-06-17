% n = 1000;
% A = -spdiags((1:n)',0,n,n);
% [n,~]=size(A);


% N = 40;
% D2 = (N+1)^2*gallery('tridiag',N);
% I = speye(N);
% D2 = kron(I,D2) + kron(D2,I);
% o = ones(N,1);
% D1 = (N+1)/2*spdiags([-o,0*o,o],-1:1,N,N);
% D1 = kron(I,D1) + kron(D1,I);
% A = D2 + 0*D1; 
% [n,~]=size(A);

% N = 50;
% e = ones(N,1);
% A = (N+1)^2*gallery('poisson',N);
% n = N*N;
% s = eigs(A,1,'smallestabs');

% load("smallLQCD_A1.mat");
% A = A1;
% [n,~] = size(A);

% load('Schw_16.mat');
% A = S;
% [n,~] = size(A);


%% Define rhs
b = rand(n,1);
b = b/norm(b);

%% define exact
exact = A\b;

%% Parameters of solve
m = 20;
k = 15;

num_quad = [3000,5000,8000];
num_tests = size(num_quad,2);


e1 = zeros(m,1); e1(1)=1;

%% Run Arnoldi
[H,V] = arnoldi( A, b , n,m, 0);

%% Create deflation subspace
%As = A + 0.001*rand(n,n);
%[Hs,Vs] = arnoldi( A, b , n, 1600, 0);

%[P,D] = eigs(Hs(1:m,1:m),k,'smallestabs');
%U = Vs(:,1:m)*P;
%norm(A*U - U*D)/norm(U)

[U,~]=eigs(A,k,'smallestabs');
C=A*U;

err_arnoldi = zeros(1,num_tests);
err_quad_arnoldi = zeros(1,num_tests);
err_rFOM = zeros(1,num_tests);
err_rGMRES = zeros(1,num_tests);


for i=1:num_tests
fprintf('Running test number '); i
arnoldi_approx = norm(b)*V(:,1:m)*(H(1:m,1:m)\e1);
err_arnoldi(i) = norm(exact - arnoldi_approx);

quad_arnoldi_Approx = quad_arnoldi(b,V,H,m,num_quad(i));
err_quad_arnoldi(i) = norm(exact - quad_arnoldi_Approx);


[rFOM_approx] = rFOM2(b,V,H,m,k,U,C,num_quad(i));
err_rFOM(i) = norm(exact - rFOM_approx);

[rGMRES_approx] = rGMRES2(b,V,H,m,k,U,C,num_quad(i));
err_rGMRES(i) = norm(exact - rGMRES_approx);
end

x=1:1:num_tests;
points = num_quad(x);
semilogy(points,err_arnoldi,'-x');
hold on;
semilogy(points,err_quad_arnoldi,'-o');
hold on;
semilogy(points,err_rFOM,'-s');
hold on;
semilogy(points, err_rGMRES);
title(' approximation accuracy vs number of quadrature nodes','interpreter','latex')
xlabel('number of quad nodes','interpreter','latex');
ylabel('$\| f(A)b - x_{m} \|_{fro}$','interpreter','latex');
grid on;
legend('Arnoldi','Arnoldi quad','rFOM','rGMRES');
xticks(points)

