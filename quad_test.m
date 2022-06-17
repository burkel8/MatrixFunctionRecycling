load("smallLQCD_A1.mat");
A = A1;
[n,~] = size(A);
%A = A - 0.6*speye(n);

%% Define rhs
b = rand(n,1);
b = b/norm(b);


%% Parameters of solve
m = 40;
k = 20;

num_quad = [100,200,300,400,500,600];
num_tests = size(num_quad,2);

problem = 'inverse';

f1_scalar = @(zx) 1.0/zx;
f1_matrix = @(Ax,bx) Ax\bx;

f2_scalar = @(zx) 1./sqrt(zx);
f2_matrix = @(Ax,bx) sqrtm(full(Ax))\bx;

if problem =="inverse"
    f_scalar = @(zx) f1_scalar(zx);
    f_matrix = @(Ax,bx) f1_matrix(Ax,bx);
elseif problem == "sign"
    f_scalar = @(zx) f2_scalar(zx);
    f_matrix = @(Ax,bx) f2_matrix(Ax,bx);
else
    error("ERROR : unknown function chosen!\n");
end

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
hold on;
semilogy(points, err_rGMRES,'-s');
title(' approximation accuracy vs number of quadrature nodes','interpreter','latex')
xlabel('number of quad nodes','interpreter','latex');
ylabel('$\| f(A)b - x_{m} \|_{fro}$','interpreter','latex');
grid on;
legend('Arnoldi','Arnoldi quad','rFOM','rGMRES');
xticks(points)

