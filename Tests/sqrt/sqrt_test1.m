% Test 6: Quadrature test: This program tests the convergence behaviour of r(FOM)2 in comparison to
% the standard Arnoldi approximation with quadrature as the number of
% quadrature points are varied.

%The code is tested on the square root function.

addpath(genpath('../'))


A = diag(1:1000);
n = size(A,1);

% choose right-hand side as normalized vector of all ones
b = ones(n,1); b = b/norm(b);

Aop = @(bx) A*bx;

p.n = size(A,1);
p.m = 20;
p.k = 15;
p.C = [];

quad_points = [1000,5000,10000,20000];
num_experiments = size(quad_points,2);

%Compute deflation subspace in advance
[U,~] = eigs(A,p.k,'smallestreal');
[p.U,~] = qr(U,0);
p.C = Aop(p.U);

p.f_scalar = @(zx) sqrt(zx);
p.f_matrix = @(Ax,bx) sqrtm(full(Ax))*bx;

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

e1 = zeros(p.m,1);
e1(1)=1;
exact = p.f_matrix(A,b);

%vectors to store results
arnoldi_err = zeros(1,num_experiments);
rfom2_err = zeros(1,num_experiments);
rfom2_err_g = zeros(1,num_experiments);

%Run Arnoldi to build basis V (of Krylov subspace) and Hessengerg matrix H
fprintf("\n Computing Arnoldi approximation fa...\n");
[V,H] = arnoldi(Aop, b, p);
  
%Call each of the methods using a new number of quadrature points for each call
for ix=1:num_experiments

    p.num_quad = quad_points(ix);
    
    %Standard Arnoldi Approximation
    fa = norm(b)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);

    %r(FOM)2 v1
    %Call all three implementations of rFOM2 for inverse square root.
    fr = rFOM2_v1(p,b,V,H);
    frg = rFOM2_v3(p,A,b,V,H);
    
    arnoldi_err(ix) = norm(fa - exact)/norm(exact);
    rfom2_err(ix) = norm(fr - exact)/norm(exact);
    rfom2_err_g(ix) = norm(frg - exact)/norm(exact);

end

%Plot Results
semilogy(arnoldi_err,'-v', 'LineWidth',1);
hold on; 
semilogy(rfom2_err,'-o', 'LineWidth',1);
hold on;
semilogy(rfom2_err_g,'-v', 'LineWidth',1);
hold off;

title('Square Root Test','interpreter','latex', 'FontSize', fontsize)
xlabel('problem index','interpreter','latex', 'FontSize', fontsize);
ylabel('Relative Error', 'FontSize',fontsize);
grid on;
lgd = legend('Arnoldi error', 'r(FOM)$^{2}$ v1','r(FOM)$^{2}$ v3','interpreter','latex');
set(lgd,'FontSize',fontsize);

