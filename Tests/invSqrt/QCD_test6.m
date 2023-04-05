% Test 6: Quadrature test: This program tests the convergence behaviour of r(FOM)2 in comparison to
% the standard Arnoldi approximation with quadrature as the number of
% quadrature points are varied.

%The code is tested on the inverse square root function.
addpath(genpath('../'))
load('../../data/smallLQCD.mat');
A = A1;

Aop = @(bx) A*bx;

p.n = size(A,1);
p.m = 40;
p.k = 20;
b = rand(p.n,1);
p.C = [];

quad_points = [1,2,3,4,5,10,20,30,50,100];
num_experiments = size(quad_points,2);

%Compute deflation subspace in advance
[U,~] = eigs(A,p.k,'smallestabs');
[p.U,~] = qr(U,0);
p.C = Aop(p.U);

p.f_scalar = @(zx) 1./sqrt(zx);
p.f_matrix = @(Ax,bx) sqrtm(full(Ax))\bx;

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

e1 = zeros(p.m,1);
e1(1)=1;
exact = p.f_matrix(A,b);

%vectors to store results
arnoldi_err = zeros(1,num_experiments);
rfom2_err = zeros(1,num_experiments);
rfom2_err_closed = zeros(1,num_experiments);
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
    fr = rFOM2_v1_invSqrt(p,b,V,H);
    frc = rFOM2_v1_invSqrt_closed(p,b,V,H);
    frg = rFOM2_v3_invSqrt(p,b,V,H);
    
    arnoldi_err(ix) = norm(fa - exact)/norm(exact);
    rfom2_err(ix) = norm(fr - exact)/norm(exact);
    rfom2_err_closed(ix) = norm(frc - exact)/norm(exact);
    rfom2_err_g(ix) = norm(frg - exact)/norm(exact);

end

%Plot Results
semilogy(arnoldi_err,'-v', 'LineWidth',1);
hold on; 
semilogy(rfom2_err,'-o', 'LineWidth',1);
hold on;
semilogy(rfom2_err_closed,'-v', 'LineWidth',1);
hold on;
semilogy(rfom2_err_g,'-v', 'LineWidth',1);
hold on;

title('Quadrature approximation Test for $f(z) = z^{-0.5}$','interpreter','latex', 'FontSize', fontsize)
xlabel('Quad Points','interpreter','latex', 'FontSize', fontsize);
ylabel('Relative Error', 'FontSize',fontsize);
grid on;
lgd = legend('Arnoldi error', 'r(FOM)$^{2}$ v1', 'r(FOM)$^{2}$ v2','r(FOM)$^{4}$ v2','interpreter','latex');
set(lgd,'FontSize',fontsize);


