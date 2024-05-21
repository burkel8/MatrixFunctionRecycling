% test_qcde_quad: Quadrature test: This program tests the convergence behaviour of r(FOM) in comparison to
% the standard Arnoldi approximation with quadrature as the number of quadrature points are varied.

%The code is tested on the inverse square root function.
addpath(genpath('../'))

set(0,...
 'defaultaxeslinewidth',1,...
'defaultaxesfontsize',18,...
'defaultlinelinewidth',2,...
'defaultpatchlinewidth',2,...
'defaultlinemarkersize',8,...
'defaulttextinterpreter','latex');

load(['data/4to4/periodic_L4_b3.55_k0.137n0_' num2str(1) '.mat'])
A = D;

Aop = @(bx) A'*(A*bx);

p.n = size(A,1);
p.m = 50;
p.k = 20;
b = rand(p.n,1);
Ab = A*b;
p.C = [];

quad_points = [1,2,3,4,5,10,20,30,50,100];
num_experiments = size(quad_points,2);

%Compute deflation subspace in advance
[U,~] = eigs(A'*A,p.k,'smallestabs');
[p.U,~] = qr(U,0);
p.C = Aop(p.U);

p.f_scalar = @(zx) 1./sqrt(zx);
p.f_matrix = @(Ax,bx) sqrtm(full(Ax))\bx;

%Paramters for fontsize and line width in plots
fontsize = 16;
linewidth = 4;


e1 = zeros(p.m,1);
e1(1)=1;
exact = p.f_matrix(A'*A,Ab);

%vectors to store results
arnoldi_err = zeros(1,num_experiments);
arnoldi_quad_err = zeros(1,num_experiments);
rfom2_err = zeros(1,num_experiments);
rfom2_err_closed = zeros(1,num_experiments);
rfom2_err_g = zeros(1,num_experiments);

%Run Arnoldi to build basis V (of Krylov subspace) and Hessengerg matrix H
fprintf("\n Computing Arnoldi approximation fa...\n");
[V,H] = arnoldi(Aop, Ab, p);
  


%Call each of the methods using a new number of quadrature points for each call
for ix=1:num_experiments

    p.num_quad = quad_points(ix);
    
    %Standard Arnoldi Approximation
    fa = norm(Ab)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);

    faa = quad_arnoldi_invSqrt(p,Ab,V,H);

    %r(FOM)2 v1
    %Call all three implementations of rFOM2 for inverse square root.
    fr = rFOM2_v1_invSqrt(p,Ab,V,H);
    frc = rFOM2_v1_invSqrt_closed(p,Ab,V,H);
    frg = rFOM2_v3_invSqrt(p,Ab,V,H);
    
    arnoldi_err(ix) = norm(fa - exact)/norm(exact);
    arnoldi_quad_err(ix) = norm(faa - exact)/norm(exact);
    rfom2_err(ix) = norm(fr - exact)/norm(exact);
    rfom2_err_closed(ix) = norm(frc - exact)/norm(exact);
    rfom2_err_g(ix) = norm(frg - exact)/norm(exact);

end

%Plot Results
semilogy(arnoldi_err,'--');
hold on; 
semilogy(arnoldi_quad_err,'-s');
hold on;
semilogy(rfom2_err,'--s');
hold on;
semilogy(rfom2_err_closed,':v');
hold on;
semilogy(rfom2_err_g,'-.o');
hold on;
xticks([1,2,3,4,5,10,20,30,50,100])
xticklabels({'1','2','3','4', '5', '10', '20', '30', '50', ' 100'});
xlabel('Number of quadrature nodes','interpreter','latex');
ylabel('Relative error');
grid on;
legend('FOM','FOM quad','rFOM v1', 'rFOM v2','rFOM v3','interpreter','latex');



